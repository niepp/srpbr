#include <windows.h>
#include <tchar.h>
#include <vector>
#include <iostream>
#include <cassert>

#include "math3d.h"
#include "texture.h"
#include "model.h"

const vector3_t cOne(1.0f, 1.0f, 1.0f);

/*
 有一边是水平的三角形
 l----------r
  \        /
   \      /
    \    /
	 \  /
	  \/
	   p
*/
struct scan_tri_t
{
	interp_vertex_t p;
	interp_vertex_t l, r;
};

enum class shading_model_t
{
	cSM_Color = 0,
	cSM_Wireframe,
	cSM_Phong,
	cSM_PBR,
	cSM_MAX,
};

bool has_spec = true;

float reci_freq;
int64_t tick_start;
uint32_t frame_count = 0, frame_rate = 0;
HWND hwnd = 0;
HDC screenDC;

int width = 0, height = 0;
float angle_speed = 1.0f;
float light_angle = cPI;
float eyedist = 3.5f;

struct uniformbuffer_t
{
	vector3_t eye;
	matrix_t world;
	matrix_t view;
	matrix_t proj;
	matrix_t mvp;
	vector3_t light_dir;
	vector3_t light_intensity;
	float specular_power;
};

uniformbuffer_t uniformbuffer;
uint32_t *framebuffer = nullptr;
float *zbuffer = nullptr;

texture2d_t albedo_tex;
texture2d_t metallic_tex;
texture2d_t roughness_tex;
texture2d_t normal_tex;

shading_model_t shading_model = shading_model_t::cSM_PBR;

model_t sphere_model;

template <typename T, int n>
int array_size(T(&)[n])
{
	return n;
}

void write_pixel(int x, int y, uint32_t color)
{
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);
	framebuffer[y * width + x] = color;
}

void write_depth(int x, int y, float z)
{
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);
	zbuffer[y * width + x] = z;
}

bool depth_test(int x, int y, float z)
{
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);
	float nowz = zbuffer[y * width + x];
	return z >= nowz;
}

float gamma(float linear_color)
{
	const float cInvg = 1.0f / 2.2f;
	return pow(linear_color, cInvg);
}

// hdr tone mapping
float aces(float value)
{
	float a = 2.51f;
	float b = 0.03f;
	float c = 2.43f;
	float d = 0.59f;
	float e = 0.14f;
	value = (value * (a * value + b)) / (value * (c * value + d) + e);
	return clamp(value, 0.0f, 1.0f);
}

void reinhard_mapping(vector3_t& color)
{
	for (int i = 0; i < 3; ++i)
	{
		color.vec[i] = aces(color.vec[i]);
		color.vec[i] = gamma(color.vec[i]);
	}
}


void phong_shading(const interp_vertex_t& p, vector4_t& out_color)
{
	texcoord_t uv = p.uv;
	uv.u /= p.pos.w;
	uv.v /= p.pos.w;

	vector3_t wnor = p.nor / p.pos.w;
	wnor.normalize();

	vector3_t wpos = p.wpos / p.pos.w;
	vector3_t l = -uniformbuffer.light_dir;

	vector4_t albedo = albedo_tex.sample(uv);

	float NoL = dot(wnor, l);
	NoL = max(NoL, 0.15f);

	vector3_t diffuse = albedo.to_vec3() * uniformbuffer.light_intensity * NoL;

	vector3_t v = uniformbuffer.eye - wpos;
	v.normalize();

	// (v + l) / 2
	vector3_t h = v + l;
	h.normalize();

	float NoH = max(dot(wnor, h), 0.0f);
	vector3_t ks(0.5f, 0.5f, 0.5f);
	vector3_t specular = ks * uniformbuffer.light_intensity * max(0.0f, pow(NoH, uniformbuffer.specular_power));

	out_color = diffuse + (has_spec ? specular : vector3_t());

}


//
vector3_t F_fresenl_schlick(float VoH, vector3_t& f0)
{
	return f0 + (cOne - f0) * pow(1.0f - VoH, 5.0f);
}

// GGX / Trowbridge-Reitz
// [Walter et al. 2007, "Microfacet models for refraction through rough surfaces"]
float D_Trowbridge_Reitz_GGX(float a2, float NoH)
{
	float d = NoH * NoH * (a2 - 1.0f) + 1.0f;
	return a2 / (cPI * d * d);
}

// Tuned to match behavior of Vis_Smith
// [Schlick 1994, "An Inexpensive BRDF Model for Physically-Based Rendering"]
float V_Schlick_GGX(float a2, float NoV, float NoL)
{
	// V = G / (NoL * NoV)
	float k = sqrt(a2) * 0.5f;
	float Vis_SchlickV = NoV * (1.0f - k) + k;
	float Vis_SchlickL = NoL * (1.0f - k) + k;
	return 0.25f / (Vis_SchlickV * Vis_SchlickL);
}

// Smith term for GGX
// [Smith 1967, "Geometrical shadowing of a random rough surface"]
float V_smith_GGX(float a2, float NoV, float NoL)
{
	// V = G / (NoL * NoV)
	float Vis_SmithV = NoV + sqrt(NoV * (NoV - NoV * a2) + a2);
	float Vis_SmithL = NoL + sqrt(NoL * (NoL - NoL * a2) + a2);
	float m = Vis_SmithV * Vis_SmithL;
	return m != 0 ? 1.0f / m : FLT_MAX;
}

void pbr_shading(const interp_vertex_t& p, vector4_t& out_color)
{
	texcoord_t uv = p.uv;
	uv.u /= p.pos.w;
	uv.v /= p.pos.w;

	vector3_t wnor = p.nor / p.pos.w;
	wnor.normalize();

	vector3_t wpos = p.wpos / p.pos.w;

	vector3_t albedo = albedo_tex.sample(uv).to_vec3();
	vector4_t metallic_texel = metallic_tex.sample(uv);
	vector4_t roughness_texel = roughness_tex.sample(uv);

	vector3_t v = uniformbuffer.eye - wpos;
	v.normalize();

	vector3_t l = -uniformbuffer.light_dir;

	// (v + l) / 2
	vector3_t h = v + l;
	h.normalize();

	float NoL = max(dot(wnor, l), 0);
	float NoH = max(dot(wnor, h), 0);
	float NoV = max(dot(wnor, v), 0);

	float metallic = metallic_texel.r;
	float roughness = roughness_texel.r;
	float a = roughness * roughness;
	float a2 = a * a;

	vector3_t temp(0.04f, 0.04f, 0.04f);
	vector3_t f0 = lerp(temp, albedo, metallic);

	vector3_t F = F_fresenl_schlick(NoL, f0);

	float D = D_Trowbridge_Reitz_GGX(a2, NoH);

	float V = V_smith_GGX(a2, NoV, NoL);

	vector3_t cook_torrance_brdf = F * D * V / 4.0f;

	vector3_t ks = F;
	vector3_t kd = (cOne - F) * (1.0f - metallic);

	vector3_t I = uniformbuffer.light_intensity * NoL;

	vector3_t diffuse_brdf = albedo;// / cPI;

	vector3_t final_color = (kd * diffuse_brdf) * I;

	if (has_spec) {
		final_color += ks * cook_torrance_brdf * I;
	}

	reinhard_mapping(final_color);

	out_color = final_color;

}

void pixel_process(int x, int y, const interp_vertex_t& p)
{
	vector4_t color;
	if (shading_model == shading_model_t::cSM_Color)
	{
		color = p.color / p.pos.w;
	}
	else if (shading_model == shading_model_t::cSM_Phong)
	{
		phong_shading(p, color);
	}
	else if (shading_model == shading_model_t::cSM_PBR)
	{
		pbr_shading(p, color);
	}

	write_pixel(x, y, makefour(color));
	write_depth(x, y, p.pos.z);
}

void scan_horizontal(const interp_vertex_t& vl, const interp_vertex_t& vr, int y)
{
	float dist = vr.pos.x - vl.pos.x;
	// left要往小取整，right要往大取整，避免三角形之间的接缝空隙！
	int left = (int)(vl.pos.x);
	int right = (int)(vr.pos.x + 0.5f);
	for (int i = left; i < right; ++i)
	{
		float w = (i - vl.pos.x) / dist;
		w = clamp(w, 0.0f, 1.0f);
		interp_vertex_t p = lerp(vl, vr, w);
		if (depth_test(i, y, p.pos.z))
		{
			pixel_process(i, y, p);
		}
	}
}

void scan_triangle(scan_tri_t *sctri)
{
	if (sctri->l.pos.x > sctri->r.pos.x)
	{
		std::swap(sctri->l, sctri->r);
	}

	assert(sctri->l.pos.y == sctri->r.pos.y);
	assert(sctri->l.pos.y != sctri->p.pos.y);

	float ymax = std::fmax(sctri->p.pos.y, sctri->l.pos.y);
	float ymin = std::fmin(sctri->p.pos.y, sctri->l.pos.y);
	float ydist = sctri->p.pos.y - sctri->l.pos.y;

	// bottom要往小取整，top要往大取整，避免三角形之间的接缝空隙！
	int bottom = (int)(ymin);
	int top = (int)(ymax + 0.5f);
	for (int i = bottom; i < top; ++i)
	{
		float cury = i + 0.0f;
		float w = (ydist > 0 ? cury - ymin : cury - ymax) / ydist;
		w = clamp(w, 0.0f, 1.0f);
		interp_vertex_t vl = lerp(sctri->l, sctri->p, w);
		interp_vertex_t vr = lerp(sctri->r, sctri->p, w);
		scan_horizontal(vl, vr, i);
	}

}

bool check_clip(vector4_t* p, int width, int height)
{
	if (p->x < 0 || p->x >= height) return false;
	if (p->y < 0 || p->y >= height) return false;
	if (p->z < 0.0f || p->z > 1.0f) return false;
	return true;
}

void perspective_divide(interp_vertex_t* p)
{
	// the point with homogeneous coordinates [x, y, z, w] corresponds to the three-dimensional Cartesian point [x/w, y/w, z/w].
	// to cvv coord x[-1, 1], y[-1, 1], z[0, 1]
	float revw = 1.0f / p->pos.w;
	p->pos.x *= revw;
	p->pos.y *= revw;
	p->pos.z *= revw;
	p->pos.w = revw; // 这里存1/w，因为透视校正原因，1/w才有线性关系： 1/p.w = lerp(1/p0.w, 1/p1.w, weight)

	p->wpos.x *= revw;
	p->wpos.y *= revw;
	p->wpos.z *= revw;

	p->nor.x *= revw;
	p->nor.y *= revw;
	p->nor.z *= revw;

	p->uv.u *= revw;
	p->uv.v *= revw;

	p->color.r *= revw;
	p->color.g *= revw;
	p->color.b *= revw;
	p->color.a *= revw;

}

void to_screen_coord(vector4_t* p)
{
	// to screen coord x[0, width], y[0, height], z[0, 1] (depth)
	p->x = (p->x + 1.0f) * 0.5f * width;
	p->y = height - 1 - (p->y + 1.0f) * 0.5f * height;
}

void vertex_process(const matrix_t& world, const matrix_t& mvp, const model_vertex_t& v, interp_vertex_t& p)
{
	p.pos = mvp * vector4_t(v.pos, 1.0f);
	p.wpos = world * v.pos;
	p.nor = world * v.nor; // suppose world contain NO no-uniform scale!
	p.color = v.color;
	p.uv = v.uv;
	perspective_divide(&p);
	to_screen_coord(&p.pos);
}

void draw_triangle(const interp_vertex_t& p0, const interp_vertex_t& p1, const interp_vertex_t& p2)
{
	// degenerate triangle
	if (p0.pos.y == p1.pos.y && p1.pos.y == p2.pos.y) return;
	if (p0.pos.x == p1.pos.x && p1.pos.x == p2.pos.x) return;

	scan_tri_t uptri, downtri;
	if (p0.pos.y == p1.pos.y) // up triangle
	{
		uptri.p = p2;
		uptri.l = p0;
		uptri.r = p1;
		scan_triangle(&uptri);
	}
	else if (p1.pos.y == p2.pos.y) // down triangle
	{
		downtri.p = p0;
		downtri.l = p1;
		downtri.r = p2;
		scan_triangle(&downtri);
	}
	else
	{
		float w = (p1.pos.y - p0.pos.y) / (p2.pos.y - p0.pos.y);
		interp_vertex_t mid = lerp(p0, p2, w);
		mid.pos.y = p1.pos.y;

		uptri.p = p2;
		uptri.l = mid;
		uptri.r = p1;

		downtri.p = p0;
		downtri.l = mid;
		downtri.r = p1;

		scan_triangle(&uptri);
		scan_triangle(&downtri);

	}

}

void draw_line(const vector4_t& p0, const vector4_t& p1, uint32_t c)
{
	vector4_t v0 = p0;
	vector4_t v1 = p1;
	float dx = std::abs(v0.x - v1.x);
	float dy = std::abs(v0.y - v1.y);

	if (dy <= dx)
	{	
		if (v0.x > v1.x) {
			std::swap(v0, v1);
		}

		int xmin = (int)(v0.x);
		int xmax = (int)(v1.x + 0.5f);
		float y = v0.y;
		write_pixel(xmin, (int)y, c);

		float delta = (v1.y - v0.y) / dx;
		for (int x = xmin; x < xmax; ++x)
		{
			y += delta;
			write_pixel(x, (int)(y + 0.5f), c);
		}
	}
	else {
		if (v0.y > v1.y) {
			std::swap(v0, v1);
		}

		int ymin = (int)(v0.y);
		int ymax = (int)(v1.y + 0.5f);
		float x = v0.x;
		write_pixel((int)x, ymin, c);

		float delta = (v1.x - v0.x) / dy;
		for (int y = ymin; y < ymax; ++y)
		{
			x += delta;
			write_pixel((int)(x + 0.5f), y, c);
		}
	}

}

void draw_cartesian_coordinate(const matrix_t& mvp)
{
	vector3_t o(0, 0, 0);
	vector3_t x(1, 0, 0);
	vector3_t y(0, 1, 0);
	vector3_t z(0, 0, 1);

	auto transform_to_screen = [](const matrix_t& mvp, const vector3_t& p) -> vector4_t
	{
		vector4_t tp = mvp * vector4_t(p, 1.0f);
		to_screen_coord(&tp);
		return tp;
	};

	vector4_t to = transform_to_screen(mvp, o);
	vector4_t tx = transform_to_screen(mvp, x);
	vector4_t ty = transform_to_screen(mvp, y);
	vector4_t tz = transform_to_screen(mvp, z);

	draw_line(to, tx, 0xff0000ff);
	draw_line(to, ty, 0x00ff00ff);
	draw_line(to, tz, 0x0000ffff);
}

void update(model_t *model)
{
	memset(framebuffer, 0, width * height * sizeof(uint32_t));
	memset(zbuffer, 0, width * height * sizeof(float));

	//angle_speed += 0.008f;
	uniformbuffer.world.set_rotate(0, 0, 1, angle_speed);

	vector3_t at(0.0f, 0.0f, 0.0f);
	vector3_t up(0.0f, 0.0f, 1.0f);
	uniformbuffer.view.set_lookat(uniformbuffer.eye, at, up);

	matrix_t mv = mul(uniformbuffer.world, uniformbuffer.view);
	uniformbuffer.mvp = mul(mv, uniformbuffer.proj);

	draw_cartesian_coordinate(uniformbuffer.mvp);

	model_vertex_vec_t& vb = model->m_model_vertex;
	interp_vertex_vec_t& vb_post = model->m_vertex_post;
	index_vec_t& ib = model->m_model_indices;

	for (int i = 0; i < vb.size(); ++i)
	{
		vertex_process(uniformbuffer.world, uniformbuffer.mvp, vb[i], vb_post[i]);
	}

	int tri_num = (int)ib.size() / 3;
	for (int i = 0; i < tri_num; ++i)
	{
		int i0 = ib[i * 3 + 0];
		int i1 = ib[i * 3 + 1];
		int i2 = ib[i * 3 + 2];

		vector4_t v01 = vb_post[i1].pos - vb_post[i0].pos;
		vector4_t v02 = vb_post[i2].pos - vb_post[i0].pos;

		float det_xy = v01.x * v02.y - v01.y * v02.x;
		if (det_xy < 0.0f) {
			continue; // backface culling
		}

		interp_vertex_t* p0 = &vb_post[i0];
		interp_vertex_t* p1 = &vb_post[i1];
		interp_vertex_t* p2 = &vb_post[i2];

		if (!check_clip(&p0->pos, width, height)) return;
		if (!check_clip(&p1->pos, width, height)) return;
		if (!check_clip(&p2->pos, width, height)) return;

		// make sure p0y <= p1y <= p2y
		if (p0->pos.y > p1->pos.y) std::swap(p0, p1);
		if (p0->pos.y > p2->pos.y) std::swap(p0, p2);
		if (p1->pos.y > p2->pos.y) std::swap(p1, p2);

		if (shading_model == shading_model_t::cSM_Wireframe)
		{
			draw_line(p0->pos, p1->pos, 0xffffffff);
			draw_line(p1->pos, p2->pos, 0xffffffff);
			draw_line(p2->pos, p0->pos, 0xffffffff);
		}
		else
		{
			draw_triangle(*p0, *p1, *p2);
		}
	}

}


void main_loop()
{
	// Show the window
	::ShowWindow(hwnd, SW_SHOWDEFAULT);
	::UpdateWindow(hwnd);

	// Enter the message loop
	MSG msg;
	ZeroMemory(&msg, sizeof(msg));
	while (msg.message != WM_QUIT)
	{
		if (::PeekMessage(&msg, NULL, 0U, 0U, PM_REMOVE))
		{
			::TranslateMessage(&msg);
			::DispatchMessage(&msg);
		}
		else
		{
			int64_t tick_now = 0;
			::QueryPerformanceCounter((LARGE_INTEGER*)&tick_now);
			float dt_ms = (tick_now - tick_start)* reci_freq;
			++frame_count;
			if (dt_ms >= 1000.0f)
			{
				frame_rate = (uint32_t)(frame_count * 1000.0f / dt_ms);
				tick_start = tick_now;
				frame_count = 0;
				TCHAR str[64];
				::wsprintfW(str, _T("sr3d %d fps"), frame_rate);
				::SetWindowText(hwnd, str);
			}
			//update(&cube_model);
			update(&sphere_model);

			HDC hDC = GetDC(hwnd);
			BitBlt(hDC, 0, 0, width, height, screenDC, 0, 0, SRCCOPY);
			ReleaseDC(hwnd, hDC);

		}
	}
}

void update_light(float new_angle)
{
	light_angle = new_angle;
	uniformbuffer.light_dir.x = sin(light_angle);
	uniformbuffer.light_dir.y = cos(light_angle);
	uniformbuffer.light_dir.z = 0;
}

LRESULT CALLBACK MsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch (msg)
	{
	case WM_CREATE:
		break;
	case WM_KEYDOWN:
		switch (wParam & 511)
		{
		case 'C':
			shading_model = (shading_model_t)((int)shading_model + 1);
			if (shading_model == shading_model_t::cSM_MAX)
				shading_model = shading_model_t::cSM_Color;
			std::cout << (int)shading_model << std::endl;
			break;
		case 'S':
			has_spec = !has_spec;
			break;
		case VK_UP:
			eyedist += 0.1f;
			break;
		case VK_DOWN:
			eyedist -= 0.1f;
			break;
		case VK_LEFT:
			update_light(light_angle + 0.1f);
			break;
		case VK_RIGHT:
			update_light(light_angle - 0.1f);
			break;
		default:
			break;
		}
		break;
	case WM_SIZE:
	case WM_EXITSIZEMOVE:
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProc(hWnd, msg, wParam, lParam);
}

HWND init_window(HINSTANCE instance, const TCHAR* title, int width, int height)
{
	const TCHAR class_name[] = _T("wndclass");

	// Register the window class
	WNDCLASSEX wc = {
		sizeof(WNDCLASSEX), CS_CLASSDC, MsgProc, 0, 0,
		instance, NULL, NULL, NULL, NULL,
		class_name, NULL
	};

	RegisterClassEx(&wc);

	// calculate client size，设置窗口客户区大小
	RECT clientSize;
	clientSize.top = 0;
	clientSize.left = 0;
	clientSize.right = width;
	clientSize.bottom = height;
	DWORD style = WS_OVERLAPPEDWINDOW | WS_BORDER | WS_CAPTION | WS_CLIPCHILDREN | WS_CLIPSIBLINGS;
	//计算客户区矩形大小
	::AdjustWindowRect(&clientSize, style, FALSE);
	int w = clientSize.right - clientSize.left;
	int h = clientSize.bottom - clientSize.top;

	// Create the application's window
	HWND hwnd = CreateWindow(class_name, title,
		style, 0, 0, w, h,
		NULL, NULL, instance, NULL);

	::ShowWindow(hwnd, SW_SHOWDEFAULT);
	::UpdateWindow(hwnd);

	return hwnd;

}

int main(void)
{
	width = 800;
	height = 600;
	hwnd = init_window(GetModuleHandle(NULL), _T("sr3d"), width, height);

	screenDC = CreateCompatibleDC(GetDC(hwnd));

	BITMAPINFO bi = {
		{ sizeof(BITMAPINFOHEADER), width, height, 1u, 32u, BI_RGB, width * height * 4u, 0, 0, 0, 0 }
	};

	HBITMAP screenBMP = CreateDIBSection(screenDC, &bi, DIB_RGB_COLORS, (void**)&framebuffer, 0, 0);
	if (screenBMP == NULL) {
		return -1;
	}

	SelectObject(screenDC, screenBMP);

	LARGE_INTEGER temp;
	::QueryPerformanceFrequency(&temp);
	reci_freq = 1000.0f / temp.QuadPart;
	::QueryPerformanceCounter((LARGE_INTEGER*)&tick_start);

	zbuffer = new float[width * height];
	memset(zbuffer, 0, width * height * sizeof(float));

	albedo_tex.load_tex("./rustediron2_basecolor.png");
	metallic_tex.load_tex("./rustediron2_metallic.png");
	roughness_tex.load_tex("./rustediron2_roughness.png");
	normal_tex.load_tex("./rustediron2_normal.png");

	uniformbuffer.eye.set(0.0f, eyedist, 0);
	uniformbuffer.light_dir.set(0.0f, -1.0f, 0.0f);
	uniformbuffer.light_intensity.set(1.0f, 1.0f, 1.0f);
	uniformbuffer.specular_power = 8.0f;

	float aspect = 1.0f * width / height;
	uniformbuffer.proj.set_perspective(cPI * 0.5f, aspect, 1.0f, 500.0f);

	sphere_model.load("mesh_sphere.obj");

	update_light(light_angle);
	main_loop();

	return 0;

}
