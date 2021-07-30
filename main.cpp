#include <windows.h>
#include <tchar.h>
#include <vector>
#include <iostream>
#include <cassert>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "math3d.h"

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
	cSM_Phong = 1,
	cSM_PBR = 2,
};


float reci_freq;
int64_t tick_start;
uint32_t frame_count = 0, frame_rate = 0;
HWND hwnd = 0;
HDC screenDC;

int width = 0, height = 0;
float angle_speed = 1.0f;
float light_angle = 0;
float eyedist = 3.5f;

struct uniformbuffer_t
{
	vector_t eye;
	matrix_t world;
	matrix_t view;
	matrix_t proj;
	matrix_t mvp;
	vector_t light_dir;
	float light_intensity;
	float specular_power;
};

uniformbuffer_t uniformbuffer;
uint32_t *framebuffer = nullptr;
float *zbuffer = nullptr;
uint32_t *texture;
int tex_width, tex_height;
shading_model_t shading_model = shading_model_t::cSM_Phong;

sphere_t sphere_model(12);
cube_t cube_model;

template <typename T, int n>
int array_size(T(&)[n])
{
	return n;
}

vector_t texture_sample(const texcoord_t& texcoord)
{
	float u = clamp(texcoord.u, 0.0f, 1.0f) * (tex_width - 1);
	float v = clamp(texcoord.v, 0.0f, 1.0f) * (tex_height - 1);
	int x0 = (int)(u);
	int y0 = (int)(v);
	int x1 = clamp(x0 + 1, 0, tex_width - 1);
	int y1 = clamp(y0 + 1, 0, tex_height - 1);
	uint32_t t00 = texture[y0 * tex_width + x0];
	uint32_t t01 = texture[y1 * tex_width + x0];
	uint32_t t10 = texture[y0 * tex_width + x1];
	uint32_t t11 = texture[y1 * tex_width + x1];

	float u_weight = u - x0;
	float v_weight = v - y0;

	vector_t c00, c10, c01, c11;
	to_color(t00, c00);	
	to_color(t10, c10);
	to_color(t01, c01);
	to_color(t11, c11);

	// bilinear interpolation
	vector_t tu0, tu1, c;
	lerp(&tu0, &c00, &c10, u_weight);
	lerp(&tu1, &c01, &c11, u_weight);
	lerp(&c, &tu0, &tu1, v_weight);

	return c;

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

void phong_shading(const interp_vertex_t& p, vector_t& out_color)
{
	texcoord_t uv = p.uv;
	uv.u /= p.pos.w;
	uv.v /= p.pos.w;

	vector_t wnor;
	vector_scale(&wnor, &p.nor, 1 / p.pos.w);
	vector_normalize(&wnor);

	vector_t wpos;
	vector_scale(&wpos, &p.wpos, 1 / p.pos.w);

	vector_t albedo = texture_sample(uv);

	float NoL = vector_dot(&wnor, &uniformbuffer.light_dir);
	NoL = max(NoL, 0.15f);

	vector_t diffuse;
	vector_scale(&diffuse, &albedo, uniformbuffer.light_intensity * NoL);

	vector_t specular;
	vector_t v;
	vector_sub(&v, &uniformbuffer.eye, &p.wpos);

	vector_t h; // (v + l) / 2
	vector_add(&h, &v, &uniformbuffer.light_dir);
	vector_normalize(&h);

	float NoH = vector_dot(&wnor, &h);
	vector_t ks(0.5f, 0.5f, 0.5f);
	vector_scale(&specular, &ks, uniformbuffer.light_intensity * max(0.0f, pow(NoH, uniformbuffer.specular_power)));

	vector_add(&out_color, &diffuse, &specular);

}

void pbr_shading(const interp_vertex_t& p, vector_t& out_color)
{

}

void pixel_process(int x, int y, const interp_vertex_t& p)
{
	vector_t color;
	if (shading_model == shading_model_t::cSM_Color)
	{		
		vector_scale(&color, &p.color, 1 / p.pos.w);
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

void scan_horizontal(interp_vertex_t* vl, interp_vertex_t* vr, int y)
{
	float dist = vr->pos.x - vl->pos.x;
	// left要往小取整，right要往大取整，避免三角形之间的接缝空隙！
	int left = (int)(vl->pos.x);
	int right = (int)(vr->pos.x + 0.5f);
	for (int i = left; i < right; ++i)
	{
		interp_vertex_t p;
		float w = (i - vl->pos.x) / dist;
		w = clamp(w, 0.0f, 1.0f);
		lerp(&p, vl, vr, w);
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

	interp_vertex_t vl, vr;
	for (int i = bottom; i < top; ++i)
	{
		float cury = i + 0.0f;
		float w = (ydist > 0 ? cury - ymin : cury - ymax) / ydist;
		w = clamp(w, 0.0f, 1.0f);
		lerp(&vl, &sctri->l, &sctri->p, w);
		lerp(&vr, &sctri->r, &sctri->p, w);
		scan_horizontal(&vl, &vr, i);
	}

}

bool check_clip(vector_t* p, int width, int height)
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
	p->wpos.w = 1.0f;

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

void to_screen_coord(vector_t* p)
{
	// to screen coord x[0, width], y[0, height], z[0, 1] (depth)
	p->x = (p->x + 1.0f) * 0.5f * width;
	p->y = height - 1 - (p->y + 1.0f) * 0.5f * height;
}

void vertex_process(const matrix_t* world, const matrix_t* mvp, const model_vertex_t& v, interp_vertex_t& p)
{
	matrix_apply(&p.pos, &v.pos, mvp);
	matrix_apply(&p.wpos, &v.pos, world);
	matrix_apply(&p.nor, &v.nor, world); // suppose world contain NO no-uniform scale!
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
		interp_vertex_t mid;
		float w = (p1.pos.y - p0.pos.y) / (p2.pos.y - p0.pos.y);
		lerp(&mid, &p0, &p2, w);
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

void update(model_base_t *model)
{
	memset(framebuffer, 0, width * height * sizeof(uint32_t));
	memset(zbuffer, 0, width * height * sizeof(float));

	//angle_speed += 0.008f;
	matrix_set_rotate(&uniformbuffer.world, 0, 0, 1, angle_speed);

	vector_t at = { 0.0f, 0.0f, 0.0f, 1.0f };
	vector_t up = { 0.0f, 0.0f, 1.0f, 1.0f };
	matrix_set_lookat(&uniformbuffer.view, &uniformbuffer.eye, &at, &up);

	matrix_t mv;
	matrix_mul(&mv, &uniformbuffer.world, &uniformbuffer.view);
	matrix_mul(&uniformbuffer.mvp, &mv, &uniformbuffer.proj);

	model_vertex_vec_t& vb = model->m_model_vertex;
	interp_vertex_vec_t& vb_post = model->m_vertex_post;
	index_vec_t& ib = model->m_model_indices;

	for (int i = 0; i < vb.size(); ++i)
	{
		vertex_process(&uniformbuffer.world, &uniformbuffer.mvp, vb[i], vb_post[i]);
	}

	int tri_num = (int)ib.size() / 3;
	for (int i = 0; i < tri_num; ++i)
	{
		int i0 = ib[i * 3 + 0];
		int i1 = ib[i * 3 + 1];
		int i2 = ib[i * 3 + 2];

		vector_t v01, v02;
		vector_sub(&v01, &vb_post[i1].pos, &vb_post[i0].pos);
		vector_sub(&v02, &vb_post[i2].pos, &vb_post[i0].pos);

		float det_xy = v01.x * v02.y - v01.y * v02.x;
		if (det_xy > 0.0f)
		{
			// backface culling
			continue;
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

		draw_triangle(*p0, *p1, *p2);

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
	uniformbuffer.light_dir.x = 0.96f * sin(light_angle);
	uniformbuffer.light_dir.y = 0.96f * cos(light_angle);
	uniformbuffer.light_dir.z = -0.2f;
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
			if (shading_model == shading_model_t::cSM_Color)
				shading_model = shading_model_t::cSM_Texture;
			else if (shading_model == shading_model_t::cSM_Texture)
				shading_model = shading_model_t::cSM_Color;
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


void load_tex(const char* tex_path, uint32_t*& src_tex, int& tex_w, int& tex_h)
{
	int components = 0;
	stbi_uc* st_img = stbi_load(tex_path, &tex_w, &tex_h, &components, STBI_rgb_alpha);
	if (st_img == nullptr)
	{
		// if we haven't returned, it's because we failed to load the file.
		printf("Failed to load image %s\nReason: %s\n", tex_path, stbi_failure_reason());
		return;
	}

	src_tex = new uint32_t[tex_w * tex_h];

	for (int i = 0; i < tex_h; ++i) {
		for (int j = 0; j < tex_w; ++j) {
			int pidx = (tex_h - 1 - i) * tex_w + j;
			vector_t c;
			c.b = st_img[pidx * 4 + 0] / 255.0f;
			c.g = st_img[pidx * 4 + 1] / 255.0f;
			c.r = st_img[pidx * 4 + 2] / 255.0f;
			c.a = st_img[pidx * 4 + 3] / 255.0f;
			src_tex[pidx] = makefour(c);
		}
	}

	stbi_image_free(st_img);
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

	load_tex("./albedo.png", texture, tex_width, tex_height);

	uniformbuffer.eye = { eyedist, 0.0f, 0.0f, 1.0f };
	uniformbuffer.light_dir = (1.0f, 0.0f, 0.0f);
	uniformbuffer.light_intensity = 2.0f;
	uniformbuffer.specular_power = 5.0f;

	float aspect = 1.0f * width / height;
	matrix_set_perspective(&uniformbuffer.proj, cPI * 0.5f, aspect, 1.0f, 500.0f);

	update_light(light_angle);
	main_loop();

	return 0;

}
