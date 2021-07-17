#include <windows.h>
#include <tchar.h>
#include <iostream>
#include <cassert>

const float cEpslion = 1e-6f;
const float cPI = 3.1415926f;
const float cRevt255 = 1.0f / 255.0f;

struct matrix_t
{
	float m[4][4];
};

struct vector_t
{
	union
	{
		float vec[4];
		struct
		{
			float x, y, z, w;
		};
		struct
		{
			float r, g, b, a;
		};
	};
};

struct texcoord_t
{
	float u, v;
};

struct vertex_t
{
	vector_t pos;
	texcoord_t uv;
	vector_t color;
};

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
typedef struct {
	vertex_t p;
	vertex_t l, r;
} scan_tri_t;

enum class shading_model_t {
	cSM_Color = 0,
	cSM_Texture = 1,
};

HWND hwnd = 0;
int width = 0, height = 0;
float angle_speed = 0.85f;
matrix_t model;
matrix_t view;
matrix_t proj;
matrix_t mvp;

uint32_t *framebuffer = nullptr;
float *zbuffer = nullptr;
HDC screenDC;
uint32_t *texture;
int tex_width, tex_height;
shading_model_t shading_model = shading_model_t::cSM_Texture;

template <typename T, int n>
int array_size(T(&)[n])
{
	return n;
}

template <typename T>
T clamp(T x, T min, T max)
{
	return (x < min) ? min : ((x > max) ? max : x);
}

uint8_t to_color_int(float c)
{
	int cint = (int)(c * 255.0f + 0.5);
	cint = clamp(cint, 0, 255);
	return (uint8_t)cint;
}

uint32_t makefour(const vector_t& color)
{
	return to_color_int(color.r)
		| to_color_int(color.g) << 8
		| to_color_int(color.b) << 16
		| to_color_int(color.a) << 24;
}

void to_color(uint32_t cint, vector_t& color)
{
	color.r = cRevt255 * (cint & 0xff);
	color.g = cRevt255 * ((cint >> 8) & 0xff);
	color.b = cRevt255 * ((cint >> 16) & 0xff);
	color.a = cRevt255 * ((cint >> 24) & 0xff);
}

uint32_t lerp(uint32_t x1, uint32_t x2, float t)
{
	return (uint32_t)(x1 * (1.0f - t) + x2 * t);
}

float lerp(float x1, float x2, float t)
{ 
	return x1 + (x2 - x1) * t;
}

void lerp(texcoord_t* p, const texcoord_t* a, const texcoord_t* b, float w)
{
	p->u = lerp(a->u, b->u, w);
	p->v = lerp(a->v, b->v, w);
}

void lerp(vector_t* p, const vector_t* a, const vector_t* b, float w)
{
	p->x = lerp(a->x, b->x, w);
	p->y = lerp(a->y, b->y, w);
	p->z = lerp(a->z, b->z, w);
	p->w = lerp(a->w, b->w, w);
}

void lerp(vertex_t* p, const vertex_t* a, const vertex_t* b, float w)
{
	lerp(&p->pos, &a->pos, &b->pos, w);
	lerp(&p->uv, &a->uv, &b->uv, w);
	lerp(&p->color, &a->color, &b->color, w);
}

// | v |
float vector_length(const vector_t* v)
{
	float sq = v->x * v->x + v->y * v->y + v->z * v->z;
	return (float)sqrt(sq);
}

// c = a + b
void vector_add(vector_t* c, const vector_t* a, const vector_t* b)
{
	c->x = a->x + b->x;
	c->y = a->y + b->y;
	c->z = a->z + b->z;
	c->w = 1.0;
}

// c = a - b
void vector_sub(vector_t* c, const vector_t* a, const vector_t* b)
{
	c->x = a->x - b->x;
	c->y = a->y - b->y;
	c->z = a->z - b->z;
	c->w = 1.0;
}

// c = a * f
void vector_scale(vector_t* c, const vector_t* a, float f)
{
	c->x = a->x * f;
	c->y = a->y * f;
	c->z = a->z * f;
	c->w = 1.0;
}

// 矢量点乘
float vector_dot(const vector_t* a, const vector_t* b)
{
	return a->x * b->x + a->y * b->y + a->z * b->z;
}

// 矢量叉乘
void vector_cross(vector_t* c, const vector_t* a, const vector_t* b)
{
	c->x = a->y * b->z - a->z * b->y;
	c->y = a->z * b->x - a->x * b->z;
	c->z = a->x * b->y - a->y * b->x;
	c->w = 1.0f;
}

// 矢量插值，t取值 [0, 1]
void vector_lerp(vector_t* c, const vector_t* a, const vector_t* b, float t)
{
	c->x = lerp(a->x, b->x, t);
	c->y = lerp(a->y, b->y, t);
	c->z = lerp(a->z, b->z, t);
	c->w = 1.0f;
}

// 矢量归一化
void vector_normalize(vector_t* v)
{
	float length = vector_length(v);
	float inv = cEpslion;
	if (abs(length) > cEpslion)
	{
		inv = 1.0f / length;
	}

	v->x *= inv;
	v->y *= inv;
	v->z *= inv;

}

// c = a + b
void matrix_add(matrix_t* c, const matrix_t* a, const matrix_t* b)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
			c->m[i][j] = a->m[i][j] + b->m[i][j];
	}
}

// c = a - b
void matrix_sub(matrix_t* c, const matrix_t* a, const matrix_t* b)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
			c->m[i][j] = a->m[i][j] + b->m[i][j];
	}
}

// c = a * b
void matrix_mul(matrix_t* c, const matrix_t* a, const matrix_t* b)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			float s = 0;
			for (int k = 0; k < 4; ++k)
			{
				s += (a->m[i][k] * b->m[k][j]);
			}
			c->m[i][j] = s;
		}
	}
}

// c = a * f
void matrix_scale(matrix_t* c, const matrix_t* a, float f)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			c->m[i][j] = a->m[i][j] * f;
		}
	}
}

// y = x * m
void matrix_apply(vector_t* y, const vector_t* x, const matrix_t* m)
{
	float X = x->x, Y = x->y, Z = x->z, W = x->w;
	y->x = X * m->m[0][0] + Y * m->m[1][0] + Z * m->m[2][0] + W * m->m[3][0];
	y->y = X * m->m[0][1] + Y * m->m[1][1] + Z * m->m[2][1] + W * m->m[3][1];
	y->z = X * m->m[0][2] + Y * m->m[1][2] + Z * m->m[2][2] + W * m->m[3][2];
	y->w = X * m->m[0][3] + Y * m->m[1][3] + Z * m->m[2][3] + W * m->m[3][3];
}

void matrix_set_identity(matrix_t* m)
{
	m->m[0][0] = m->m[1][1] = m->m[2][2] = m->m[3][3] = 1.0f;
	m->m[0][1] = m->m[0][2] = m->m[0][3] = 0.0f;
	m->m[1][0] = m->m[1][2] = m->m[1][3] = 0.0f;
	m->m[2][0] = m->m[2][1] = m->m[2][3] = 0.0f;
	m->m[3][0] = m->m[3][1] = m->m[3][2] = 0.0f;
}

void matrix_set_zero(matrix_t* m) {
	m->m[0][0] = m->m[0][1] = m->m[0][2] = m->m[0][3] = 0.0f;
	m->m[1][0] = m->m[1][1] = m->m[1][2] = m->m[1][3] = 0.0f;
	m->m[2][0] = m->m[2][1] = m->m[2][2] = m->m[2][3] = 0.0f;
	m->m[3][0] = m->m[3][1] = m->m[3][2] = m->m[3][3] = 0.0f;
}

// 平移变换
void matrix_set_translate(matrix_t* m, float x, float y, float z)
{
	matrix_set_identity(m);
	m->m[3][0] = x;
	m->m[3][1] = y;
	m->m[3][2] = z;
}

// 缩放变换
void matrix_set_scale(matrix_t* m, float x, float y, float z)
{
	matrix_set_identity(m);
	m->m[0][0] = x;
	m->m[1][1] = y;
	m->m[2][2] = z;
}

// 旋转矩阵
void matrix_set_rotate(matrix_t* m, float x, float y, float z, float theta)
{
	float qsin = (float)sin(theta * 0.5f);
	float qcos = (float)cos(theta * 0.5f);
	vector_t vec = { x, y, z, 1.0f };
	float w = qcos;
	vector_normalize(&vec);
	x = vec.x * qsin;
	y = vec.y * qsin;
	z = vec.z * qsin;
	m->m[0][0] = 1 - 2 * y * y - 2 * z * z;
	m->m[1][0] = 2 * x * y - 2 * w * z;
	m->m[2][0] = 2 * x * z + 2 * w * y;
	m->m[0][1] = 2 * x * y + 2 * w * z;
	m->m[1][1] = 1 - 2 * x * x - 2 * z * z;
	m->m[2][1] = 2 * y * z - 2 * w * x;
	m->m[0][2] = 2 * x * z - 2 * w * y;
	m->m[1][2] = 2 * y * z + 2 * w * x;
	m->m[2][2] = 1 - 2 * x * x - 2 * y * y;
	m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.0f;
	m->m[3][0] = m->m[3][1] = m->m[3][2] = 0.0f;
	m->m[3][3] = 1.0f;
}

// view matrix
void matrix_set_lookat(matrix_t* m, const vector_t* eye, const vector_t* at, const vector_t* up)
{
	vector_t xaxis, yaxis, zaxis;

	vector_sub(&zaxis, at, eye);
	vector_normalize(&zaxis);
	vector_cross(&xaxis, up, &zaxis);
	vector_normalize(&xaxis);
	vector_cross(&yaxis, &zaxis, &xaxis);

	m->m[0][0] = xaxis.x;
	m->m[1][0] = xaxis.y;
	m->m[2][0] = xaxis.z;
	m->m[3][0] = -vector_dot(&xaxis, eye);

	m->m[0][1] = yaxis.x;
	m->m[1][1] = yaxis.y;
	m->m[2][1] = yaxis.z;
	m->m[3][1] = -vector_dot(&yaxis, eye);

	m->m[0][2] = zaxis.x;
	m->m[1][2] = zaxis.y;
	m->m[2][2] = zaxis.z;
	m->m[3][2] = -vector_dot(&zaxis, eye);

	m->m[0][3] = m->m[1][3] = m->m[2][3] = 0.0f;
	m->m[3][3] = 1.0f;
}

// projection matrix (ref to D3DXMatrixPerspectiveFovLH)
void matrix_set_perspective(matrix_t* m, float fovy, float aspect, float zn, float zf)
{
	float fax = 1.0f / (float)tan(fovy * 0.5f);
	matrix_set_zero(m);
	m->m[0][0] = (float)(fax / aspect);
	m->m[1][1] = (float)(fax);
	m->m[2][2] = zf / (zf - zn);
	m->m[3][2] = -zn * zf / (zf - zn);
	m->m[2][3] = 1;
}

uint32_t texture_sample(const texcoord_t& texcoord)
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
	
	return makefour(c);

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

void pixel_process(int x, int y, const vertex_t& p)
{
	uint32_t c = 0;
	if (shading_model == shading_model_t::cSM_Color)
	{
		vector_t color;
		vector_scale(&color, &p.color, p.pos.w);
		c = makefour(color);
	}
	else if (shading_model == shading_model_t::cSM_Texture)
	{
		texcoord_t uv = p.uv;
		uv.u *= p.pos.w;
		uv.v *= p.pos.w;
		c = texture_sample(uv);
	}

	write_pixel(x, y, c);
	write_depth(x, y, p.pos.z);
}

void scan_horizontal(vertex_t* vl, vertex_t* vr, int y)
{
	float dist = vr->pos.x - vl->pos.x;
	// left要往小取整，right要往大取整，避免三角形之间的接缝空隙！
	int left = (int)(vl->pos.x);
	int right = (int)(vr->pos.x + 0.5f);
	for (int i = left; i < right; ++i)
	{
		vertex_t p;
		float w = (i - vl->pos.x) / dist;
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

	float ymax = std::fmax(sctri->p.pos.y, sctri->l.pos.y);
	float ymin = std::fmin(sctri->p.pos.y, sctri->l.pos.y);
	float ydist = sctri->p.pos.y - sctri->l.pos.y;

	// bottom要往小取整，top要往大取整，避免三角形之间的接缝空隙！
	int bottom = (int)(ymin + 0.0f);
	int top = (int)(ymax + 0.5f);

	vertex_t vl, vr;
	for (int i = bottom; i < top; ++i)
	{
		float cury = i + 0.0f;
		float w = ydist > 0 ? (cury - ymin) / ydist : (cury - ymax) / ydist;
		w = clamp(w, 0.0f, 1.0f);
		lerp(&vl, &sctri->l, &sctri->p, w);
		lerp(&vr, &sctri->r, &sctri->p, w);
		scan_horizontal(&vl, &vr, i);
	}

}


bool check_clip(vector_t* p)
{
	if (p->x < -1.0f || p->x > 1.0f) return false;
	if (p->y < -1.0f || p->y > 1.0f) return false;
	if (p->z < 0.0f || p->z > 1.0f) return false;
	return true;
}

void perspective_divide(vertex_t* p)
{
	// the point with homogeneous coordinates [x, y, z, w] corresponds to the three-dimensional Cartesian point [x/w, y/w, z/w].
	// to cvv coord x[-1, 1], y[-1, 1], z[0, 1]
	float revw = 1.0f / p->pos.w;
	p->pos.x *= revw;
	p->pos.y *= revw;
	p->pos.z *= revw;

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

void vertex_process(const matrix_t* mvp, const vertex_t& v, vertex_t& p)
{
	p = v;
	matrix_apply(&p.pos, &v.pos, mvp);
	perspective_divide(&p);
	to_screen_coord(&p.pos);
}

void draw_triangle(const matrix_t* mvp, const vertex_t& p0, const vertex_t& p1, const vertex_t& p2)
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
		vertex_t mid;
		mid.pos.x = mid.pos.y = mid.pos.z = mid.pos.w = 0;
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

vertex_t box_vb[12] = {
	{ { -1, -1, -1, 1 }, { 0, 0 }, { 1, 0, 0 } },
	{ { -1,  1, -1, 1 }, { 0, 1 }, { 0, 1, 0 } },

	{ { -1, -1,  1, 1 }, { 1, 0 }, { 0, 0, 1 } },
	{ { -1,  1,  1, 1 }, { 1, 1 }, { 1, 0, 0 } },

	{ { 1, -1,  1, 1 }, { 0, 0 }, { 1, 0, 0 } },
	{ { 1,  1,  1, 1 }, { 0, 1 }, { 0, 1, 0 } },

	{ { 1, -1, -1, 1 }, { 1, 0 }, { 0, 0, 1 } },
	{ { 1,  1, -1, 1 }, { 1, 1 }, { 1, 0, 0 } },

	{ { -1,  1, -1, 1 }, { 1, 0 }, { 0, 1, 0 } },
	{ { 1,  1, -1, 1 }, { 0, 0 }, { 1, 0, 0 } },

	{ { -1, -1, -1, 1 }, { 1, 1 }, { 1, 0, 0 } },
	{ { 1, -1, -1, 1 }, { 0, 1 }, { 0, 0, 1 } },

};

int box_ib1[] = {
	0, 1, 3,  3, 2, 0,
	2, 3, 5,  5, 4, 2,
	4, 5, 7,  7, 6, 4,
	6, 7, 1,  1, 0, 6,
	3, 8, 9,  9, 5, 3,
	2, 4, 11, 11, 10, 2,
};

int box_ib[] = {
	3, 8, 9, 
	9, 5, 3,
};

vertex_t vb_post[12]; // vertex buffer post transform

void update()
{
	memset(framebuffer, 0, width * height * sizeof(uint32_t));
	memset(zbuffer, 0, width * height * sizeof(float));

//	angle_speed += 0.0001f;
	matrix_set_rotate(&model, 0.5f, 0.5f, 0.5f, cPI *  1.0f);

	matrix_t mv;
	matrix_mul(&mv, &model, &view);
	matrix_mul(&mvp, &mv, &proj);

	for (int i = 0; i < array_size(box_vb); ++i)
	{	
		vertex_process(&mvp, box_vb[i], vb_post[i]);
	}

	int tri_num = array_size(box_ib) / 3;
	for (int i = 0; i < tri_num; ++i)
	{
		int i0 = box_ib[i * 3 + 0];
		int i1 = box_ib[i * 3 + 1];
		int i2 = box_ib[i * 3 + 2];

		//if (!check_clip(&p0.pos)) return;
		//if (!check_clip(&p1.pos)) return;
		//if (!check_clip(&p2.pos)) return;

		vector_t v01, v02;
		vector_sub(&v01, &vb_post[i1].pos, &vb_post[i0].pos);
		vector_sub(&v02, &vb_post[i2].pos, &vb_post[i0].pos);

		float det_xy = v01.x * v02.y - v01.y * v02.x;
		if (det_xy > 0.0f)
		{
			// backface culling
			continue;
		}

		vertex_t* p0 = &vb_post[i0];
		vertex_t* p1 = &vb_post[i1];
		vertex_t* p2 = &vb_post[i2];

		// make sure p0y <= p1y <= p2y
		if (p0->pos.y > p1->pos.y) std::swap(p0, p1);
		if (p0->pos.y > p2->pos.y) std::swap(p0, p2);
		if (p1->pos.y > p2->pos.y) std::swap(p1, p2);

		draw_triangle(&mvp, *p0, *p1, *p2);

	}

	HDC hDC = GetDC(hwnd);
	BitBlt(hDC, 0, 0, width, height, screenDC, 0, 0, SRCCOPY);
	ReleaseDC(hwnd, hDC);
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
			update();
		}
	}
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

	zbuffer = new float[width * height];
	memset(zbuffer, 0, width * height * sizeof(float));

	tex_width = 512;
	tex_height = 512;
	texture = new uint32_t[tex_width * tex_height];
	for (int i = 0; i < tex_height; ++i)
	{
		for (int j = 0; j < tex_width; ++j)
		{
			int x = i / 32, y = j / 32;
			texture[i * tex_width + j] = ((x + y) & 1) ? 0xffffff : 0x3fbcef;
		}
	}

	screenDC = CreateCompatibleDC(GetDC(hwnd));

	BITMAPINFO bi = {
		{ sizeof(BITMAPINFOHEADER), width, -height, 1, 32, BI_RGB, width * height * 4, 0, 0, 0, 0 }
	};

	HBITMAP screenBMP = CreateDIBSection(screenDC, &bi, DIB_RGB_COLORS, (void**)&framebuffer, 0, 0);
	if (screenBMP == NULL)
	{
		return -1;
	}

	SelectObject(screenDC, screenBMP);

	float aspect = 1.0f * width / height;
	matrix_set_perspective(&proj, cPI * 0.5f, aspect, 1.0f, 500.0f);

	vector_t eye = { 3.0f, 0.0f, 0.0f, 1.0f };
	vector_t at = { 0.0f, 0.0f, 0.0f, 1.0f };
	vector_t up = { 0.0f, 0.0f, 1.0f, 1.0f };
	matrix_set_lookat(&view, &eye, &at, &up);

	main_loop();

	return 0;

}
