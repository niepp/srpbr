#ifndef __MATH3D_H
#define __MATH3D_H

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
	vector_t(float _x = 0, float _y = 0, float _z = 0, float _w = 1.0f):
		x(_x), y(_y), z(_z), w(_w) 
	{ }
};

struct texcoord_t
{
	float u, v;
	texcoord_t(float _u = 0, float _v = 0) :
		u(_u), v(_v)
	{ }
};

struct vertex_t
{
	vector_t pos;
	vector_t nor;
	texcoord_t uv;
	vector_t color;
};

typedef std::vector<vertex_t> vertex_vec_t;
typedef std::vector<uint16_t> index_vec_t;

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
	lerp(&p->nor, &a->nor, &b->nor, w);
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

struct model_base_t {
	vertex_vec_t m_model_vertex;
	index_vec_t m_model_indices;
	vertex_vec_t m_vertex_post;
};


/*
	unit sphere with center at (0,0,0) and radius of 1.0
*/
struct sphere_t : public model_base_t {

	sphere_t(uint32_t seg) {
		gen_sphere(seg, seg * 2);
	}
private:
	void assign_vertex(uint32_t v_idx, const vector_t& p);
	void gen_sphere(uint32_t seg1, uint32_t seg2);
};

inline void sphere_t::assign_vertex(uint32_t v_idx, const vector_t& p)
{
	m_model_vertex[v_idx].pos = { p.x, p.y, p.z, 1.0f };
	m_model_vertex[v_idx].nor = { p.x, p.y, p.z, 1.0f };
	m_model_vertex[v_idx].uv = { p.x * 0.5f + 0.5f, p.y * 0.5f + 0.5f };
	m_model_vertex[v_idx].color = { 1, 0, 0, 1 };
}

inline void sphere_t::gen_sphere(uint32_t seg1, uint32_t seg2)
{
	uint32_t v_count = (seg1 - 1) * seg2 + 2;
	m_model_vertex.resize(v_count);
	m_vertex_post.resize(v_count);

	assign_vertex(0, vector_t(0, 0, 1.0f));
	for (uint32_t i = 0; i < seg1 - 1; ++i)
	{
		float alpha = cPI * (i + 1) / seg1;
		float radius = sin(alpha);
		float z = cos(alpha);
		for (uint32_t j = 0; j < seg2; ++j)
		{
			float beta = 2.0f * cPI * j / seg2;
			float x = radius * cos(beta);
			float y = radius * sin(beta);
			uint32_t v_idx = i * seg2 + j + 1;
			assign_vertex(v_idx, vector_t(x, y, z));
		}
	}

	assign_vertex(v_count - 1, vector_t(0, 0, -1.0f));
	m_model_indices.resize((seg1 * 2 - 2) * 3 * seg2);

	uint32_t tri = 0;
	for (uint32_t j = 1; j < seg2 + 1; ++j)
	{
		m_model_indices[tri * 3 + 0] = 0;
		m_model_indices[tri * 3 + 1] = j;
		m_model_indices[tri * 3 + 2] = (j % seg2) + 1;
		++tri;
	}

	for (uint32_t i = 1; i < seg1 - 1; ++i)
	{
		for (uint32_t j = 1; j < seg2 + 1; ++j)
		{
			m_model_indices[tri * 3 + 0] = (i - 1) * seg2 + j;
			m_model_indices[tri * 3 + 1] = i * seg2 + j;
			m_model_indices[tri * 3 + 2] = i * seg2 + (j % seg2) + 1;
			++tri;
			m_model_indices[tri * 3 + 0] = (i - 1) * seg2 + (j % seg2) + 1;
			m_model_indices[tri * 3 + 1] = (i - 1) * seg2 + j;
			m_model_indices[tri * 3 + 2] = i * seg2 + (j % seg2) + 1;
			++tri;
		}
	}

	for (uint32_t j = 1; j < seg2 + 1; ++j)
	{
		m_model_indices[tri * 3 + 0] = (seg1 - 2) * seg2 + (j % seg2) + 1;
		m_model_indices[tri * 3 + 1] = (seg1 - 2) * seg2 + j;
		m_model_indices[tri * 3 + 2] = v_count - 1;
		++tri;
	}

}



struct cube_t : public model_base_t {

	cube_t() {

		m_model_vertex.resize(24);
		m_vertex_post.resize(24);
		m_model_indices.resize(36);

		m_model_vertex[0] = { vector_t(-1, -1, -1), vector_t(-1, 0, 0), texcoord_t(0, 0), vector_t(1, 0, 0) };
		m_model_vertex[1] = { vector_t(-1, 1, -1), vector_t(-1, 0, 0), texcoord_t(0, 1), vector_t(0, 1, 0) };
		m_model_vertex[2] = { vector_t(-1, 1, 1), vector_t(-1, 0, 0), texcoord_t(1, 1), vector_t(1, 0, 0) };
		m_model_vertex[3] = { vector_t(-1, -1, 1), vector_t(-1, 0, 0), texcoord_t(1, 0), vector_t(0, 0, 1) };
		m_model_vertex[4] = { vector_t(-1, -1, 1), vector_t(0, 0, 1), texcoord_t(1, 0), vector_t(0, 0, 1) };
		m_model_vertex[5] = { vector_t(-1, 1, 1), vector_t(0, 0, 1), texcoord_t(1, 1), vector_t(1, 0, 0) };
		m_model_vertex[6] = { vector_t(1, 1, 1), vector_t(0, 0, 1), texcoord_t(0, 1), vector_t(0, 1, 0) };
		m_model_vertex[7] = { vector_t(1, -1, 1), vector_t(0, 0, 1), texcoord_t(0, 0), vector_t(1, 0, 0) };
		m_model_vertex[8] = { vector_t(1, -1, 1), vector_t(1, 0, 0), texcoord_t(0, 0), vector_t(1, 0, 0) };
		m_model_vertex[9] = { vector_t(1, 1, 1), vector_t(1, 0, 0), texcoord_t(0, 1), vector_t(0, 1, 0) };
		m_model_vertex[10] = { vector_t(1, 1, -1), vector_t(1, 0, 0), texcoord_t(1, 1), vector_t(1, 0, 0) };
		m_model_vertex[11] = { vector_t(1, -1, -1), vector_t(1, 0, 0), texcoord_t(1, 0), vector_t(0, 0, 1) };
		m_model_vertex[12] = { vector_t(1, -1, -1), vector_t(0, 0, -1), texcoord_t(1, 0), vector_t(0, 0, 1) };
		m_model_vertex[13] = { vector_t(1, 1, -1), vector_t(0, 0, -1), texcoord_t(1, 1), vector_t(1, 0, 0) };
		m_model_vertex[14] = { vector_t(-1, 1, -1), vector_t(0, 0, -1), texcoord_t(0, 1), vector_t(0, 1, 0) };
		m_model_vertex[15] = { vector_t(-1, -1, -1), vector_t(0, 0, -1), texcoord_t(0, 0), vector_t(1, 0, 0) };
		m_model_vertex[16] = { vector_t(-1, 1, 1), vector_t(0, 1, 0), texcoord_t(1, 1), vector_t(1, 0, 0) };
		m_model_vertex[17] = { vector_t(-1, 1, -1), vector_t(0, 1, 0), texcoord_t(1, 0), vector_t(0, 1, 0) };
		m_model_vertex[18] = { vector_t(1, 1, -1), vector_t(0, 1, 0), texcoord_t(0, 0), vector_t(1, 0, 0) };
		m_model_vertex[19] = { vector_t(1, 1, 1), vector_t(0, 1, 0), texcoord_t(0, 1), vector_t(0, 1, 0) };
		m_model_vertex[20] = { vector_t(-1, -1, 1), vector_t(0, -1, 0), texcoord_t(1, 0), vector_t(0, 0, 1) };
		m_model_vertex[21] = { vector_t(1, -1, 1), vector_t(0, -1, 0), texcoord_t(0, 0), vector_t(1, 0, 0) };
		m_model_vertex[22] = { vector_t(1, -1, -1), vector_t(0, -1, 0), texcoord_t(0, 1), vector_t(0, 0, 1) };
		m_model_vertex[23] = { vector_t(-1, -1, -1), vector_t(0, -1, 0), texcoord_t(1, 1), vector_t(1, 0, 0) };

		int indexs[] = {
			0, 1, 2, 2, 3, 0,
			4, 5, 6, 6, 7, 4,
			8, 9, 10, 10, 11, 8,
			12, 13, 14, 14, 15, 12,
			16, 17, 18, 18, 19, 16,
			20, 21, 22, 22, 23, 20,
		};
		for (int i = 0; i < 36; ++i)
		{
			m_model_indices[i] = indexs[i];
		}
	};

};

#endif