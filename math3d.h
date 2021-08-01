#ifndef __MATH3D_H__
#define __MATH3D_H__

const float cEpslion = 1e-6f;
const float cPI = 3.1415926f;
const float cRevt255 = 1.0f / 255.0f;

struct texcoord_t
{
	float u, v;
	texcoord_t(float _u = 0, float _v = 0) :
		u(_u), v(_v)
	{ }
	void set(float _u = 0, float _v = 0)
	{
		u = _u;
		v = _v;	
	}
};

struct vector3_t
{
	union
	{
		float vec[3];
		struct
		{
			float x, y, z;
		};
	};

	vector3_t() :
		x(0.0f), y(0.0f), z(0.0f)
	{ }

	vector3_t(float _x, float _y, float _z) :
		x(_x), y(_y), z(_z)
	{ }

	vector3_t& operator =(const vector3_t& vec)
	{
		x = vec.x;
		y = vec.y;
		z = vec.z;
		return *this;
	}
	
	void set(float _x = 0, float _y = 0, float _z = 0)
	{
		x = _x;
		y = _y;
		z = _z;
	}

	vector3_t operator -() const
	{
		return vector3_t(-x, -y, -z);
	}

	vector3_t& operator +=(const vector3_t& vec)
	{
		x += vec.x;
		y += vec.y;
		z += vec.z;
		return *this;
	}

	vector3_t& operator -=(const vector3_t& vec)
	{
		x -= vec.x;
		y -= vec.y;
		z -= vec.z;
		return *this;
	}

	vector3_t operator +(const vector3_t& vec) const
	{
		return vector3_t(x + vec.x, y + vec.y, z + vec.z);
	}

	vector3_t operator -(const vector3_t& vec) const
	{
		return vector3_t(x - vec.x, y - vec.y, z - vec.z);
	}

	vector3_t& operator *=(const float& scale)
	{
		x *= scale;
		y *= scale;
		z *= scale;
		return *this;
	}

	vector3_t& operator /=(const float& scale)
	{
		const float reci = 1.0f / scale;
		x *= reci;
		y *= reci;
		z *= reci;
		return *this;
	}

	vector3_t operator *(const float& scale) const
	{
		return vector3_t(x * scale, y * scale, z * scale);
	}

	vector3_t operator /(const float& scale) const
	{
		const float reci = 1.0f / scale;
		return vector3_t(x * reci, y * reci, z * reci);
	}

	vector3_t operator *(const vector3_t& vec)
	{
		return vector3_t(x * vec.x, y * vec.y, z * vec.z);
	}

	float length() const
	{
		float sq = x * x + y * y + z * z;
		return (float)sqrt(sq);
	}

	void normalize()
	{
		float len = length();
		float inv = cEpslion;
		if (abs(len) > cEpslion)
		{
			inv = 1.0f / len;
		}
		x *= inv;
		y *= inv;
		z *= inv;
	}
};

struct vector4_t
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

	vector4_t() :
		x(0.0f), y(0.0f), z(0.0f), w(1.0f)
	{ }

	vector4_t(float _x, float _y, float _z, float _w = 1.0f):
		x(_x), y(_y), z(_z), w(_w) 
	{ }

	vector4_t(const vector3_t& v3, float _w = 1.0f) :
		x(v3.x), y(v3.y), z(v3.z), w(_w)
	{ }

	vector4_t& operator =(const vector4_t& vec)
	{
		x = vec.x;
		y = vec.y;
		z = vec.z;
		w = vec.w;
		return *this;
	}

	void set(float _x = 0, float _y = 0, float _z = 0, float _w = 1.0f)
	{
		x = _x;
		y = _y;
		z = _z;
		w = _w;
	}

	vector4_t operator -() const
	{
		return vector4_t(-x, -y, -z, -w);
	}

	vector4_t& operator +=(const vector4_t& vec)
	{
		x += vec.x;
		y += vec.y;
		z += vec.z;
		w += vec.w;
		return *this;
	}

	vector4_t& operator -=(const vector4_t& vec)
	{
		x -= vec.x;
		y -= vec.y;
		z -= vec.z;
		w -= vec.w;
		return *this;
	}

	vector4_t operator +(const vector4_t& vec) const
	{
		return vector4_t(x + vec.x, y + vec.y, z + vec.z, w + vec.w);
	}

	vector4_t operator -(const vector4_t& vec) const
	{
		return vector4_t(x - vec.x, y - vec.y, z - vec.z, w - vec.w);
	}

	vector4_t& operator *=(const float& scale)
	{
		x *= scale;
		y *= scale;
		z *= scale;
		w *= scale;
		return *this;
	}

	vector4_t& operator /=(const float& scale)
	{
		const float reci = 1.0f / scale;
		x *= reci;
		y *= reci;
		z *= reci;
		w *= reci;
		return *this;
	}

	vector4_t operator *(const float& scale) const
	{
		return vector4_t(x * scale, y * scale, z * scale, w * scale);
	}

	vector4_t operator /(const float& scale) const
	{
		const float reci = 1.0f / scale;
		return vector4_t(x * reci, y * reci, z * reci, w * reci);
	}

	vector4_t operator *(const vector4_t& vec)
	{
		return vector4_t(x * vec.x, y * vec.y, z * vec.z, w * vec.w);
	}

	float length() const
	{
		float sq = x * x + y * y + z * z + w * w;
		return (float)sqrt(sq);
	}

	float length3d() const
	{
		float sq = x * x + y * y + z * z;
		return (float)sqrt(sq);
	}

	void normalize()
	{
		float len = length();
		float inv = cEpslion;
		if (abs(len) > cEpslion)
		{
			inv = 1.0f / len;
		}
		x *= inv;
		y *= inv;
		z *= inv;
		w *= inv;
	}

	vector3_t to_vec3() const
	{
		return vector3_t(x, y, z);
	}

};

inline float dot(const vector3_t& a, const vector3_t& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline vector3_t cross(const vector3_t& a, const vector3_t& b)
{
	vector3_t c;
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	return c;
}

inline float dot(const vector4_t& a, const vector4_t& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline vector4_t cross(const vector4_t& a, const vector4_t& b)
{
	vector4_t c;
	c.x = a.y * b.z - a.z * b.y;
	c.y = a.z * b.x - a.x * b.z;
	c.z = a.x * b.y - a.y * b.x;
	c.w = 1.0f;
	return c;
}

struct matrix_t
{
	float m[4][4];

	matrix_t()
	{
		m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
		m[0][1] = m[0][2] = m[0][3] = 0.0f;
		m[1][0] = m[1][2] = m[1][3] = 0.0f;
		m[2][0] = m[2][1] = m[2][3] = 0.0f;
		m[3][0] = m[3][1] = m[3][2] = 0.0f;
	}

	vector3_t operator *(const vector3_t& vec) const
	{
		vector3_t ret;
		float X = vec.x, Y = vec.y, Z = vec.z;
		ret.x = X * m[0][0] + Y * m[1][0] + Z * m[2][0];
		ret.y = X * m[0][1] + Y * m[1][1] + Z * m[2][1];
		ret.z = X * m[0][2] + Y * m[1][2] + Z * m[2][2];
		return ret;
	}

	vector4_t operator *(const vector4_t& vec) const
	{
		vector4_t ret;
		float X = vec.x, Y = vec.y, Z = vec.z, W = vec.w;
		ret.x = X * m[0][0] + Y * m[1][0] + Z * m[2][0] + W * m[3][0];
		ret.y = X * m[0][1] + Y * m[1][1] + Z * m[2][1] + W * m[3][1];
		ret.z = X * m[0][2] + Y * m[1][2] + Z * m[2][2] + W * m[3][2];
		ret.w = X * m[0][3] + Y * m[1][3] + Z * m[2][3] + W * m[3][3];
		return ret;
	}

	void set_identity()
	{
		m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
		m[0][1] = m[0][2] = m[0][3] = 0.0f;
		m[1][0] = m[1][2] = m[1][3] = 0.0f;
		m[2][0] = m[2][1] = m[2][3] = 0.0f;
		m[3][0] = m[3][1] = m[3][2] = 0.0f;
	}

	void set_translate(float x, float y, float z)
	{
		set_identity();
		m[3][0] = x;
		m[3][1] = y;
		m[3][2] = z;
	}

	void set_scale(float x, float y, float z)
	{
		set_identity();
		m[0][0] = x;
		m[1][1] = y;
		m[2][2] = z;
	}

	void set_rotate(float x, float y, float z, float theta)
	{
		float qsin = (float)sin(theta * 0.5f);
		float qcos = (float)cos(theta * 0.5f);
		vector3_t vec(x, y, z);
		float w = qcos;
		vec.normalize();
		x = vec.x * qsin;
		y = vec.y * qsin;
		z = vec.z * qsin;
		m[0][0] = 1 - 2 * y * y - 2 * z * z;
		m[1][0] = 2 * x * y - 2 * w * z;
		m[2][0] = 2 * x * z + 2 * w * y;
		m[0][1] = 2 * x * y + 2 * w * z;
		m[1][1] = 1 - 2 * x * x - 2 * z * z;
		m[2][1] = 2 * y * z - 2 * w * x;
		m[0][2] = 2 * x * z - 2 * w * y;
		m[1][2] = 2 * y * z + 2 * w * x;
		m[2][2] = 1 - 2 * x * x - 2 * y * y;
		m[0][3] = m[1][3] = m[2][3] = 0.0f;
		m[3][0] = m[3][1] = m[3][2] = 0.0f;
		m[3][3] = 1.0f;
	}

	// view matrix
	void set_lookat(const vector3_t& eye, const vector3_t& at, const vector3_t& up)
	{
		vector3_t xaxis, yaxis, zaxis;
		zaxis = at - eye;
		zaxis.normalize();

		xaxis = cross(up, zaxis);
		xaxis.normalize();

		yaxis = cross(zaxis, xaxis);

		m[0][0] = xaxis.x;
		m[1][0] = xaxis.y;
		m[2][0] = xaxis.z;
		m[3][0] = -dot(xaxis, eye);

		m[0][1] = yaxis.x;
		m[1][1] = yaxis.y;
		m[2][1] = yaxis.z;
		m[3][1] = -dot(yaxis, eye);

		m[0][2] = zaxis.x;
		m[1][2] = zaxis.y;
		m[2][2] = zaxis.z;
		m[3][2] = -dot(zaxis, eye);

		m[0][3] = m[1][3] = m[2][3] = 0.0f;
		m[3][3] = 1.0f;
	}

	// projection matrix (ref to D3DXMatrixPerspectiveFovLH)
	void set_perspective(float fovy, float aspect, float zn, float zf)
	{
		float fax = 1.0f / (float)tan(fovy * 0.5f);
		m[0][0] = (float)(fax / aspect);
		m[1][1] = (float)(fax);
		m[2][2] = zf / (zf - zn);
		m[3][2] = -zn * zf / (zf - zn);
		m[2][3] = 1;
		m[3][3] = 0;
	}

};


inline matrix_t mul(const matrix_t& a, const matrix_t& b)
{
	matrix_t c;
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			float s = 0;
			for (int k = 0; k < 4; ++k)
			{
				s += (a.m[i][k] * b.m[k][j]);
			}
			c.m[i][j] = s;
		}
	}
	return c;
}

struct model_vertex_t
{
	vector3_t pos;
	vector3_t nor;
	texcoord_t uv;
	vector4_t color;
};

struct interp_vertex_t
{
	vector3_t wpos; // world position
	vector4_t pos;  // screen position
	vector3_t nor;  // world normal
	texcoord_t uv;
	vector4_t color;
};

typedef std::vector<model_vertex_t> model_vertex_vec_t;
typedef std::vector<interp_vertex_t> interp_vertex_vec_t;
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

unsigned int makefour(const vector4_t& color)
{
	return to_color_int(color.r)
		| to_color_int(color.g) << 8
		| to_color_int(color.b) << 16
		| to_color_int(color.a) << 24;
}

void to_color(unsigned int cint, vector4_t& color)
{
	color.r = cRevt255 * (cint & 0xff);
	color.g = cRevt255 * ((cint >> 8) & 0xff);
	color.b = cRevt255 * ((cint >> 16) & 0xff);
	color.a = cRevt255 * ((cint >> 24) & 0xff);
}

unsigned int lerp(unsigned int x1, unsigned int x2, float t)
{
	return (unsigned int)(x1 * (1.0f - t) + x2 * t + 0.5f);
}

float lerp(float x1, float x2, float t)
{
	return x1 + (x2 - x1) * t;
}

texcoord_t lerp(const texcoord_t& a, const texcoord_t& b, float w)
{
	texcoord_t p;
	p.u = lerp(a.u, b.u, w);
	p.v = lerp(a.v, b.v, w);
	return p;
}

vector3_t lerp(const vector3_t& a, const vector3_t& b, float w)
{
	vector3_t p;
	p.x = lerp(a.x, b.x, w);
	p.y = lerp(a.y, b.y, w);
	p.z = lerp(a.z, b.z, w);
	return p;
}

vector4_t lerp(const vector4_t& a, const vector4_t& b, float w)
{
	vector4_t p;
	p.x = lerp(a.x, b.x, w);
	p.y = lerp(a.y, b.y, w);
	p.z = lerp(a.z, b.z, w);
	p.w = lerp(a.w, b.w, w);
	return p;
}

interp_vertex_t lerp( const interp_vertex_t& a, const interp_vertex_t& b, float w)
{
	interp_vertex_t p;
	p.wpos = lerp(a.wpos, b.wpos, w);
	p.pos = lerp(a.pos, b.pos, w);
	p.nor = lerp(a.nor, b.nor, w);
	p.uv = lerp(a.uv, b.uv, w);
	p.color = lerp(a.color, b.color, w);
	return p;
}


struct model_base_t 
{
	model_vertex_vec_t m_model_vertex;
	index_vec_t m_model_indices;
	interp_vertex_vec_t m_vertex_post;
};


/*
	unit sphere with center at (0,0,0) and radius of 1.0
*/
struct sphere_t : public model_base_t
{
	sphere_t(int seg)
	{
		gen_sphere(seg, seg * 2);
	}
private:
	void assign_vertex(int v_idx, const vector3_t& p);
	void gen_sphere(int seg1, int seg2);
};

inline void sphere_t::assign_vertex(int v_idx, const vector3_t& p)
{
	m_model_vertex[v_idx].pos.set(p.x, p.y, p.z);
	m_model_vertex[v_idx].nor.set(p.x, p.y, p.z);
	m_model_vertex[v_idx].uv.set(p.x * 0.5f + 0.5f, p.y * 0.5f + 0.5f);
	m_model_vertex[v_idx].color.set(1, 0, 0, 1);
}

inline void sphere_t::gen_sphere(int seg1, int seg2)
{
	int v_count = (seg1 - 1) * seg2 + 2;
	m_model_vertex.resize(v_count);
	m_vertex_post.resize(v_count);

	assign_vertex(0, vector3_t(0, 0, 1.0f));
	for (int i = 0; i < seg1 - 1; ++i)
	{
		float alpha = cPI * (i + 1) / seg1;
		float radius = sin(alpha);
		float z = cos(alpha);
		for (int j = 0; j < seg2; ++j)
		{
			float beta = 2.0f * cPI * j / seg2;
			float x = radius * cos(beta);
			float y = radius * sin(beta);
			int v_idx = i * seg2 + j + 1;
			assign_vertex(v_idx, vector3_t(x, y, z));
		}
	}

	assign_vertex(v_count - 1, vector3_t(0, 0, -1.0f));
	m_model_indices.resize((seg1 * 2 - 2) * 3 * seg2);

	int tri = 0;
	for (int j = 1; j < seg2 + 1; ++j)
	{
		m_model_indices[tri * 3 + 0] = 0;
		m_model_indices[tri * 3 + 1] = j;
		m_model_indices[tri * 3 + 2] = (j % seg2) + 1;
		++tri;
	}

	for (int i = 1; i < seg1 - 1; ++i)
	{
		for (int j = 1; j < seg2 + 1; ++j)
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

	for (int j = 1; j < seg2 + 1; ++j)
	{
		m_model_indices[tri * 3 + 0] = (seg1 - 2) * seg2 + j;
		m_model_indices[tri * 3 + 1] = (seg1 - 2) * seg2 + (j % seg2) + 1;
		m_model_indices[tri * 3 + 2] = v_count - 1;
		++tri;
	}

	for (int i = 0; i < tri; ++i)
	{
		std::swap(m_model_indices[i * 3 + 0], m_model_indices[i * 3 + 1]);
	}

}


struct cube_t : public model_base_t {

	cube_t() {

		m_model_vertex.resize(24);
		m_vertex_post.resize(24);
		m_model_indices.resize(36);

		m_model_vertex[0] = { vector3_t(-1, -1, -1), vector3_t(-1, 0, 0), texcoord_t(0, 0), vector4_t(1, 0, 0) };
		m_model_vertex[1] = { vector3_t(-1, 1, -1), vector3_t(-1, 0, 0), texcoord_t(0, 1), vector4_t(0, 1, 0) };
		m_model_vertex[2] = { vector3_t(-1, 1, 1), vector3_t(-1, 0, 0), texcoord_t(1, 1), vector4_t(1, 0, 0) };
		m_model_vertex[3] = { vector3_t(-1, -1, 1), vector3_t(-1, 0, 0), texcoord_t(1, 0), vector4_t(0, 0, 1) };
		m_model_vertex[4] = { vector3_t(-1, -1, 1), vector3_t(0, 0, 1), texcoord_t(1, 0), vector4_t(0, 0, 1) };
		m_model_vertex[5] = { vector3_t(-1, 1, 1), vector3_t(0, 0, 1), texcoord_t(1, 1), vector4_t(1, 0, 0) };
		m_model_vertex[6] = { vector3_t(1, 1, 1), vector3_t(0, 0, 1), texcoord_t(0, 1), vector4_t(0, 1, 0) };
		m_model_vertex[7] = { vector3_t(1, -1, 1), vector3_t(0, 0, 1), texcoord_t(0, 0), vector4_t(1, 0, 0) };
		m_model_vertex[8] = { vector3_t(1, -1, 1), vector3_t(1, 0, 0), texcoord_t(0, 0), vector4_t(1, 0, 0) };
		m_model_vertex[9] = { vector3_t(1, 1, 1), vector3_t(1, 0, 0), texcoord_t(0, 1), vector4_t(0, 1, 0) };
		m_model_vertex[10] = { vector3_t(1, 1, -1), vector3_t(1, 0, 0), texcoord_t(1, 1), vector4_t(1, 0, 0) };
		m_model_vertex[11] = { vector3_t(1, -1, -1), vector3_t(1, 0, 0), texcoord_t(1, 0), vector4_t(0, 0, 1) };
		m_model_vertex[12] = { vector3_t(1, -1, -1), vector3_t(0, 0, -1), texcoord_t(1, 0), vector4_t(0, 0, 1) };
		m_model_vertex[13] = { vector3_t(1, 1, -1), vector3_t(0, 0, -1), texcoord_t(1, 1), vector4_t(1, 0, 0) };
		m_model_vertex[14] = { vector3_t(-1, 1, -1), vector3_t(0, 0, -1), texcoord_t(0, 1), vector4_t(0, 1, 0) };
		m_model_vertex[15] = { vector3_t(-1, -1, -1), vector3_t(0, 0, -1), texcoord_t(0, 0), vector4_t(1, 0, 0) };
		m_model_vertex[16] = { vector3_t(-1, 1, 1), vector3_t(0, 1, 0), texcoord_t(1, 1), vector4_t(1, 0, 0) };
		m_model_vertex[17] = { vector3_t(-1, 1, -1), vector3_t(0, 1, 0), texcoord_t(1, 0), vector4_t(0, 1, 0) };
		m_model_vertex[18] = { vector3_t(1, 1, -1), vector3_t(0, 1, 0), texcoord_t(0, 0), vector4_t(1, 0, 0) };
		m_model_vertex[19] = { vector3_t(1, 1, 1), vector3_t(0, 1, 0), texcoord_t(0, 1), vector4_t(0, 1, 0) };
		m_model_vertex[20] = { vector3_t(-1, -1, 1), vector3_t(0, -1, 0), texcoord_t(1, 0), vector4_t(0, 0, 1) };
		m_model_vertex[21] = { vector3_t(1, -1, 1), vector3_t(0, -1, 0), texcoord_t(0, 0), vector4_t(1, 0, 0) };
		m_model_vertex[22] = { vector3_t(1, -1, -1), vector3_t(0, -1, 0), texcoord_t(0, 1), vector4_t(0, 0, 1) };
		m_model_vertex[23] = { vector3_t(-1, -1, -1), vector3_t(0, -1, 0), texcoord_t(1, 1), vector4_t(1, 0, 0) };

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

#endif //__MATH3D_H__
