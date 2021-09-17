#ifndef __MATH3D_H__
#define __MATH3D_H__

/*
	 z
	 ^
	 |
	 |
	 |
	 |
	 o--------------> x
    /
   /
 |/_
 y

the coordinate system same to UE4: left-hand z-up
*/

const float cEpslion = 1e-6f;
const float cPI = 3.1415926f;
const float cRevt255 = 1.0f / 255.0f;
const float cRevt65535 = 1.0f / 65535.0f;


inline bool is_valid(float f)
{
	return f == f;
}

inline float fast_inv_sqrt(float x)
{
	float y = x;
	float x2 = x * 0.5f;
	int i = *(int*)&y;  // evil floating point bit level hacking
	i = 0x5f3759df - (i >> 1); // what the fuck?
	y = *(float*)&i;
	y = y * (1.5f - (x2 * y * y)); // 1st iteration
	// y = y * (1.5f - (x2 * y * y)); // 2nd iteration, this can be removed
	return y;
}

inline double fast_inv_sqrt(double x)
{
	double y = x;
	double x2 = y * 0.5;
	int64_t i = *(int64_t*)&y;
	// The magic number is for doubles is from https://cs.uwaterloo.ca/~m32rober/rsqrt.pdf
	i = 0x5fe6eb50c7b537a9 - (i >> 1);
	y = *(double*)&i;
	y = y * (1.5 - (x2 * y * y));   // 1st iteration
	// y  = y * ( 1.5 - ( x2 * y * y ) ); // 2nd iteration, this can be removed
	return y;
}

inline float fast_exp(float x)
{
	int a = 185 * (int)x + 16249;
	a <<= 16;
	float f = *(reinterpret_cast<float*>(&a));
	return f;
}

inline double fast_exp(double x)
{
	double d;
	*(reinterpret_cast<int*>(&d) + 0) = 0;
	*(reinterpret_cast<int*>(&d) + 1) = static_cast<int>(1512775 * x + 1072632447);
	return d;
}



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

	static vector3_t one() {
		static vector3_t cOne(1.0f, 1.0f, 1.0f);
		return cOne;
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

	static vector4_t one() {
		static vector4_t cOne(1.0f, 1.0f, 1.0f, 1.0f);
		return cOne;
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

template <typename T>
T clamp(const T& x, const T& min, const T& max)
{
	return (x < min) ? min : ((x > max) ? max : x);
}

template <typename T>
void clamp(T& x, T min, T max)
{
	x = (x < min) ? min : ((x > max) ? max : x);
}

template <typename T>
bool appro_equal(T t, T c)
{
	return abs(t - c) < cEpslion;
}

uint8_t to_color_int(float c)
{
	int cint = (int)(c * 255.0f + 0.5);
	clamp(cint, 0, 255);
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

vector3_t minv(const vector3_t& a, const vector3_t& b)
{
	vector3_t p;
	p.x = min(a.x, b.x);
	p.y = min(a.y, b.y);
	p.z = min(a.z, b.z);
	return p;
}

vector3_t maxv(const vector3_t& a, const vector3_t& b)
{
	vector3_t p;
	p.x = max(a.x, b.x);
	p.y = max(a.y, b.y);
	p.z = max(a.z, b.z);
	return p;
}

vector4_t minv(const vector4_t& a, const vector4_t& b)
{
	vector4_t p;
	p.x = min(a.x, b.x);
	p.y = min(a.y, b.y);
	p.z = min(a.z, b.z);
	p.w = min(a.w, b.w);
	return p;
}

vector4_t maxv(const vector4_t& a, const vector4_t& b)
{
	vector4_t p;
	p.x = max(a.x, b.x);
	p.y = max(a.y, b.y);
	p.z = max(a.z, b.z);
	p.w = max(a.w, b.w);
	return p;
}

vector3_t reflect(const vector3_t& n, const vector3_t& l)
{
	vector3_t r = n * 2.0f * dot(n, l) - l;
	return r;
}


#endif //__MATH3D_H__
