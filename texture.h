#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

/* standard conversion from rgbe to float pixels */
/* note: Ward uses ldexp(rgb, e-(128+8)) */
void rgbe2float(float& red, float& green, float& blue, unsigned char rgbe[4])
{
	if (rgbe[3]) {   /*nonzero pixel*/
		float f = ldexp(1.0f, rgbe[3] - (int)(128 + 8));
		red = rgbe[0] * f;
		green = rgbe[1] * f;
		blue = rgbe[2] * f;
	}
	else {
		red = green = blue = 0.0f;
	}
}

vector4_t* load_tex_impl(const char* tex_path, int& width, int& height, bool is_srgb, bool vertical_flip = true)
{
	int components = 0;
	stbi_us* st_img = stbi_load_16(tex_path, &width, &height, &components, STBI_rgb_alpha);
	if (st_img == nullptr) {
		// failed to load the file.
		printf("Failed to load image %s\nReason: %s\n", tex_path, stbi_failure_reason());
		return nullptr;
	}

	vector4_t* texture = new vector4_t[width * height];
	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{
			int src_idx = i * width + j;
			int dst_idx = (vertical_flip ? height - 1 - i : i) * width + j;
			vector4_t c;
			c.r = st_img[src_idx * 4 + 0] * cRevt65535;
			c.g = st_img[src_idx * 4 + 1] * cRevt65535;
			c.b = st_img[src_idx * 4 + 2] * cRevt65535;
			c.a = st_img[src_idx * 4 + 3] * cRevt65535;
			c = is_srgb ? srgb_to_linear(c) : c;
			texture[dst_idx] = c;
		}
	}
	stbi_image_free(st_img);
	return texture;
}

void save_tex_impl(const char* tex_path, const int& width, const int& height, unsigned int* data)
{
	stbi_write_png(tex_path, width, height, 4, data, width * 4);
}

vector4_t sample_bilinear(vector4_t* texture, int width, int height, const texcoord_t& texcoord)
{
	float u = clamp(texcoord.u, 0.0f, 1.0f) * (width - 1);
	float v = clamp(texcoord.v, 0.0f, 1.0f) * (height - 1);
	int x0 = (int)(u);
	int y0 = (int)(v);
	int x1 = clamp(x0 + 1, 0, width - 1);
	int y1 = clamp(y0 + 1, 0, height - 1);
	const vector4_t& t00 = texture[y0 * width + x0];
	const vector4_t& t01 = texture[y1 * width + x0];
	const vector4_t& t10 = texture[y0 * width + x1];
	const vector4_t& t11 = texture[y1 * width + x1];

	float u_weight = u - x0;
	float v_weight = v - y0;

	// bilinear interpolation
	vector4_t tu0 = lerp(t00, t10, u_weight);
	vector4_t tu1 = lerp(t01, t11, u_weight);
	return lerp(tu0, tu1, v_weight);
}

/*

		！！！！！！
	   |  +z  |
	   |face2 |
 ！！！！！！ ！！！！！！ ！！！！！！ ！！！！！！
|  -x  |  -y  |  +x  |  +y  |
|face1 |face4 |face0 |face5 |
 ！！！！！！ ！！！！！！ ！！！！！！ ！！！！！！
	   |  -z  |
	   |face3 |
		！！！！！！

		^z
		|
		|！！！！！！
		|  +x  |
		|face0 |
		！！！！！！！！！！>y

		^z
		|
		|！！！！！！
		|  -x  |
		|face1 |
		！！！！！！！！！！>-y

		^y
		|
		|！！！！！！
		|  +z  |
		|face2 |
		！！！！！！！！！！>+x

		^-y
		|
		|！！！！！！
		|  -z  |
		|face3 |
		！！！！！！！！！！>+x

		^z
		|
		|！！！！！！
		|  -y  |
		|face4 |
		！！！！！！！！！！>+x

		^z
		|
		|！！！！！！
		|  +y  |
		|face5 |
		！！！！！！！！！！>-x

*/
int direction_to_cubeuv(const vector3_t& dir, texcoord_t& tc)
{
	int face_index = 0;
	float absX = abs(dir.x);
	float absY = abs(dir.y);
	float absZ = abs(dir.z);
	float ma;
	if (absZ >= absX && absZ >= absY)
	{
		face_index = dir.z > 0.0f ? 2 : 3;
		tc.set(dir.x, dir.z > 0.0f ? dir.y : -dir.y);
		ma = absZ;
	}
	else if (absY >= absX)
	{
		face_index = dir.y > 0.0f ? 5 : 4;
		tc.set(dir.y > 0.0f ? -dir.x : dir.x, dir.z);
		ma = absY;
	}
	else
	{
		face_index = dir.x > 0.0f ? 0 : 1;
		tc.set(dir.x > 0.0f ? dir.y : -dir.y, dir.z);
		ma = absX;
	}
	tc.u = (tc.u / ma + 1.0f) * 0.5f;
	tc.v = (tc.v / ma + 1.0f) * 0.5f;
	return face_index;
}



void cube_uv_to_direction(int index, const texcoord_t& tc, vector3_t& dir)
{
	// convert range [0, 1] to [-1, 1]
	float u = 2.0f * tc.u - 1.0f;
	float v = 2.0f * tc.v - 1.0f;
	switch (index)
	{
	case 0: dir.set(1.0f, u, v); break;		// +X
	case 1: dir.set(-1.0f, -u, v); break;	// -X
	case 2: dir.set(u, v, 1.0f); break;		// +Z
	case 3: dir.set(u, -v, -1.0f); break;	// -Z
	case 4: dir.set(u, -1.0f, v); break;	// -Y
	case 5: dir.set(-u, 1.0f, v); break;	// +Y
	}
	dir.normalize();
}

class texture2d_t
{
	vector4_t* texture;
	int width, height;
public:
	texture2d_t() :
		texture(nullptr), width(0), height(0)
	{ }

	~texture2d_t()
	{
		delete[] texture;
		texture = nullptr;
	}

	int tex_width() const
	{
		return width;
	}

	int tex_height() const
	{
		return width;
	}

	void init(int w, int h)
	{
		if (width != w || height != h) {
			width = w;
			height = h;
			if (texture != nullptr) {
				delete texture;
			}
			texture = new vector4_t[width * height];
		}
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				int idx = i * width + j;
				texture[idx].set(0.0f, 0.0f, 0.0f, 1.0f);
			}
		}
	}

	void load_tex(const char* tex_path, bool is_srgb, bool vertical_flip = true)
	{
		if (texture) {
			delete texture;
		}
		texture = load_tex_impl(tex_path, width, height, is_srgb, vertical_flip);
	}

	void save_tex(const char* tex_path, bool is_srgb, bool vertical_flip = true)
	{
		unsigned int* data = new unsigned int[width * height];
		for (int i = 0; i < height; ++i) {
			for (int j = 0; j < width; ++j) {
				int src_idx = i * width + j;
				int dst_idx = (vertical_flip ? height - 1 - i : i) * width + j;
				vector4_t color = is_srgb ? gamma_correction(texture[src_idx]) : texture[src_idx];
				data[dst_idx] = makefour(color);
			}
		}
		save_tex_impl(tex_path, width, height, data);
		delete[] data;
	}

	void write_at(int x, int y, const vector4_t& color)
	{
		assert(x >= 0 && x < width);
		assert(y >= 0 && y < height);
		texture[y * width + x] = color;
	}

	vector4_t read_at(int x, int y) const
	{
		assert(x >= 0 && x < width);
		assert(y >= 0 && y < height);
		return texture[y * width + x];
	}

	vector4_t sample(const texcoord_t& texcoord) const
	{
		return sample_bilinear(texture, width, height, texcoord);
	}

};

class cube_texture_t
{
	const int base_xy[6][2] = {
		{2, 1}, {0, 1}, {1, 2},
		{1, 0}, {1, 1}, {3, 1}
	};
	texture2d_t faces[6];
public:
	cube_texture_t()
	{
	}

	cube_texture_t(int siz)
	{
		for (int i = 0; i < 6; ++i) {
			faces[i].init(siz, siz);
		}
	}

	int size() const
	{
		return faces[0].tex_width();
	}

	const texture2d_t& get_face(int face_id) const
	{
		return faces[face_id];
	}

	texture2d_t& get_face(int face_id)
	{
		return faces[face_id];
	}

	void load_tex(const std::string& tex_path, const std::string& ext_name, bool is_srgb, bool vertical_flip = true)
	{
		const std::string face_tags[6] = { "px", "nx", "pz", "nz", "ny", "py" };
		for (int i = 0; i < 6; ++i) {
			std::string path = tex_path + "_" + face_tags[i] + ext_name;
			faces[i].load_tex(path.c_str(), is_srgb, vertical_flip);
		}
	}

	void load_tex(const std::string& tex_path, bool is_srgb, bool vertical_flip = true)
	{
		int w = 0, h = 0;
		vector4_t *texture = load_tex_impl(tex_path.c_str(), w, h, is_srgb, vertical_flip);
		assert(w * 3 == h * 4);
		int siz = w / 4;
		assert(w == siz * 4);

		for (int k = 0; k < 6; ++k) {
			faces[k].init(siz, siz);
			int bx = base_xy[k][0] * siz;
			int by = base_xy[k][1] * siz;
			for (int i = 0; i < siz; ++i) {
				for (int j = 0; j < siz; ++j) {
					int idx = (i + by) * w + (j + bx);
					faces[k].write_at(j, i, texture[idx]);
				}
			}
		}
		delete[] texture;
	}

	void save_tex(const std::string& tex_path, bool is_srgb)
	{
		int siz = this->size();
		int w = siz * 4;
		int h = siz * 3;
		unsigned int* data = new unsigned int[w * h];
		for (int i = 0; i < w *h; ++i) {
			data[i] = 0xFFFFFFFF;
		}

		for (int k = 0; k < 6; ++k)	{
			int bx = base_xy[k][0] * siz;
			int by = base_xy[k][1] * siz;
			for (int i = 0; i < siz; ++i) {
				for (int j = 0; j < siz; ++j) {
					int idx = (h - 1 - (i + by)) * w + (j + bx);
					vector4_t color = faces[k].read_at(j, i);
					color = is_srgb ? gamma_correction(color) : color;
					data[idx] = makefour(color);
				}
			}
		}
		save_tex_impl(tex_path.c_str(), w, h, data);
		delete[] data;
	}

	vector4_t sample(const vector3_t& dir) const
	{
		assert(appro_equal(dir.length(), 1.0f));
		texcoord_t tc;
		int face_id = direction_to_cubeuv(dir, tc);
		return faces[face_id].sample(tc);
	}

};


#endif //__TEXTURE_H__
