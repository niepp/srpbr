#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


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

vector4_t* load_tex_impl(const char* tex_path, int& width, int& height)
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
			int pidx = (height - 1 - i) * width + j;
			vector4_t c;
			c.b = st_img[pidx * 4 + 0] * cRevt65535;
			c.g = st_img[pidx * 4 + 1] * cRevt65535;
			c.r = st_img[pidx * 4 + 2] * cRevt65535;
			c.a = st_img[pidx * 4 + 3] * cRevt65535;
			texture[pidx] = c;
		}
	}
	stbi_image_free(st_img);
	return texture;
}

vector4_t sample_bilinear_interpolation(vector4_t* texture, int width, int height, const texcoord_t& texcoord)
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

void cube_uv_to_direction(int index, float u, float v, vector3_t& dir)
{
	// convert range [0, 1] to [-1, 1]
	u = 2.0f * u - 1.0f;
	v = 2.0f * v - 1.0f;
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
		delete texture;
		texture = nullptr;
	}

	void load_tex(const char* tex_path)
	{
		texture = load_tex_impl(tex_path, width, height);
	}

	vector4_t sample(const texcoord_t& texcoord)
	{
		return sample_bilinear_interpolation(texture, width, height, texcoord);
	}

};

class cube_texture_t
{
	texture2d_t surfaces[6];
public:
	cube_texture_t()
	{
	}

	void load_tex(const std::string& tex_path, const std::string& ext_name)
	{
		const std::string faces[6] = { "px", "nx", "py", "ny", "pz", "nz" };
		for (int i = 0; i < 6; ++i) {
			std::string path = tex_path + "_" + faces[i] + ext_name;
			surfaces[i].load_tex(path.c_str());
		}
	}

	vector4_t sample(const vector3_t& dir)
	{
		assert(appro_equal(dir.length(), 1.0f));
		texcoord_t tc;
		int face_index = direction_to_cubeuv(dir, tc);
		return surfaces[face_index].sample(tc);
	}

};

#endif //__TEXTURE_H__
