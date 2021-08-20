#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

struct texture2d_t
{
	unsigned int *texture;
	int width, height;

	texture2d_t() :
		texture(nullptr), width(0), height(0)
	{ }

	void load_tex(const char* tex_path)
	{
		int components = 0;
		stbi_uc* st_img = stbi_load(tex_path, &width, &height, &components, STBI_rgb_alpha);
		if (st_img == nullptr) {
			// failed to load the file.
			printf("Failed to load image %s\nReason: %s\n", tex_path, stbi_failure_reason());
			return;
		}

		texture = new unsigned int[width * height];

		for (int i = 0; i < height; ++i)
		{
			for (int j = 0; j < width; ++j)
			{
				int pidx = (height - 1 - i) * width + j;
				vector4_t c;
				c.b = st_img[pidx * 4 + 0] / 255.0f;
				c.g = st_img[pidx * 4 + 1] / 255.0f;
				c.r = st_img[pidx * 4 + 2] / 255.0f;
				c.a = st_img[pidx * 4 + 3] / 255.0f;
				texture[pidx] = makefour(c);
			}
		}
		stbi_image_free(st_img);
	}

	vector4_t sample(const texcoord_t& texcoord)
	{
		float u = clamp(texcoord.u, 0.0f, 1.0f) * (width - 1);
		float v = clamp(texcoord.v, 0.0f, 1.0f) * (height - 1);
		int x0 = (int)(u);
		int y0 = (int)(v);
		int x1 = clamp(x0 + 1, 0, width - 1);
		int y1 = clamp(y0 + 1, 0, height - 1);
		unsigned int t00 = texture[y0 * width + x0];
		unsigned int t01 = texture[y1 * width + x0];
		unsigned int t10 = texture[y0 * width + x1];
		unsigned int t11 = texture[y1 * width + x1];

		float u_weight = u - x0;
		float v_weight = v - y0;

		vector4_t c00, c10, c01, c11;
		to_color(t00, c00);
		to_color(t10, c10);
		to_color(t01, c01);
		to_color(t11, c11);

		// bilinear interpolation
		vector4_t tu0 = lerp(c00, c10, u_weight);
		vector4_t tu1 = lerp(c01, c11, u_weight);
		return lerp(tu0, tu1, v_weight);
	}
};

#endif //__TEXTURE_H__
