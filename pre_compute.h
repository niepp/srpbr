#ifndef __PRE_COMPUTE_H__
#define __PRE_COMPUTE_H__

#include <thread>
#include <future>
#include <vector>
#include <string>
#include <iomanip>

float radical_inverse(unsigned int bits) {
	bits = (bits << 16u) | (bits >> 16u);
	bits = ((bits & 0x55555555u) << 1u) | ((bits & 0xAAAAAAAAu) >> 1u);
	bits = ((bits & 0x33333333u) << 2u) | ((bits & 0xCCCCCCCCu) >> 2u);
	bits = ((bits & 0x0F0F0F0Fu) << 4u) | ((bits & 0xF0F0F0F0u) >> 4u);
	bits = ((bits & 0x00FF00FFu) << 8u) | ((bits & 0xFF00FF00u) >> 8u);
	return float(bits) * 2.3283064365386963e-10f; // / 0x100000000
}

void hammersley(unsigned int i, unsigned int N, float &e1, float &e2)
{
	e1 = 1.0f * i / N;
	e2 = radical_inverse(i);
}

vector3_t hemisphere_sample_uniform(float u, float v)
{
	float phi = v * 2.0f * cPI;
	float cosTheta = 1.0f - u;
	float sinTheta = sqrt(1.0f - cosTheta * cosTheta);
	return vector3_t(cos(phi) * sinTheta, sin(phi) * sinTheta, cosTheta);
}


// [ tan_x ]
// | tan_y |
// [ n     ]
// makeup a rotation matrix
void makeup_rot(const vector3_t& n, vector3_t& tan_x, vector3_t& tan_y)
{
	vector3_t up = fabs(n.z) < 0.999f ? vector3_t(0.0f, 0.0f, 1.0f) : vector3_t(1.0f, 0.0f, 0.0f);
	tan_x = cross(up, n);
	tan_x.normalize();
	tan_y = cross(n, tan_x);
	tan_y.normalize();
}

void generate_irradiance_map(const std::string& src_texpath, const std::string& dst_texpath)
{
	cube_texture_t src_cube_tex;
	src_cube_tex.load_tex(src_texpath, false);

	int siz = src_cube_tex.size();
	cube_texture_t dst_cube_tex(siz);

	auto face_filter_func = [&](texture2d_t& face, int face_id) {
		for (int i = 0; i < siz; ++i)
		{
			for (int j = 0; j < siz; ++j)
			{
				texcoord_t tc;
				tc.u = 1.0f * i / (siz - 1);
				tc.v = 1.0f * j / (siz - 1);
				vector3_t n;
				cube_uv_to_direction(face_id, tc, n);

				vector3_t tan_x, tan_y;
				makeup_rot(n, tan_x, tan_y);

				vector3_t irradiance(0, 0, 0);

				const int sample_num = 8192;
				for (int k = 0; k < sample_num; ++k)
				{
					float e1, e2;
					hammersley(k, sample_num, e1, e2);

					// h is hemisphere direction in tangent space
					vector3_t h = hemisphere_sample_uniform(e1, e2);

					// convert h to world space
					vector3_t l = tan_x * h.x + tan_y * h.y + n * h.z;
					l.normalize();
					vector3_t color = src_cube_tex.sample(l).to_vec3();

					float NoL = max(dot(n, l), 0);
					irradiance += color * NoL;
					
				}
				irradiance = irradiance * 2.0f / float(sample_num);
				face.write_at(i, j, irradiance);
			}
		}
	};

	std::vector<std::future<void>> futures;
	for (int face_id = 0; face_id < 6; ++face_id) {
		texture2d_t& face = dst_cube_tex.get_face(face_id);
		futures.push_back(std::move(std::async(std::launch::async, [&face_filter_func, &face, face_id]() {
			face_filter_func(face, face_id);
		})));
	}

	for (auto& future : futures) {
		future.get();
	}

	dst_cube_tex.save_tex(dst_texpath, false);

}

vector3_t importance_sample_GGX(float e1, float e2, float roughness)
{
	float a = roughness * roughness;
	float a2 = a * a;

	float phi = 2 * cPI * e1;
	float cos_theta = sqrt((1 - e2) / (1 + (a2 - 1) * e2));
	float sin_theta = sqrt(1 - cos_theta* cos_theta);

	vector3_t h;
	h.x = sin_theta * cos(phi);
	h.y = sin_theta * sin(phi);
	h.z = cos_theta;

	return vector3_t(h.x, h.y, h.z);

}

// reference to "Real Shading in Unreal Engine 4. by Brian Karis, Epic Games"
void generate_prefilter_envmap(const std::string& src_texpath, const std::string& dst_texpath)
{
	cube_texture_t src_cube_tex;
	src_cube_tex.load_tex(src_texpath, false);

	auto face_filter_func = [&](float roughness, int siz, texture2d_t& face, int face_id) {
		for (int i = 0; i < siz; ++i)
		{
			for (int j = 0; j < siz; ++j)
			{
				texcoord_t tc;
				tc.u = 1.0f * i / (siz - 1);
				tc.v = 1.0f * j / (siz - 1);
				vector3_t n;
				cube_uv_to_direction(face_id, tc, n);

				vector3_t tan_x, tan_y;
				makeup_rot(n, tan_x, tan_y);

				vector3_t v = n; // specular is dependent on view direction, so assume v == n == r for precompute!

				vector3_t filtered_color(0, 0, 0);
				float weight = 0;

				const int sample_num = 1024;
				for (int k = 0; k < sample_num; ++k)
				{
					float e1, e2;
					hammersley(k, sample_num, e1, e2);

					vector3_t s = importance_sample_GGX(e1, e2, roughness);

					// convert normal direction s to world space
					vector3_t h = tan_x * s.x + tan_y * s.y + n * s.z;
					h.normalize();

					vector3_t l = reflect(h, v);
					l.normalize();

					float NoL = max(dot(n, l), 0.0f);
					if (NoL > 0)
					{
						vector3_t color = src_cube_tex.sample(l).to_vec3();
						filtered_color += color * NoL; // here the weight NoL is not present in Equation
						weight += NoL;
					}
				}

				filtered_color /= max(weight, 0.001f);
				face.write_at(i, j, filtered_color);
			}
		}
	};

	int max_siz = src_cube_tex.size();
	int mip_num = (int)ceil(std::log2(max_siz)) + 1;
	for (int i = 0; i < mip_num; ++i)
	{
		int siz = (max_siz >> i);
		if (siz <= 0) {
			break;
		}
		cube_texture_t dst_cube_tex(siz);
		float roughness = 1.0f * i / (mip_num - 1);

		std::vector<std::future<void>> futures;

		for (int face_id = 0; face_id < 6; ++face_id) {
			texture2d_t& face = dst_cube_tex.get_face(face_id);
			futures.push_back(std::move(std::async(std::launch::async, [&face_filter_func, roughness, siz, &face, face_id]() {
				face_filter_func(roughness, siz, face, face_id);
			})));
		}
		for (auto& future : futures) {
			future.get();
		}

		const char* ext = ::strrchr(src_texpath.c_str(), '.');
		dst_cube_tex.save_tex(dst_texpath + "_mip_" + std::to_string(i) + ext, false);

	}

}

float schlickGGX_geometry_ibl(float n_dot_v, float roughness)
{
	float k = roughness * roughness / 2.0f;
	return n_dot_v / (n_dot_v * (1 - k) + k);
}

float geometry_smith_ibl(float n_dot_v, float n_dot_l, float roughness)
{
	float g1 = schlickGGX_geometry_ibl(n_dot_v, roughness);
	float g2 = schlickGGX_geometry_ibl(n_dot_l, roughness);
	return g1 * g2;
}

vector3_t IntegrateBRDF(float NdotV, float roughness)
{
	vector3_t v;
	v.x = sqrt(1.0f - NdotV * NdotV);
	v.y = 0.0f;
	v.z = NdotV;

	float A = 0.0f;
	float B = 0.0f;

	vector3_t n(0.0f, 0.0f, 1.0f);
	vector3_t tan_x, tan_y;
	makeup_rot(n, tan_x, tan_y);

	const int SAMPLE_COUNT = 1024;
	for (int i = 0u; i < SAMPLE_COUNT; ++i)
	{
		float e1, e2;
		hammersley(i, SAMPLE_COUNT, e1, e2);
		vector3_t s = importance_sample_GGX(e1, e2, roughness);

		// convert h to world space
		vector3_t h = tan_x * s.x + tan_y * s.y + n * s.z;
		h.normalize();

		vector3_t l = reflect(h, v);
		l.normalize();

		float NoL = max(l.z, 0.0f);
		float NoV = max(v.z, 0.0f);
		float NoH = max(h.z, 0.0f);
		float VoH = max(dot(v, h), 0.0f);

		if (NoL > 0.0)
		{
			float G = geometry_smith_ibl(NoV, NoL, roughness);
			float G_Vis = (G * VoH) / (NoH * NoV);
			float Fc = pow(1.0f - VoH, 5.0f);
			A += (1.0f - Fc) * G_Vis;
			B += Fc * G_Vis;
		}
	}
	A /= float(SAMPLE_COUNT);
	B /= float(SAMPLE_COUNT);
	return vector3_t(A, B, 0.0f);
}

void generate_BRDF_LUT(const std::string& dst_texpath, int siz = 256)
{
	texture2d_t lut;
	lut.init(siz, siz);
	for (int i = 0; i < siz; ++i)
	{
		float NoV = 1.0f * i / (siz - 1.0f);
		for (int j = 0; j < siz; ++j)
		{
			float roughness = 1.0f * j / (siz - 1.0f);
			vector3_t c = IntegrateBRDF(NoV, roughness);
			lut.write_at(i, j, c);
		}
	}
	lut.save_tex(dst_texpath.c_str(), false);
}

#endif // __PRE_COMPUTE_H__
