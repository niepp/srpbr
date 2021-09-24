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

//vec2 IntegrateBRDF(float NdotV, float roughness)
//{
//	vec3 V;
//	V.x = sqrt(1.0 - NdotV * NdotV);
//	V.y = 0.0;
//	V.z = NdotV;
//
//	float A = 0.0;
//	float B = 0.0;
//
//	vec3 N = vec3(0.0, 0.0, 1.0);
//
//	const uint SAMPLE_COUNT = 1024u;
//	for (uint i = 0u; i < SAMPLE_COUNT; ++i)
//	{
//		vec2 Xi = Hammersley(i, SAMPLE_COUNT);
//		vec3 H = ImportanceSampleGGX(Xi, N, roughness);
//		vec3 L = normalize(2.0 * dot(V, H) * H - V);
//
//		float NdotL = max(L.z, 0.0);
//		float NdotH = max(H.z, 0.0);
//		float VdotH = max(dot(V, H), 0.0);
//
//		if (NdotL > 0.0)
//		{
//			float G = GeometrySmith(N, V, L, roughness);
//			//(DF/(4*NdotL*NdotV))/DPF
//			float G_Vis = (G * VdotH) / (NdotH * NdotV); //DG/pdf的感觉
//			float Fc = pow(1.0 - VdotH, 5.0);//是公式里的
//
//			A += (1.0 - Fc) * G_Vis;
//			B += Fc * G_Vis;
//		}
//	}
//	A /= float(SAMPLE_COUNT);
//	B /= float(SAMPLE_COUNT);
//	return vec2(A, B);
//}

void generate_irradiance_map(const std::string& src_texpath, const std::string& dst_texpath)
{
	cube_texture_t src_cube_tex;
	src_cube_tex.load_tex(src_texpath, ".png", true);
	
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

				// [ tan_x ]
				// | tan_y |
				// [ n     ]
				// makeup a rotation matrix
				vector3_t up = fabs(n.z) < 0.999f ? vector3_t(0.0f, 0.0f, 1.0f) : vector3_t(1.0f, 0.0f, 0.0f);
				vector3_t tan_x = cross(up, n);
				tan_x.normalize();
				vector3_t tan_y = cross(n, tan_x);
				tan_y.normalize();

				vector3_t irradiance(0, 0, 0);
				const int sample_num = 100;
				for (int u = 0; u < sample_num; ++u)
				{
					for (int v = 0; v < sample_num; ++v)
					{
						// h is hemisphere direction in tangent space
						vector3_t h = hemisphere_sample_uniform(1.0f * u / (sample_num - 1), 1.0f * v / (sample_num - 1));

						// convert h to world space
						vector3_t sample_dir = tan_x * h.x + tan_y * h.y + n * h.z;
						sample_dir.normalize();
						vector3_t color = src_cube_tex.sample(sample_dir).to_vec3();

						float n_dot_l = max(dot(sample_dir, n), 0);
						irradiance += color * n_dot_l;

					}
				}

				irradiance = irradiance * cPI / (sample_num * sample_num);
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

	dst_cube_tex.save_tex(dst_texpath, true);

}


vector4_t importance_sample_GGX(float e1, float e2, const vector3_t& n, float roughness)
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

	float d = (cos_theta * a2 - cos_theta) * cos_theta + 1.0f;
	float D = a2 / (cPI * d * d);
	float PDF = D * cos_theta;

	return vector4_t(h.x, h.y, h.z, PDF);

}

void generate_prefilter_envmap(const std::string& src_texpath, const std::string& dst_texpath)
{
	cube_texture_t src_cube_tex;
	src_cube_tex.load_tex(src_texpath, ".png", true);

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

				// [ tan_x ]
				// | tan_y |
				// [ n     ]
				// makeup a rotation matrix
				vector3_t up = fabs(n.z) < 0.999f ? vector3_t(0.0f, 0.0f, 1.0f) : vector3_t(1.0f, 0.0f, 0.0f);
				vector3_t tan_x = cross(up, n);
				tan_x.normalize();
				vector3_t tan_y = cross(n, tan_x);
				tan_y.normalize();

				vector3_t r = n;
				vector3_t v = r;

				vector3_t filtered_color(0, 0, 0);
				float weight = 0;

				const int sample_num = 1024;
				for (int i = 0; i < sample_num; ++i)
				{
					float e1, e2;
					hammersley(i, sample_num, e1, e2);

					vector4_t s = importance_sample_GGX(e1, e2, n, roughness);

					// convert h to world space
					vector3_t h = tan_x * s.x + tan_y * s.y + n * s.z;
					h.normalize();

					vector3_t l = h * 2.0f * dot(v, h) - v;
					l.normalize();

					float NoL = max(dot(n, l), 0.0f);
					if (NoL > 0)
					{
						vector3_t color = src_cube_tex.sample(l).to_vec3();
						filtered_color += color * NoL;
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

		dst_cube_tex.save_tex(dst_texpath + "_mip_" + std::to_string(i), true);

	}

}

#endif // __PRE_COMPUTE_H__
