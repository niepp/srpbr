#ifndef __PRE_COMPUTE_H__
#define __PRE_COMPUTE_H__

#include <thread>
#include <future>
#include <vector>
#include <string>


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
//// ----------------------------------------------------------------------------
//void main()
//{
//	vec2 integratedBRDF = IntegrateBRDF(TexCoords.x, TexCoords.y);
//	FragColor = integratedBRDF;
//}



cube_texture_t* calc_irradiance_map(const cube_texture_t *src_cube_tex)
{
	int w = src_cube_tex->size();
	cube_texture_t* dst_cube_tex = new cube_texture_t(w);

	auto face_filter_func = [&](texture2d_t &face, int face_id) {
		for (int i = 0; i < w; ++i)
		{
			for (int j = 0; j < w; ++j)
			{
				float u = 1.0f * j / (w - 1);
				float v = 1.0f * i / (w - 1);
				vector3_t n;
				cube_uv_to_direction(face_id, u, v, n);

				// [ right ]
				// | up    |
				// [ n     ]
				// makeup a rotation matrix
				vector3_t up = fabs(n.y) < 0.999f ? vector3_t(0.0f, 1.0f, 0.0f) : vector3_t(0.0f, 0.0f, 1.0f);
				vector3_t right = cross(up, n);
				up = cross(n, right);

				vector4_t irradiance(0, 0, 0, 0);
				float sample_step = 0.125f;
				int sample_num = 0;
				for (float phi = 0.0f; phi < 2.0f * cPI; phi += sample_step)
				{
					for (float theta = 0.0f; theta < 0.5f * cPI; theta += sample_step)
					{
						// spherical to cartesian (in tangent space)
						vector3_t tangent_dir(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
						// tangent space to world
						vector3_t sample_dir = right * tangent_dir.x + up * tangent_dir.y + n * tangent_dir.z;
						sample_dir.normalize();
						vector4_t color = src_cube_tex->sample(sample_dir);

						irradiance += color * sin(theta) * cos(theta);
						++sample_num;
					}
				}
				irradiance = irradiance * (cPI / sample_num);
				face.write_at(j, i, irradiance);
			}
		}
	};

	std::vector<std::future<void>> futures;
	for (int face_id = 0; face_id < 6; ++face_id) {
		texture2d_t& face = dst_cube_tex->get_face(face_id);
		futures.push_back(std::move(std::async(std::launch::async, [&face_filter_func, &face, face_id]() {
			face_filter_func(face, face_id);
		})));
	}

	for (auto& future : futures) {
		future.get();
	}

	return dst_cube_tex;

}

void generate_irradiance_map(const std::string& src_texpath, const std::string& dst_texpath)
{
	cube_texture_t src_cube;
	src_cube.load_tex(src_texpath, ".png");
	cube_texture_t *dst_cube = calc_irradiance_map(&src_cube);
	dst_cube->save_tex(dst_texpath);	
}

#endif // __PRE_COMPUTE_H__