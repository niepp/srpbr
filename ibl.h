#ifndef __IBL_H__
#define __IBL_H__

#include <vector>
#include <string>
#include <sstream>

class texture2d_t;
class cube_texture_t;

class ibl_t
{
	cube_texture_t *irradiance_map; // diffuse environment map
	std::vector<cube_texture_t*> prefilter_maps; // specular environment mipmaps
	texture2d_t *brdf_lut;			// brdf lookup texture
public:

	ibl_t() {
		irradiance_map = new cube_texture_t;
		brdf_lut = new texture2d_t;

		/* brdf lookup texture */
		std::string lut_path = "./resource/brdf_lut.png";
		brdf_lut->load_tex(lut_path.c_str(), false);
	}

	virtual ~ibl_t() {
		delete irradiance_map;
		irradiance_map = nullptr;

		for (auto m : prefilter_maps) {
			delete m;
		}
		prefilter_maps.clear();

		delete brdf_lut;
		brdf_lut = nullptr;
	}

	void load(const std::string& tex_dir_path)
	{
		if (!prefilter_maps.empty()) {
			for (auto m : prefilter_maps) {
				delete m;
			}
			prefilter_maps.clear();
		}

		irradiance_map->load_tex(tex_dir_path + "irradiance.png", false);
		for (int i = 0; i < 10; ++i)
		{
			std::ostringstream ss;
			ss << tex_dir_path;
			ss << "prefilter_mip_";
			ss << i;
			ss << ".png";
			cube_texture_t* pf_map = new cube_texture_t;
			pf_map->load_tex(ss.str(), false);
			prefilter_maps.push_back(pf_map);
		}
	}

	vector3_t calc_lighting(const pbr_param_t& pbr_param, const vector3_t& albedo)
	{
		vector3_t F = F_fresenl_schlick_roughness(pbr_param.NoV, pbr_param.f0, pbr_param.roughness);
		vector3_t ks = F;
		vector3_t kd = (vector3_t::one() - ks) * (1.0f - pbr_param.metallic);
		vector3_t irradiance = irradiance_map->sample(pbr_param.n).to_vec3();
		vector3_t diffuse = irradiance * albedo;

		vector3_t r = reflect(pbr_param.n, pbr_param.v);
		r.normalize();

		vector3_t prefiltered_part = sample_trilinear(r, pbr_param.roughness);

		texcoord_t lut_uv(pbr_param.NoV, pbr_param.roughness);

		vector3_t envBRDF = brdf_lut->sample(lut_uv).to_vec3();

		vector3_t specular = prefiltered_part * (F * envBRDF.x + vector3_t::one() * envBRDF.y);

		return kd * diffuse + specular;
	}

private:
	vector3_t sample_trilinear(const vector3_t& r, float roughness)
	{
		assert(prefilter_maps.size() > 1);
		assert(roughness >= 0.0f && roughness <= 1.0f);

		int max_level = (int)(prefilter_maps.size() - 1);
		float lod = roughness * max_level;
		int miplevel0 = (int)lod;
		int miplevel1 = min(miplevel0 + 1, max_level);
		float mip_weight = lod - miplevel0;

		vector3_t prefiltered_part0 = prefilter_maps[miplevel0]->sample(r).to_vec3();
		vector3_t prefiltered_part1 = prefilter_maps[miplevel1]->sample(r).to_vec3();
		return lerp(prefiltered_part0, prefiltered_part1, mip_weight);
	}

};

#endif //__IBL_H__

