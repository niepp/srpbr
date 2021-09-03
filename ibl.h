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
	std::vector<cube_texture_t*> prefilter_maps; // specular environment maps
	texture2d_t *brdf_lut;			// brdf lookup texture
public:

	ibl_t() {
		irradiance_map = new cube_texture_t;
		brdf_lut = new texture2d_t;
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

	void load(const std::string& tex_path)
	{
		std::string irmap_path = tex_path + "ibl_textures/irradiance_map";
		irradiance_map->load_tex(irmap_path, ".tga");

		for (int i = 0; i < 10; ++i) // todo
		{
			std::ostringstream stringStream;
			stringStream << tex_path;
			stringStream << "ibl_textures/prefilter_map_mip";
			stringStream << i;
			std::string copyOfStr = stringStream.str();

			cube_texture_t* pf_map = new cube_texture_t;
			pf_map->load_tex(stringStream.str(), ".tga");
			prefilter_maps.push_back(pf_map);
		}

		/* brdf lookup texture */
		std::string lut_path = tex_path + "brdf_lut.tga";
		brdf_lut->load_tex(lut_path.c_str());

	}

	vector3_t calc_lighting(const pbr_param_t& pbr_param, const vector3_t& albedo)
	{
		vector3_t ks = F_fresenl_schlick_roughness(pbr_param.NoV, pbr_param.f0, pbr_param.roughness);
		vector3_t kd = vector3_t::one() - ks;
		vector3_t irradiance = irradiance_map->sample(pbr_param.n).to_vec3();
		vector3_t diffuse = kd * irradiance * albedo;
		return diffuse;
	}

};

#endif //__IBL_H__

