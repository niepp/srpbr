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
		
		delete brdf_lut;
		brdf_lut = nullptr;
	}

	void load(const std::string& tex_path)
	{
		std::string irmap_path = tex_path + "ibl_textures/irradiance_map";
		irradiance_map->load_tex(irmap_path, ".tga");

		for (int i = 0; i < 10; ++i) // todo
		{
			std::ostringstream stringStream("ibl_textures/prefilter_map_mip");
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

};

#endif //__IBL_H__

