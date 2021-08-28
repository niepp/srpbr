#ifndef __IBL_H__
#define __IBL_H__

class texture2d_t;
class cube_texture_t;

class ibl_t
{
	cube_texture_t* irradiance_map;
	cube_texture_t* prefilter_maps;
	texture2d_t* brdf_lut;
public:
	void load(const char* tex_path)
	{
	}


};

#endif //__IBL_H__

