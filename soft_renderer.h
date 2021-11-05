#ifndef __SOFT_RENDERER_H__
#define __SOFT_RENDERER_H__
#include <memory>

/*
 有一边是水平的三角形
 l----------r
  \        /
   \      /
	\    /
	 \  /
	  \/
	   p
*/
struct scan_tri_t
{
	interp_vertex_t p;
	interp_vertex_t l, r;
};

enum class shading_model_t
{
	eSM_Wireframe = 0,
	eSM_Color,
	eSM_Phong,
	eSM_PBR,
	eSM_Skybox,
	eSM_MAX,
};

enum class cull_mode_t
{
	eCM_None,
	eCM_CW,
	eCM_CCW,
};

struct uniformbuffer_t
{
	vector3_t eye;
	matrix_t view;
	matrix_t proj;
	matrix_t viewproj;
	vector3_t light_dir;
	vector3_t light_intensity;
	float specular_power;
};

struct shader_resource_t
{
	texture2d_t albedo_tex;
	texture2d_t metallic_tex;
	texture2d_t roughness_tex;
	texture2d_t normal_tex;
};

struct camera_t
{
	vector3_t eye;
	vector3_t lookat;
	float fovy;
	float aspect;
	float pnear;
	float pfar;
};

class soft_renderer_t
{
private:
	int width;
	int height;
	uniformbuffer_t uniformbuffer;
	uint32_t* framebuffer;
	float* zbuffer;
	pbr_param_t pbr_param;
	ibl_t ibl;
	shading_model_t shading_model;
	cull_mode_t cull_mode;

	cube_texture_t *sky_env_map; // environment map
	shader_resource_t* shader_resource;

	camera_t camera;

public:
	soft_renderer_t(int w, int h, uint32_t* fb, float infovy = cPI * 0.5f, float innear = 0.1f, float infar = 2000.0f) :
		width(w), height(h),
		framebuffer(fb),
		zbuffer(nullptr),
		shading_model(shading_model_t::eSM_PBR),
		cull_mode(cull_mode_t::eCM_CW)
	{
		zbuffer = new float[width * height];
		for (int i = 0; i < width * height; ++i) {
			zbuffer[i] = 0;
		}

		ibl.load("./resource/ibl_textures/");

		camera.fovy = infovy;
		camera.aspect = 1.0f * width / height;
		camera.pnear = innear;
		camera.pfar = infar;

		// setup perspective matrix
		uniformbuffer.proj.set_perspective(camera.fovy, camera.aspect, camera.pnear, camera.pfar);

	}

	~soft_renderer_t()
	{
		delete[] zbuffer;
		framebuffer = nullptr;
		zbuffer = nullptr;
	}

	void on_change_size(uint32_t* newfb, int neww, int newh);

	shading_model_t get_shading_model() const { return shading_model; }
	void set_shading_model(shading_model_t in_shading_model) { shading_model = in_shading_model; }

	cull_mode_t get_cull_mode() const { return cull_mode; }
	void set_cull_mode(cull_mode_t in_cull_mode) { cull_mode = in_cull_mode; }
	void set_sky_env_map(cube_texture_t* skymap) { sky_env_map = skymap; }
	void set_shader_resource(shader_resource_t* srv) { shader_resource = srv; }

	void clear(bool bclear_fb, bool bclear_zb);
	void save_framebuffer(const std::string& fb_texpath);
	void save_depthbuffer(const std::string& depth_texpath);
	void set_light_direction(const vector3_t& light_angle);
	void set_light_intensity(const vector4_t& light_inten);
	void set_eye_lookat(const vector3_t& eye, const vector3_t& at);
	void render_model(mesh_t* mesh, const matrix_t& world);
	void draw_cartesian_coordinate();

private:
	void write_pixel(int x, int y, uint32_t color);
	void write_depth(int x, int y, float z);
	bool depth_test(int x, int y, float z);
	vector3_t normal_mapping(const vector3_t& vt_n, const vector3_t& ts_n);
	void phong_shading(const interp_vertex_t& p, vector4_t& out_color);
	void pbr_shading(const interp_vertex_t& p, vector4_t& out_color);
	void pixel_process(int x, int y, const interp_vertex_t& p);

	void scan_horizontal(const interp_vertex_t& vl, const interp_vertex_t& vr, int y);
	void scan_triangle(scan_tri_t* sctri);
	bool check_clip(const vector4_t* p, int width, int height);
	bool check_clip_triangle(const vector4_t* p0, const vector4_t* p1, const vector4_t* p2);
	void perspective_divide(interp_vertex_t* p);
	void to_screen_coord(vector4_t* p);
	void vertex_process(const matrix_t& world, const matrix_t& mvp, const mesh_vertex_t& v, interp_vertex_t& p);
	void draw_triangle(const interp_vertex_t& p0, const interp_vertex_t& p1, const interp_vertex_t& p2);
	void draw_line(const vector4_t& p0, const vector4_t& p1, uint32_t c);
};

void soft_renderer_t::on_change_size(uint32_t* newfb, int neww, int newh)
{
	framebuffer = newfb;
	width = neww;
	height = newh;

	delete[] zbuffer;

	zbuffer = new float[width * height];
	for (int i = 0; i < width * height; ++i) {
		zbuffer[i] = 0;
	}

	camera.aspect = 1.0f * width / height;

	// setup perspective matrix
	uniformbuffer.proj.set_perspective(camera.fovy, camera.aspect, camera.pnear, camera.pfar);

	// cache view X proj
	uniformbuffer.viewproj = mul(uniformbuffer.view, uniformbuffer.proj);
}

void soft_renderer_t::write_pixel(int x, int y, uint32_t color)
{
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);
	framebuffer[(height - y - 1) * width + x] = color;
}

void soft_renderer_t::write_depth(int x, int y, float z)
{
	// zbuffer store the NDC(Normalized Device Coordinates) reverse-z value
	// https://developer.nvidia.com/content/depth-precision-visualized
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);
	zbuffer[(height - y - 1) * width + x] = z;
}

bool soft_renderer_t::depth_test(int x, int y, float z)
{
	assert(x >= 0 && x < width);
	assert(y >= 0 && y < height);
	float nowz = zbuffer[(height - y - 1) * width + x];
	return z > nowz;
}

vector3_t soft_renderer_t::normal_mapping(const vector3_t& vt_n, const vector3_t& ts_n)
{
	// [ tan_x ]
	// | tan_y |
	// [ n     ]
	// makeup a rotation matrix
	vector3_t up = fabs(vt_n.z) < 0.999f ? vector3_t(0.0f, 0.0f, 1.0f) : vector3_t(1.0f, 0.0f, 0.0f);
	vector3_t tan_x = cross(up, vt_n);
	tan_x.normalize();
	vector3_t tan_y = cross(vt_n, tan_x);
	tan_y.normalize();

	// convert h to world space
	vector3_t ws_n = tan_x * ts_n.x + tan_y * ts_n.y + vt_n * ts_n.z;
	ws_n.normalize();
	return ws_n;
}

void soft_renderer_t::phong_shading(const interp_vertex_t& p, vector4_t& out_color)
{
	texcoord_t uv = p.uv;
	uv.u /= p.pos.w;
	uv.v /= p.pos.w;

	vector3_t vt_n = p.nor / p.pos.w;
	vt_n.normalize();

	vector3_t wpos = p.wpos / p.pos.w;
	vector3_t l = -uniformbuffer.light_dir;

	vector4_t albedo = shader_resource->albedo_tex.sample(uv);
	vector3_t tangent_space_n = shader_resource->normal_tex.sample(uv).to_vec3();

	tangent_space_n = tangent_space_n * 2.0f - vector3_t::one();
	tangent_space_n.normalize();
	vector3_t n = normal_mapping(vt_n, tangent_space_n);
	n.normalize();

	float NoL = max(dot(n, l), 0.0f);

	vector3_t diffuse = albedo.to_vec3() * uniformbuffer.light_intensity * NoL;

	vector3_t v = uniformbuffer.eye - wpos;
	v.normalize();

	// (v + l) / 2
	vector3_t h = v + l;
	h.normalize();

	float NoH = max(dot(n, h), 0.0f);
	vector3_t ks(0.5f, 0.5f, 0.5f);
	vector3_t specular = ks * uniformbuffer.light_intensity * max(0.0f, pow(NoH, uniformbuffer.specular_power));

	vector3_t radiance = diffuse + specular;

	vector3_t ldr = reinhard_mapping(radiance);

	out_color = gamma_correction(ldr);

}

void soft_renderer_t::pbr_shading(const interp_vertex_t& p, vector4_t& out_color)
{
	texcoord_t uv = p.uv;
	uv.u /= p.pos.w;
	uv.v /= p.pos.w;

	vector3_t vt_n = p.nor / p.pos.w;
	vt_n.normalize();

	vector3_t wpos = p.wpos / p.pos.w;

	vector3_t tangent_space_n = shader_resource->normal_tex.sample(uv).to_vec3();
	vector3_t albedo = shader_resource->albedo_tex.sample(uv).to_vec3();
	vector4_t metallic_texel = shader_resource->metallic_tex.sample(uv);
	vector4_t roughness_texel = shader_resource->roughness_tex.sample(uv);

	tangent_space_n = tangent_space_n * 2.0f - vector3_t::one();
	tangent_space_n.normalize();
//	pbr_param.n = normal_mapping(vt_n, tangent_space_n);
	pbr_param.n = vt_n;

	pbr_param.n.normalize();

	pbr_param.v = uniformbuffer.eye - wpos;
	pbr_param.v.normalize();

	pbr_param.l = -uniformbuffer.light_dir;

	// (v + l) / 2
	pbr_param.h = pbr_param.v + pbr_param.l;
	pbr_param.h.normalize();

	pbr_param.NoL = max(dot(pbr_param.n, pbr_param.l), 0);
	pbr_param.NoH = max(dot(pbr_param.n, pbr_param.h), 0);
	pbr_param.NoV = max(dot(pbr_param.n, pbr_param.v), 0);
	pbr_param.HoV = max(dot(pbr_param.h, pbr_param.v), 0);

	pbr_param.metallic = metallic_texel.r;
	pbr_param.roughness = roughness_texel.r;

	pbr_param.roughness = clamp(pbr_param.roughness, 0.02f, 1.0f); // roughness为0，对应着完全的反射，对于非面积光源，会产生无穷大的反射值（能量全部集中到了面积为0的一点上）

	float a = pbr_param.roughness * pbr_param.roughness;

	// F：菲涅尔系数，是光线发生反射与折射的比例，当光线垂直进入表面时的菲涅尔系数记为F0，F0是可以实际测量出来的
	// 不同的材质，这个F0是不同的，F0越大，感觉是越明亮，金属材质通常F0比较大，而且是有颜色的（不同频率的光，对应菲涅尔系数不同），非金属材质最低也有一点点反射，就取(0.04, 0.04, 0.04)，
	// 渲染中通常用一个权重系数（金属度）在非金属材质的最小菲涅尔系数和具有最强金属性的材质菲涅尔系数（自身漫反射颜色）之间进行线性插值得到F0.
	vector3_t temp(0.04f, 0.04f, 0.04f);
	pbr_param.f0 = lerp(temp, albedo, pbr_param.metallic);

	vector3_t F = F_fresenl_schlick(pbr_param.HoV, pbr_param.f0);

	float D = D_Trowbridge_Reitz_GGX(a, pbr_param.NoH);

	float V = V_Schlick_GGX(pbr_param.roughness, pbr_param.NoV, pbr_param.NoL);

	if (!is_valid(V)) {
		assert(0);
	}

	vector3_t cook_torrance_brdf = F * D * V;

	//kd和ks分别是漫反射和镜面反射系数，表征了入射能量在镜面反射和漫反射之间进行分配，以满足能量守恒定律。 因此kd + ks < 1
	//F代表了材质表面反射率，那么我们可以直接让ks = F，然后令kd = (1 - ks) * (1 - metallic)。这实际上是将入射能量在镜面反射和漫反射之间进行了分配，以满足能量守恒定律
	vector3_t ks = F;
	vector3_t kd = (vector3_t::one() - F) * (1.0f - pbr_param.metallic);

	vector3_t directIrradiance = uniformbuffer.light_intensity * pbr_param.NoL;

	vector3_t diffuse_brdf = albedo / cPI;

	vector3_t brdf = kd * diffuse_brdf + cook_torrance_brdf;

	vector3_t radiance = brdf * directIrradiance;

	vector3_t indirect_radiance = ibl.calc_lighting(pbr_param, albedo);

	radiance += indirect_radiance;

	vector3_t ldr = reinhard_mapping(radiance);

	out_color = gamma_correction(ldr);

}

void soft_renderer_t::pixel_process(int x, int y, const interp_vertex_t& p)
{
	vector4_t color;
	switch (shading_model)
	{
	case shading_model_t::eSM_Color:
		color = p.color / p.pos.w;
		break;
	case shading_model_t::eSM_Phong:
		phong_shading(p, color);
		break;
	case shading_model_t::eSM_PBR:
		pbr_shading(p, color);
		break;
	case shading_model_t::eSM_Skybox:
	{
		vector3_t n = p.nor / p.pos.w;
		n.normalize();
		vector3_t radiance = sky_env_map->sample(n).to_vec3();
		color = gamma_correction(radiance);
	}
	break;
	default:
		break;
	}
	std::swap(color.r, color.b); // windows HDC is bgra format!
	write_pixel(x, y, makefour(color));
	write_depth(x, y, p.pos.z);
}

void soft_renderer_t::scan_horizontal(const interp_vertex_t& vl, const interp_vertex_t& vr, int y)
{
	float dist = vr.pos.x - vl.pos.x;
	// left要往小取整，right要往大取整，避免三角形之间的接缝空隙！
	int left = (int)(vl.pos.x);
	int right = (int)(vr.pos.x + 0.5f);
	left = clamp(left, 0, width - 1);
	right = clamp(right, 0, width - 1);
	for (int i = left; i < right; ++i)
	{
		float w = (i - vl.pos.x) / dist;
		w = clamp(w, 0.0f, 1.0f);
		interp_vertex_t p = lerp(vl, vr, w);
		if (depth_test(i, y, p.pos.z)) {
			pixel_process(i, y, p);
		}
	}
}

void soft_renderer_t::scan_triangle(scan_tri_t* sctri)
{
	if (sctri->l.pos.x > sctri->r.pos.x) {
		std::swap(sctri->l, sctri->r);
	}

	assert(sctri->l.pos.y == sctri->r.pos.y);
	assert(sctri->l.pos.y != sctri->p.pos.y);

	float ymax = std::fmax(sctri->p.pos.y, sctri->l.pos.y);
	float ymin = std::fmin(sctri->p.pos.y, sctri->l.pos.y);
	float ydist = sctri->p.pos.y - sctri->l.pos.y;

	// bottom要往小取整，top要往大取整，避免三角形之间的接缝空隙！
	int bottom = (int)(ymin);
	int top = (int)(ymax + 0.5f);
	bottom = clamp(bottom, 0, height - 1);
	top = clamp(top, 0, height - 1);

	for (int i = bottom; i < top; ++i)
	{
		float cury = i + 0.0f;
		float w = (ydist > 0 ? cury - ymin : cury - ymax) / ydist;
		w = clamp(w, 0.0f, 1.0f);
		interp_vertex_t vl = lerp(sctri->l, sctri->p, w);
		interp_vertex_t vr = lerp(sctri->r, sctri->p, w);
		scan_horizontal(vl, vr, i);
	}
}

bool soft_renderer_t::check_clip(const vector4_t* p, int width, int height)
{
	if (p->x < 0 || p->x >= width) return false;
	if (p->y < 0 || p->y >= height) return false;
	if (p->z < 0.0f || p->z > 1.0f) return false;
	return true;
}

bool soft_renderer_t::check_clip_triangle(const vector4_t* p0, const vector4_t* p1, const vector4_t* p2)
{
	return check_clip(p0, width, height)
		|| check_clip(p1, width, height)
		|| check_clip(p2, width, height);
}

void soft_renderer_t::perspective_divide(interp_vertex_t* p)
{
	// the point with homogeneous coordinates [x, y, z, w] corresponds to the three-dimensional Cartesian point [x/w, y/w, z/w].
	// to cvv coord x[-1, 1], y[-1, 1], z[0, 1]
	float revw = 1.0f / p->pos.w;
	p->pos.x *= revw;
	p->pos.y *= revw;
	p->pos.z *= revw;
	p->pos.w = revw; // 这里存1/w，因为透视校正原因，1/w才有线性关系： 1/p.w = lerp(1/p0.w, 1/p1.w, weight)

	p->wpos.x *= revw;
	p->wpos.y *= revw;
	p->wpos.z *= revw;

	p->nor.x *= revw;
	p->nor.y *= revw;
	p->nor.z *= revw;

	p->uv.u *= revw;
	p->uv.v *= revw;

	p->color.r *= revw;
	p->color.g *= revw;
	p->color.b *= revw;
	p->color.a *= revw;

}

void soft_renderer_t::to_screen_coord(vector4_t* p)
{
	// to screen coord x[0, width], y[0, height], z[0, 1] (depth)
	p->x = (p->x + 1.0f) * 0.5f * width;
	p->y = height - 1 - (p->y + 1.0f) * 0.5f * height;
}

void soft_renderer_t::vertex_process(const matrix_t& world, const matrix_t& mvp, const mesh_vertex_t& v, interp_vertex_t& p)
{
	p.pos = mvp * vector4_t(v.pos, 1.0f);
	p.wpos = world.mul_point(v.pos);
	p.nor = world.mul_vector(v.nor); // suppose world contain NO no-uniform scale!
	p.color = v.color;
	p.uv = v.uv;
	perspective_divide(&p);
	to_screen_coord(&p.pos);
}

void soft_renderer_t::draw_triangle(const interp_vertex_t& p0, const interp_vertex_t& p1, const interp_vertex_t& p2)
{
	// degenerate triangle
	if (p0.pos.y == p1.pos.y && p1.pos.y == p2.pos.y) return;
	if (p0.pos.x == p1.pos.x && p1.pos.x == p2.pos.x) return;

	scan_tri_t uptri, downtri;
	if (p0.pos.y == p1.pos.y) // up triangle
	{
		uptri.p = p2;
		uptri.l = p0;
		uptri.r = p1;
		scan_triangle(&uptri);
	}
	else if (p1.pos.y == p2.pos.y) // down triangle
	{
		downtri.p = p0;
		downtri.l = p1;
		downtri.r = p2;
		scan_triangle(&downtri);
	}
	else
	{
		float w = (p1.pos.y - p0.pos.y) / (p2.pos.y - p0.pos.y);
		interp_vertex_t mid = lerp(p0, p2, w);
		mid.pos.y = p1.pos.y;

		uptri.p = p2;
		uptri.l = mid;
		uptri.r = p1;

		downtri.p = p0;
		downtri.l = mid;
		downtri.r = p1;

		scan_triangle(&uptri);
		scan_triangle(&downtri);

	}

}

void soft_renderer_t::draw_line(const vector4_t& p0, const vector4_t& p1, uint32_t c)
{
	vector4_t v0 = p0;
	vector4_t v1 = p1;

	if (!(check_clip(&v0, width, height) || check_clip(&v0, width, height))) {
		return;
	}

	float dx = std::abs(v0.x - v1.x);
	float dy = std::abs(v0.y - v1.y);

	if (dy <= dx)
	{
		if (v0.x > v1.x) {
			std::swap(v0, v1);
		}

		int xmin = clamp((int)(v0.x), 0, width - 1);
		int xmax = clamp((int)(v1.x + 0.5f), 0, width - 1);
		float y = v0.y;
		int iy = (int)y;
		if (iy >= 0 && iy < height) {
			write_pixel(xmin, iy, c);
		}

		if (xmax > xmin) {
			float delta = (v1.y - v0.y) / dx;
			for (int x = xmin; x < xmax; ++x)
			{
				y += delta;
				iy = (int)y;
				if (iy < 0 || iy >= height) {
					break;
				}
				write_pixel(x, iy, c);
			}
		}
	}
	else {
		if (v0.y > v1.y) {
			std::swap(v0, v1);
		}

		int ymin = clamp((int)(v0.y), 0, height - 1);
		int ymax = clamp((int)(v1.y + 0.5f), 0, height - 1);
		float x = v0.x;
		int ix = (int)x;
		if (ix >= 0 && ix < width) {
			write_pixel(ix, ymin, c);
		}

		if (ymax > ymin) {
			float delta = (v1.x - v0.x) / dy;
			for (int y = ymin; y < ymax; ++y)
			{
				x += delta;
				ix = (int)x;
				if (ix < 0 || ix >= width) {
					break;
				}
				write_pixel(ix, y, c);
			}
		}
	}

}

void soft_renderer_t::draw_cartesian_coordinate()
{
	static const vector3_t o(0, 0, 0);
	static const vector3_t x(1, 0, 0);
	static const vector3_t y(0, 1, 0);
	static const vector3_t z(0, 0, 1);

	auto transform_to_screen = [this](const matrix_t& vp, const vector3_t& p) -> vector4_t
	{
		vector4_t tp = vp * vector4_t(p, 1.0f);
		float revw = std::abs(1.0f / tp.w);
		tp.x *= revw;
		tp.y *= revw;
		tp.z *= revw;
		tp.w = revw;
		to_screen_coord(&tp);
		tp.x = clamp(tp.x, 0.0f, (width - 1) * 1.0f);
		tp.y = clamp(tp.y, 0.0f, (height - 1) * 1.0f);
		return tp;
	};

	vector4_t to = transform_to_screen(uniformbuffer.viewproj, o);
	vector4_t tx = transform_to_screen(uniformbuffer.viewproj, x);
	vector4_t ty = transform_to_screen(uniformbuffer.viewproj, y);
	vector4_t tz = transform_to_screen(uniformbuffer.viewproj, z);

	draw_line(to, tx, 0xffff0000);
	draw_line(to, ty, 0xff00ff00);
	draw_line(to, tz, 0xff0000ff);
}

void soft_renderer_t::clear(bool bclear_fb, bool bclear_zb)
{
	if (bclear_fb) {
		memset(framebuffer, 0, width * height * sizeof(uint32_t));
	}

	if (bclear_zb) {
		for (int i = 0; i < width * height; ++i) {
			zbuffer[i] = 0;
		}
	}
}

void soft_renderer_t::save_framebuffer(const std::string& fb_texpath)
{
	unsigned int* data = new unsigned int[width * height];
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			int src_idx = i * width + j;
			int dst_idx = (height - 1 - i) * width + j;
			vector4_t color = to_color(framebuffer[src_idx]);
			std::swap(color.r, color.b);
			data[dst_idx] = makefour(color);
		}
	}
	stbi_write_png(fb_texpath.c_str(), width, height, 4, data, width * 4);
	delete[] data;
}

void soft_renderer_t::save_depthbuffer(const std::string& depth_texpath)
{
	unsigned int* data = new unsigned int[width * height];
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			int src_idx = i * width + j;
			int dst_idx = (height - 1 - i) * width + j;
			float z = zbuffer[src_idx];
			data[dst_idx] = makefour(vector4_t(z, z, z, 1.0f));
		}
	}
	stbi_write_png(depth_texpath.c_str(), width, height, 4, data, width * 4);
	delete[] data;
}

void soft_renderer_t::set_light_direction(const vector3_t& light_angle)
{
	float angle_y = clamp(light_angle.y, 0.0f, cPI);
	float cosw = cos(angle_y);
	float sinw = sin(angle_y);
	uniformbuffer.light_dir.x = sin(light_angle.x) * sinw;
	uniformbuffer.light_dir.y = cos(light_angle.x) * sinw;
	uniformbuffer.light_dir.z = cosw;
}

void soft_renderer_t::set_light_intensity(const vector4_t& light_inten)
{
	uniformbuffer.light_intensity = light_inten.to_vec3();
	uniformbuffer.specular_power = light_inten.w;
}

void soft_renderer_t::set_eye_lookat(const vector3_t &eye, const vector3_t &at)
{
	camera.eye = eye;
	camera.lookat = at;

	uniformbuffer.eye = eye;

	// setup view matrix
	const vector3_t up(0.0f, 0.0f, 1.0f); // z axis
	uniformbuffer.view.set_lookat(uniformbuffer.eye, at, up);

	// cache view X proj
	uniformbuffer.viewproj = mul(uniformbuffer.view, uniformbuffer.proj);
}

void soft_renderer_t::render_model(mesh_t* mesh, const matrix_t &world)
{
	mesh_vertex_vec_t& vb = mesh->m_mesh_vertex;
	interp_vertex_vec_t& vb_post = mesh->m_vertex_post;
	index_vec_t& ib = mesh->m_mesh_indices;

	matrix_t mvp = mul(world, uniformbuffer.viewproj);

	for (int i = 0; i < vb.size(); ++i)
	{
		vertex_process(world, mvp, vb[i], vb_post[i]);
	}

	int tri_num = (int)ib.size() / 3;
	for (int i = 0; i < tri_num; ++i)
	{
		int i0 = ib[i * 3 + 0];
		int i1 = ib[i * 3 + 1];
		int i2 = ib[i * 3 + 2];

		vector4_t v01 = vb_post[i1].pos - vb_post[i0].pos;
		vector4_t v02 = vb_post[i2].pos - vb_post[i0].pos;

		float det_xy = v01.x * v02.y - v01.y * v02.x;
		if (cull_mode == cull_mode_t::eCM_CW && det_xy >= 0.0f) {
			continue; // backface culling
		}
		else if (cull_mode == cull_mode_t::eCM_CCW && det_xy < 0.0f) {
			continue; // forward face culling
		}

		const interp_vertex_t* p0 = &vb_post[i0];
		const interp_vertex_t* p1 = &vb_post[i1];
		const interp_vertex_t* p2 = &vb_post[i2];

		if (!check_clip_triangle(&p0->pos, &p1->pos, &p2->pos)) {
			continue;
		}

		// make sure p0y <= p1y <= p2y
		if (p0->pos.y > p1->pos.y) std::swap(p0, p1);
		if (p0->pos.y > p2->pos.y) std::swap(p0, p2);
		if (p1->pos.y > p2->pos.y) std::swap(p1, p2);

		if (shading_model == shading_model_t::eSM_Wireframe)
		{
			draw_line(p0->pos, p1->pos, 0xffffffff);
			draw_line(p1->pos, p2->pos, 0xffffffff);
			draw_line(p2->pos, p0->pos, 0xffffffff);
		}
		else
		{
			draw_triangle(*p0, *p1, *p2);
		}
	}

}

#endif //__SOFT_RENDERER_H__

