#ifndef __SCENE_H__
#define __SCENE_H__

#include <cassert>
#include <vector>
#include <memory>

struct model_t
{
	shading_model_t shading_model;
	cull_mode_t cull_mode;
	shader_resource_t shader_resource;
	mesh_t mesh;
	matrix_t local_to_world;
	model_t(shading_model_t shading_md, cull_mode_t cull_md) :
		shading_model(shading_md),
		cull_mode(cull_md)
	{
		local_to_world.set_identity();
	}
};


class scene_t
{
	std::vector<model_t*> scn_models;

	float view_angle;
	vector3_t eye_pos;
	vector3_t look_at;

	vector3_t light_angle;
	vector4_t light_intensity;

	cube_texture_t sky_env_map;
	mesh_t skybox_mesh;
	matrix_t skybox_world;

	shading_model_t debug_sm = shading_model_t::eSM_MAX;

	bool auto_rotate = true;
	int current_ibl_idx = 0;
	int current_model_idx = 0;

public:
	scene_t() :
		view_angle(0),
		eye_pos(0.4f, 4.5f, 0.25f),
		light_angle(cPI, cPI / 2.0f, 0.0f),
		current_model_idx(0)
	{
	}

	void load()
	{
		view_angle = 0.0f;
		eye_pos.set(0.4f, 4.5f, 0.25f);
		look_at.set(0.0f, 0.0f, 0.0f);

		light_angle.set(cPI, cPI / 2.0f, 0.0f);
		light_intensity.set(1.0f, 1.0f, 1.0f, 8.0f);

		sky_env_map.load_tex("./resource/ibl_textures/env.png", true);
		skybox_mesh.load("./resource/mesh_sphere.obj");
		skybox_world.set_scale(1000.0f, 1000.0f, 1000.0f);
		skybox_world.translate(0, 0, -500.0f);

		model_t* ball = new model_t(shading_model_t::eSM_PBR, cull_mode_t::eCM_CW);
		ball->mesh.load("./resource/mesh_sphere.obj");
		ball->shader_resource.albedo_tex.load_tex("./resource/rustediron2_basecolor.png", true);
		ball->shader_resource.metallic_tex.load_tex("./resource/rustediron2_metallic.png", true);
		ball->shader_resource.roughness_tex.load_tex("./resource/rustediron2_roughness.png", true);
		ball->shader_resource.normal_tex.load_tex("./resource/rustediron2_normal.png", false);
		ball->local_to_world.set_scale(1.6f, 1.6f, 1.6f);
		ball->local_to_world.translate(0, 0, -0.8f);
		scn_models.push_back(ball);

		model_t* gun = new model_t(shading_model_t::eSM_PBR, cull_mode_t::eCM_CW);
		gun->mesh.load("./resource/gun/Air_Gun.obj");
		gun->shader_resource.albedo_tex.load_tex("./resource/gun/Air_Gun_Default_color.jpg", true);
		gun->shader_resource.metallic_tex.load_tex("./resource/gun/Air_Gun_Default_metalness.jpg", true);
		gun->shader_resource.roughness_tex.load_tex("./resource/gun/Air_Gun_Default_roughness.jpg", true);
		gun->shader_resource.normal_tex.load_tex("./resource/gun/Air_Gun_Default_nmap.jpg", false);

		gun->local_to_world.rotate(0.0f, 0.0f, 1.0f, -cPI * 0.25f);
		gun->local_to_world.rotate(0.0f, 1.0f, 0.0f, cPI * 0.25f);
		gun->local_to_world.scale(500.0f, 500.0f, 500.0f);
		scn_models.push_back(gun);

	}

	void render_sky(soft_renderer_t* soft_renderer)
	{
		soft_renderer->set_shading_model(shading_model_t::eSM_Skybox);
		soft_renderer->set_cull_mode(cull_mode_t::eCM_CCW);
		soft_renderer->set_sky_env_map(&sky_env_map);
		soft_renderer->render_model(&skybox_mesh, skybox_world);
	}

	void render(soft_renderer_t* soft_renderer)
	{
		soft_renderer->clear(true, true);

		if (auto_rotate) {
			view_angle += cPI / 128.0f;
		}

		matrix_t mrot;
		mrot.set_rotate(0, 0, 1, view_angle);

		vector3_t eye = mrot.mul_point(eye_pos);
		soft_renderer->set_eye_lookat(eye, look_at);

		soft_renderer->set_light_intensity(light_intensity);

		soft_renderer->set_light_direction(light_angle);

		soft_renderer->draw_cartesian_coordinate();

		// render sky box
		render_sky(soft_renderer);

		auto& m = scn_models[current_model_idx];
		{
			soft_renderer->set_shading_model(m->shading_model);
			soft_renderer->set_cull_mode(m->cull_mode);
			soft_renderer->set_shader_resource(&m->shader_resource);
			soft_renderer->render_model(&m->mesh, m->local_to_world);
		}
	}

	void toggle_auto_rotate() {
		auto_rotate = !auto_rotate;
	}

	void next_ibl(soft_renderer_t* soft_renderer) {
		static const std::string ibl_paths[] = { "./resource/ibl_textures/", "./resource/epic_quad/" };
		++current_ibl_idx;
		if (current_ibl_idx >= array_size(ibl_paths)) {
			current_ibl_idx = 0;
		}
		sky_env_map.load_tex(ibl_paths[current_ibl_idx] + "env.png", true);
		soft_renderer->load_ibl(ibl_paths[current_ibl_idx]);
	}

	void next_model() {
		++current_model_idx;
		if (current_model_idx >= scn_models.size()) {
			current_model_idx = 0;
		}
	}
};

#endif //__SCENE_H__

