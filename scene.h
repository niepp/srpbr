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


struct scene_t
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

	scene_t() :
		view_angle(0),
		eye_pos(0.4f, 4.5f, 0.25f),
		light_angle(cPI, cPI / 2.0f, 0.0f)
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
		skybox_world.apply_translate(0, 0, -500.0f);

		model_t* gun = new model_t(shading_model_t::eSM_PBR, cull_mode_t::eCM_CW);
		gun->mesh.load("./resource/mesh_sphere.obj");
		gun->shader_resource.albedo_tex.load_tex("./resource/rustediron2_basecolor.png", true);
		gun->shader_resource.metallic_tex.load_tex("./resource/rustediron2_metallic.png", true);
		gun->shader_resource.roughness_tex.load_tex("./resource/rustediron2_roughness.png", true);
		gun->shader_resource.normal_tex.load_tex("./resource/rustediron2_normal.png", false);
		gun->local_to_world.set_scale(1.6f, 1.6f, 1.6f);
		gun->local_to_world.apply_translate(0, 0, -0.8f);

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

		view_angle += cPI / 128.0f;

		matrix_t mrot;
		mrot.set_rotate(0, 0, 1, view_angle);

		vector3_t eye = mrot.mul_point(eye_pos);
		soft_renderer->set_eye_lookat(eye, look_at);

		soft_renderer->set_light_intensity(light_intensity);

		soft_renderer->set_light_direction(light_angle);

		soft_renderer->draw_cartesian_coordinate();

		// render sky box
		render_sky(soft_renderer);

		for (auto &m : scn_models)
		{
			soft_renderer->set_shading_model(m->shading_model);
			soft_renderer->set_cull_mode(m->cull_mode);
			soft_renderer->set_shader_resource(&m->shader_resource);
			soft_renderer->render_model(&m->mesh, m->local_to_world);
		}
	}
};

#endif //__SCENE_H__

