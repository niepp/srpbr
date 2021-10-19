#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <tchar.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <cassert>
#include <chrono>

#include "math3d.h"
#include "pbr_common.h"
#include "texture.h"
#include "ibl.h"
#include "model.h"

#include "pre_compute.h"

#include "soft_renderer.h"


float reci_freq;
int64_t tick_start;
uint32_t frame_count = 0, frame_rate = 0;
HWND hwnd = 0;
HDC screenDC;

int width = 800;
int height = 600;
model_t sphere_model;
float view_angle = 0;
vector3_t eye_pos(0.4f, 4.5f, 0.25f);
vector3_t light_angle(cPI, cPI / 2.0f, 0.0f);
soft_renderer_t *soft_renderer = nullptr;

template <typename T, int n>
int array_size(T(&)[n])
{
	return n;
}

inline long long get_now_ms()
{
	auto tp_now = std::chrono::high_resolution_clock::now();
	return std::chrono::duration_cast<std::chrono::milliseconds>(tp_now.time_since_epoch()).count();
}


void render_scene()
{
	soft_renderer->clear(true, true);

	view_angle += cPI / 128.0f;

	matrix_t mrot;
	mrot.set_rotate(0, 0, 1, view_angle);

	vector3_t eye = mrot.mul_point(eye_pos);
	soft_renderer->update_eye_position(eye);

	soft_renderer->draw_cartesian_coordinate();

	matrix_t world;
	world.set_scale(1000.0f, 1000.0f, 1000.0f);
	world.apply_translate(0, 0, -500.0f);

	shading_model_t old_shading_model = soft_renderer->get_shading_model();
	if (old_shading_model != shading_model_t::eSM_Wireframe) {
		soft_renderer->set_shading_model(shading_model_t::eSM_Skybox);
	}
	soft_renderer->set_cull_mode(cull_mode_t::eCM_CCW);
	soft_renderer->render_model(&sphere_model, world);
	soft_renderer->set_shading_model(old_shading_model);

	world.set_scale(1.6f, 1.6f, 1.6f);
	world.apply_translate(0, 0, -0.8f);

	soft_renderer->set_cull_mode(cull_mode_t::eCM_CW);
	soft_renderer->render_model(&sphere_model, world);

}


void main_loop()
{
	// Show the window
	::ShowWindow(hwnd, SW_SHOWDEFAULT);
	::UpdateWindow(hwnd);

	// Enter the message loop
	MSG msg;
	ZeroMemory(&msg, sizeof(msg));
	while (msg.message != WM_QUIT)
	{
		if (::PeekMessage(&msg, NULL, 0U, 0U, PM_REMOVE))
		{
			::TranslateMessage(&msg);
			::DispatchMessage(&msg);
		}
		else
		{
			int64_t tick_now = 0;
			::QueryPerformanceCounter((LARGE_INTEGER*)&tick_now);
			float dt_ms = (tick_now - tick_start) * reci_freq;
			if (++frame_count >= 90) // limit to fps 90
			{
				if (dt_ms < 1000.0f) {
					::Sleep(1000 - (DWORD)dt_ms);
				}
			}
			if (dt_ms >= 1000.0f)
			{
				frame_rate = (uint32_t)(frame_count * 1000.0f / dt_ms);
				tick_start = tick_now;
				frame_count = 0;
				TCHAR str[64];
				::wsprintfW(str, _T("srpbr %d fps"), frame_rate);
				::SetWindowText(hwnd, str);
			}

			render_scene();

			HDC hDC = GetDC(hwnd);
			BitBlt(hDC, 0, 0, width, height, screenDC, 0, 0, SRCCOPY);
			ReleaseDC(hwnd, hDC);

		}
	}
}

LRESULT CALLBACK MsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch (msg)
	{
	case WM_CREATE:
		break;
	case WM_KEYDOWN:
		switch (wParam & 0x1FF)
		{
		case 'C':
		{
			shading_model_t shading_model = (shading_model_t)((int)soft_renderer->get_shading_model() + 1);
			if (shading_model == shading_model_t::eSM_Skybox) {
				shading_model = (shading_model_t)((int)shading_model + 1);
			}
			if (shading_model == shading_model_t::eSM_MAX) {
				shading_model = shading_model_t::eSM_Color;
			}
			soft_renderer->set_shading_model(shading_model);
			std::cout << (int)shading_model << std::endl;
		}
			break;
		case VK_UP:
			light_angle.y += 0.1f;
			soft_renderer->update_light(light_angle);
			break;
		case VK_DOWN:
			light_angle.y -= 0.1f;
			soft_renderer->update_light(light_angle);
			break;
		case VK_LEFT:
			light_angle.x -= 0.1f;
			soft_renderer->update_light(light_angle);
			break;
		case VK_RIGHT:
			light_angle.x += 0.1f;
			soft_renderer->update_light(light_angle);
			break;
		case 'P':
			soft_renderer->save_framebuffer("./result/framebuffer_" + std::to_string(get_now_ms()) + ".png");
			break;
		default:
			break;
		}
		break;
	case WM_SIZE:
	case WM_EXITSIZEMOVE:
		break;
	case WM_DESTROY:
		PostQuitMessage(0);
		return 0;
	}
	return DefWindowProc(hWnd, msg, wParam, lParam);
}

HWND init_window(HINSTANCE instance, const TCHAR* title, int width, int height)
{
	const TCHAR class_name[] = _T("wndclass");

	// Register the window class
	WNDCLASSEX wc = {
		sizeof(WNDCLASSEX), CS_CLASSDC, MsgProc, 0, 0,
		instance, NULL, NULL, NULL, NULL,
		class_name, NULL
	};

	RegisterClassEx(&wc);

	// calculate client size，设置窗口客户区大小
	RECT clientSize;
	clientSize.top = 0;
	clientSize.left = 0;
	clientSize.right = width;
	clientSize.bottom = height;
	DWORD style = WS_OVERLAPPEDWINDOW | WS_BORDER | WS_CAPTION | WS_CLIPCHILDREN | WS_CLIPSIBLINGS;
	//计算客户区矩形大小
	::AdjustWindowRect(&clientSize, style, FALSE);
	int w = clientSize.right - clientSize.left;
	int h = clientSize.bottom - clientSize.top;

	// Create the application's window
	HWND hwnd = CreateWindow(class_name, title,
		style, 0, 0, w, h,
		NULL, NULL, instance, NULL);

	::ShowWindow(hwnd, SW_SHOWDEFAULT);
	::UpdateWindow(hwnd);

	return hwnd;

}

int main(void)
{
	//generate_irradiance_map("./resource/ibl_textures/env.png", "./resource/ibl_textures/irradiance.png");
	////generate_prefilter_envmap("./resource/ibl_textures/env.png", "./resource/ibl_textures/prefilter");
	////generate_BRDF_LUT("./resource/brdf_lut.png");
	//return 0;

	hwnd = init_window(GetModuleHandle(NULL), _T(""), width, height);

	screenDC = CreateCompatibleDC(GetDC(hwnd));

	BITMAPINFO bi = {
		{ sizeof(BITMAPINFOHEADER), width, height, 1u, 32u, BI_RGB, width * height * 4u, 0, 0, 0, 0 }
	};

	uint32_t* framebuffer = nullptr;
	HBITMAP screenBMP = CreateDIBSection(screenDC, &bi, DIB_RGB_COLORS, (void**)&framebuffer, 0, 0);
	if (screenBMP == NULL) {
		return -1;
	}

	SelectObject(screenDC, screenBMP);

	LARGE_INTEGER temp;
	::QueryPerformanceFrequency(&temp);
	reci_freq = 1000.0f / temp.QuadPart;
	::QueryPerformanceCounter((LARGE_INTEGER*)&tick_start);

	sphere_model.load("./resource/mesh_sphere.obj");

	soft_renderer = new soft_renderer_t(width, height, framebuffer, eye_pos);

	main_loop();

	delete soft_renderer;
	return 0;

}

