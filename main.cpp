#define _CRT_SECURE_NO_WARNINGS
#include <windows.h>
#include <tchar.h>
#include <math.h>
#include <vector>
#include <atomic>
#include <thread>
#include <iostream>
#include <cassert>
#include <chrono>

#include "math3d.h"
#include "pbr_common.h"
#include "texture.h"
#include "ibl.h"
#include "mesh.h"

#include "pre_compute.h"
#include "soft_renderer.h"
#include "scene.h"

enum class command_t {
	cNone = 0,
	cSaveFramebuffer,
	cSaveDepthbuffer,
};

float reci_freq;
int64_t tick_start;
uint32_t frame_count = 0, frame_rate = 0;
HWND hwnd = 0;
HDC screenDC;

int width = 800;
int height = 600;

scene_t scn;
soft_renderer_t* soft_renderer = nullptr;

std::atomic<command_t> cmd(command_t::cNone);

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

void read_console(std::atomic<command_t>& readcmd)
{
	std::cout << " --------------------------------- " << std::endl;	
	std::cout << " c.save_fb : save framebuffer to result " << std::endl;
	std::cout << " c.save_depth : save depthbuffer to result " << std::endl;
	std::cout << " --------------------------------- \n" << std::endl;
	std::string buffer;
	while (true) {
		std::cin >> buffer;
		if (buffer == "c.save_fb") {
			readcmd.store(command_t::cSaveFramebuffer);
		}
		else if (buffer == "c.save_depth") {
			readcmd.store(command_t::cSaveDepthbuffer);
		}
	}
}

void on_console_cmd()
{
	switch (cmd.load())
	{
	case command_t::cSaveFramebuffer:
		soft_renderer->save_framebuffer("./result/framebuffer_" + std::to_string(get_now_ms()) + ".png");
		break;
	case command_t::cSaveDepthbuffer:
		soft_renderer->save_depthbuffer("./result/depthbuffer_" + std::to_string(get_now_ms()) + ".png");		
		break;
	}
	cmd.store(command_t::cNone);
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

			on_console_cmd();
		
			scn.render(soft_renderer);

			HDC hDC = GetDC(hwnd);
			BitBlt(hDC, 0, 0, width, height, screenDC, 0, 0, SRCCOPY);
			ReleaseDC(hwnd, hDC);

		}
	}
}

HWND init_window(HINSTANCE instance, const TCHAR* title, int width, int height)
{
	const TCHAR class_name[] = _T("wndclass");

	// Register the window class
	WNDCLASSEX wc = {
		sizeof(WNDCLASSEX), CS_CLASSDC, DefWindowProc, 0, 0,
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

	std::thread readcmd_thread(read_console, std::ref(cmd));

	scn.load();
	soft_renderer = new soft_renderer_t(width, height, framebuffer);

	main_loop();

	delete soft_renderer;

	return 0;

}

