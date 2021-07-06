#include <windows.h>
#include <tchar.h>
#include <iostream>

class Timer
{
public:
	Timer();
	void reset();
	float get_elapse_milliseconds() const;
private:
	uint64_t get_tick() const;
private:
	uint64_t m_start_count;
	static double m_reci_freq;
};

double Timer::m_reci_freq = 0;

Timer::Timer()
{
	LARGE_INTEGER temp;
	::QueryPerformanceFrequency(&temp);
	m_reci_freq = (double)(1000.0 / temp.QuadPart);
	reset();
}

uint64_t Timer::get_tick() const
{
	uint64_t tick = 0;
	::QueryPerformanceCounter((LARGE_INTEGER*)&tick);
	return tick;
}

void Timer::reset()
{
	m_start_count = get_tick();
}

float Timer::get_elapse_milliseconds() const
{
	uint64_t tick = get_tick();
	return (float)((tick - m_start_count) * m_reci_freq);
}


HWND hwnd = 0;
float frame_rate = 30.0f;
Timer timer;

uint32_t display_count = 0;
float fps_refresh = 0.0f;


LRESULT CALLBACK MsgProc(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam)
{
	switch (msg)
	{
	case WM_CREATE:
		break;
	case WM_KEYDOWN:
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

void update()
{

}

void render()
{
	++display_count;
}

void tick()
{
	float interv = 1000.0f / frame_rate;
	float elapse_time = timer.get_elapse_milliseconds();

	if (elapse_time < interv) // 逻辑帧时间未到
	{
		render();
	}
	else
	{
		timer.reset();
		update();
		render();

		const float cRefresh = 500.0f;
		fps_refresh += elapse_time;
		if (fps_refresh > cRefresh)
		{
			float render_rate = 1000.0f * display_count / fps_refresh;
			fps_refresh -= cRefresh;
			display_count = 0;
			TCHAR str[164];
			swprintf(str, 164, _T("fps: %.1f"), render_rate);
			::SetWindowText(hwnd, str);
		}
	}
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
			tick();
		}
	}
}

int main(void)
{
	hwnd = init_window(GetModuleHandle(NULL), _T("sr3d"), 800, 600);

	main_loop();


	return 0;
}
