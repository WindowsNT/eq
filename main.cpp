#include "stdafx.h"

HWND MainWindow = 0;
HINSTANCE hAppInstance = 0;
HICON hIcon1 = 0;

#include "f:\\TOOLS\\XML\\x3\\x3\\xml3all.h"
#include "eq.hpp"

// ----
const TCHAR* ttitle = _T("EQ App");
// ----



unsigned long GetNextMP3Frame(unsigned char* dat, unsigned int ds, int& st)
{
	// Find a 0xFF with one with 3 
	for (unsigned int i = st; i < (ds - 1); i++)
	{
		unsigned char d0 = dat[i];
		unsigned char d1 = dat[i + 1];
		//		unsigned char d2 = dat[i + 2];
		//		unsigned char d3 = dat[i + 2];
		if (d0 != 0xFF)
			continue;
		if ((d1 & 0xE0) != 0xE0)
			continue;
		int Ly = ((d1 >> 1) & 3);

		if (Ly != 1)
			continue;

		st = i;

		unsigned int j = 0;
		memcpy(&j, dat + i, 4);
		return j;
	}
	return 0;
}


void ReadMP3Header(FILE* fp, int& sr, int& nCh)
{
	int CP = ftell(fp);
	std::vector<unsigned char> dat(100000);
	size_t ds = fread(dat.data(), 1, 100000, fp);
	fseek(fp, CP, SEEK_SET);

	if (ds == 0)
		return;


	int st = 0;

	for (int jF = 0; jF < 1; jF++)
	{
		size_t hdr = GetNextMP3Frame(dat.data(), (unsigned int)ds, st);
		if (hdr == 0)
			break; // duh

		// This is it
//		DWORD a = dat[st];
//		DWORD b = dat[st + 1];
		DWORD c = dat[st + 2];
		DWORD d = dat[st + 3];

		st++; // To go to next


		unsigned int SRT = (c & 0xC) >> 2;
		switch (SRT)
		{
			case 0:
				sr = 44100;
				break;
			case 1:
				sr = 48000;
				break;
			case 2:
				sr = 32000;
				break;
		}

		unsigned int CHT = (d & 0xC0) >> 6;
		switch (CHT)
		{
			case 3:
				nCh = 1;
				break;
			default:
				nCh = 2;
				break;
		}
	}

}

EQ::PARAMETRICEQ prx;
EQ::GRAPHICEQ prx2;


#include "f:\\tools\\uwl\\rw\\rw.hpp"
bool Ending = true;
using namespace std;
#include "..\\VSTX\\test\\wave.h"
#define MP3_BLOCK_SIZE 522
CONV c;
std::shared_ptr<WOUT> wout;
int SR = 48000;
int BS = 100;
int BR = 32;
std::recursive_mutex mu;

void StartMP3(std::wstring fi)
{
	using namespace std;
	FILE* fp = 0;
	_wfopen_s(&fp, fi.c_str(), L"rb");
	if (!fp)
		return;
	int nCh = 0;
	ReadMP3Header(fp, SR, nCh);
	if (!SR)
		return;
	BS = SR / 10;
	BR = 32;

	wout = 0;
	wout = make_shared<WOUT>();
	wout->Open(WAVE_MAPPER, SR, 16);

	vector<char> ww(1000);
	WAVEFORMATEX* w2 = (WAVEFORMATEX*)ww.data();
	WAVEFORMATEX& w = *w2;
	w.wFormatTag = WAVE_FORMAT_MPEGLAYER3;
	//		SR = 44100;


	w.cbSize = MPEGLAYER3_WFX_EXTRA_BYTES;
	w.wFormatTag = WAVE_FORMAT_MPEGLAYER3;
	w.nChannels = (WORD)nCh;
	w.nAvgBytesPerSec = 128 * (1024 / 8);  // not really used but must be one of 64, 96, 112, 128, 160kbps
	w.wBitsPerSample = 0;                  // MUST BE ZERO
	w.nBlockAlign = 1;                     // MUST BE ONE
	w.nSamplesPerSec = SR;              // 44.1kHz
	MPEGLAYER3WAVEFORMAT* mp3format = (MPEGLAYER3WAVEFORMAT*)&w;
	mp3format->fdwFlags = MPEGLAYER3_FLAG_PADDING_OFF;
	mp3format->nBlockSize = MP3_BLOCK_SIZE;             // voodoo value #1
	mp3format->nFramesPerBlock = 1;                     // MUST BE ONE
	mp3format->nCodecDelay = 1393;                      // voodoo value #2
	mp3format->wID = MPEGLAYER3_ID_MPEG;

	// This is a raw stream, based on ACM, so start decompression
	vector<char> e(1000);
	WAVEFORMATEX* wDest = (WAVEFORMATEX*)e.data();
	wDest->wFormatTag = WAVE_FORMAT_PCM;
	wDest->nSamplesPerSec = w.nSamplesPerSec;
	wDest->nChannels = 1;
	wDest->wBitsPerSample = 16;
	auto f = acmFormatSuggest(0, &w, wDest, 1000, ACM_FORMATSUGGESTF_NSAMPLESPERSEC | ACM_FORMATSUGGESTF_WFORMATTAG | ACM_FORMATSUGGESTF_NCHANNELS | ACM_FORMATSUGGESTF_WBITSPERSAMPLE);
	if (f)
		return;
	SR = wDest->nSamplesPerSec;
	BS = SR / 10;

	HACMSTREAM has = 0;
	f = acmStreamOpen(&has, 0, &w, wDest, 0, 0, 0, 0);
	if (f)
		return;

	DWORD InputSize = 3145728;
	InputSize = MP3_BLOCK_SIZE;//*2000;
	InputSize *= 100;

	DWORD OutputSize = 0;
	acmStreamSize(has, InputSize, &OutputSize, ACM_STREAMSIZEF_SOURCE);

	vector<char> in(InputSize + 1);
	vector<char> out(OutputSize + 1);
	ACMSTREAMHEADER ash = { 0 };
	ash.cbStruct = sizeof(ash);
	ash.pbSrc = (LPBYTE)in.data();
	ash.cbSrcLength = InputSize;
	ash.pbDst = (LPBYTE)out.data();
	ash.cbDstLength = OutputSize;
	f = acmStreamPrepareHeader(has, &ash, 0);

	vector<char> dx;
	size_t BytesRead = 0;
	Ending = false;
	if (true)
	{
		std::lock_guard<std::recursive_mutex> lg(mu);
		prx.Build(SR);
	
	}
	for (;;)
	{
		size_t FBR = fread(in.data(), 1, InputSize, fp);
		if (FBR == 0)
			break;
		BytesRead += FBR;
		ash.cbSrcLength = (DWORD)FBR;
		f = acmStreamConvert(has, &ash, 0);

		if (ash.cbSrcLengthUsed != ash.cbSrcLength)
			fseek(fp, (signed)(ash.cbSrcLengthUsed - ash.cbSrcLength), SEEK_CUR);

		auto ds = dx.size();
		dx.resize(ds + ash.cbDstLengthUsed);
		memcpy(dx.data() + ds, out.data(), ash.cbDstLengthUsed);


		for (;;)
		{
			if (dx.size() < (size_t)(BS * 2))
				break;

			if (Ending)
				break;
			// Convert
			int fsz = BS;
			vector<float> d2(fsz);
			c.f16t32((const short*)dx.data(), fsz, d2.data());

			vector<float> de = d2;

			if (true)
			{
				std::lock_guard<std::recursive_mutex> lg(mu);
				prx.Run(SR, d2.data(), fsz, de.data());
			}


			vector<short> d4(fsz * 2);
			memcpy(d2.data(), de.data(), fsz * sizeof(float));
			c.f32t16(d2.data(), fsz, d4.data());
			if (wout)
				wout->Write((const char*)d4.data(), fsz * 2);

			dx.erase(dx.begin(), dx.begin() + BS * 2);
		}

	}

	acmStreamUnprepareHeader(has, &ash, 0);
	acmStreamClose(has, 0);
	fclose(fp);

}


std::wstring OpenSingleFile(HWND hh, const wchar_t* filter, int fidx, const wchar_t* initf, const wchar_t* dir, const wchar_t* title)
{
	OPENFILENAME of = { 0 };
	of.lStructSize = sizeof(of);
	of.hwndOwner = hh;
	of.lpstrFilter = filter;
	of.lpstrInitialDir = dir;
	of.nFilterIndex = fidx;
	std::vector<wchar_t> fnx(10000);
	of.lpstrFile = fnx.data();
	if (initf)
		wcscpy_s(fnx.data(), 10000, initf);
	of.nMaxFile = 10000;
	of.lpstrTitle = title;
	of.Flags = OFN_FILEMUSTEXIST | OFN_EXPLORER;
	if (!GetOpenFileName(&of))
		return std::wstring(L"");
	return std::wstring(fnx.data());
}

//#include "curve.hpp"

CComPtr<ID2D1HwndRenderTarget> d;
CComPtr<ID2D1Factory> fa;


class M : public EQ::EQCALLBACK
{
public:

	M()
	{
	}

	virtual void RedrawRequest(EQ::EQ*)
	{
		InvalidateRect(MainWindow, 0, TRUE);
		UpdateWindow(MainWindow);
	}
	virtual void Dirty(EQ::EQ* e,bool)
	{

	}

};

std::shared_ptr<M> cxt;


LRESULT CALLBACK Main_DP(HWND hh, UINT mm, WPARAM ww, LPARAM ll)
{
	switch (mm)
	{
		case WM_CREATE:
		{
			MainWindow = hh;
			prx.SetWindow(hh);
			cxt = std::make_shared<M>();
			prx.AddCallback(cxt);
			XML3::XML x("eq.xml");
			prx.Unser(x.GetRootElement());



			x.Save();
			break;
		}

		case WM_COMMAND:
		{
			int LW = LOWORD(ww);
			UNREFERENCED_PARAMETER(LW);

			if (LW == 101)
			{

			}
			return 0;
		}

		case WM_CLOSE:
		{
			XML3::XML x("eq.xml");
			prx.Ser(x.GetRootElement());
			x.Save();
			DestroyWindow(hh);
			return 0;
		}

		case WM_KEYDOWN:
		case WM_SYSKEYDOWN:
		{
			prx.KeyDown(ww, ll);
			if (ww == VK_SPACE)
			{
				if (!Ending)
				{
					Ending = true;
				}
				else
				{
					auto thr = []()
					{
#ifdef _DEBUG
						wstring fi = L"G:\\mp3\\- 09- James Blunt - You Are Beautiful(1).mp3";
#else
						auto fi = OpenSingleFile(0, L"*.mp3", 0, 0, 0, 0);
#endif
						StartMP3(fi.c_str());
					};
					std::thread t(thr);
					t.detach();

				}
			}
			return 0;
		}
		case WM_MOUSEMOVE:
		{
			prx.MouseMove(ww, ll);
			return 0;
		}
case WM_MOUSEWHEEL:
{
	prx.MouseWheel(ww, ll);
	return 0;
}
case WM_LBUTTONDOWN:
		{
			prx.LeftDown(ww, ll);
			return 0;
		}
		case WM_RBUTTONDOWN:
		{
			prx.RightDown(ww, ll);
			return 0;
		}
		case WM_LBUTTONUP:
		{
			prx.LeftUp(ww, ll);
			return 0;
		}
		case WM_LBUTTONDBLCLK:
		{
			prx.LeftDoubleClick(ww, ll);
			return 0;
		}

		case WM_ERASEBKGND:
		{return 1;
		}
		case WM_PAINT:
		{
			PAINTSTRUCT ps = {};
			BeginPaint(hh, &ps);

			RECT rc;
			GetClientRect(hh, &rc);
			if (!fa)
				D2D1CreateFactory(D2D1_FACTORY_TYPE::D2D1_FACTORY_TYPE_MULTI_THREADED, &fa);
			if (!d)
			{
				D2D1_HWND_RENDER_TARGET_PROPERTIES hp;
				hp.hwnd = hh;
				hp.pixelSize.width = rc.right;
				hp.pixelSize.height = rc.bottom;
				d.Release();

				fa->CreateHwndRenderTarget(D2D1::RenderTargetProperties(), D2D1::HwndRenderTargetProperties(hh, D2D1::SizeU(rc.right - rc.left, rc.bottom - rc.top)), &d);
			}
			d->BeginDraw();
			if (true)
			{
				std::lock_guard<std::recursive_mutex> lg(mu);
				prx.Build(SR);
			}
			prx.Paint(fa,d, rc);
			d->EndDraw();
			EndPaint(hh, &ps);
			return 0;
		}

		case WM_SIZE:
		{
			if (!d)
				return 0;

			RECT rc;
			GetClientRect(hh, &rc);
			D2D1_SIZE_U u;
			u.width = rc.right;
			u.height = rc.bottom;
			d->Resize(u);
			return 0;
		}

		case WM_DESTROY:
		{
			PostQuitMessage(0);
			return 0;
		}
	}
	return DefWindowProc(hh, mm, ww, ll);
}



int __stdcall WinMain(HINSTANCE h, HINSTANCE, LPSTR, int)
{
	WSADATA wData;
	WSAStartup(MAKEWORD(2, 2), &wData);
	CoInitializeEx(0, COINIT_APARTMENTTHREADED);
	INITCOMMONCONTROLSEX icex = { 0 };
	icex.dwICC = ICC_LISTVIEW_CLASSES | ICC_DATE_CLASSES | ICC_WIN95_CLASSES;
	icex.dwSize = sizeof(icex);
	InitCommonControlsEx(&icex);
	InitCommonControls();
	hIcon1 = LoadIcon(0, IDI_APPLICATION);
	hAppInstance = h;

	WNDCLASSEX wClass = { 0 };
	wClass.cbSize = sizeof(wClass);

	wClass.style = CS_DBLCLKS | CS_HREDRAW | CS_VREDRAW | CS_PARENTDC;
	wClass.lpfnWndProc = (WNDPROC)Main_DP;
	wClass.hInstance = h;
	wClass.hIcon = hIcon1;
	wClass.hCursor = LoadCursor(0, IDC_ARROW);
	wClass.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	wClass.lpszClassName = _T("CLASS");
	wClass.hIconSm = hIcon1;
	RegisterClassEx(&wClass);

	MainWindow = CreateWindowEx(0,
		_T("CLASS"),
		ttitle,
		WS_OVERLAPPEDWINDOW | WS_CLIPSIBLINGS |
		WS_CLIPCHILDREN, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
		0, 0, h, 0);

	ShowWindow(MainWindow, SW_SHOW);


	MSG msg;

	HACCEL acc = LoadAccelerators(h, L"MENU_1");
	while (GetMessage(&msg, 0, 0, 0))
	{
		if (TranslateAccelerator(msg.hwnd, acc, &msg))
			continue;
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	ExitProcess(0);
}