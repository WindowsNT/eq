//

extern "C" {
#include "biquad.h"
}

#ifdef TURBO_PLAY

#include "..\\\AED\\\fft.hpp"
#include "..\\fft-real-master\\\FFTReal.h"
#include "..\\DSPFilters\\DSPFilters\\include\\DspFilters\\Elliptic.h"
#include "..\\DSPFilters\\DSPFilters\\include\\DspFilters\\ChebyshevII.h"
#include "..\\DSPFilters\\DSPFilters\\include\\DspFilters\\ChebyshevI.h"
#include "..\\DSPFilters\\DSPFilters\\include\\DspFilters\\Butterworth.h"
#else
#include "alldspfilters.hpp"
#endif

#pragma warning(disable:4100)

namespace EQ
{
#ifdef _WIN64
	typedef signed long long ssize_t;
#else
	typedef signed long ssize_t;
#endif

	class EQ;
	using namespace std;


	class EQCALLBACK
	{
	public:

		virtual void RedrawRequest(EQ* pr) = 0;
		virtual void Dirty(EQ* e,bool) = 0;

	};

	class MMCB : public EQCALLBACK
	{
	public:
		HWND hC = 0;

		virtual void RedrawRequest(EQ* p)
		{
			if (!IsWindow(hC))
				return;
			InvalidateRect(hC, 0, TRUE);
			UpdateWindow(hC);
		}
		virtual void Dirty(EQ* q,bool)
		{

		}
	};


	class EQ
	{
	public:
		HWND hParent;

		void SetWindow(HWND hh)
		{
			hParent = hh;
		}


		struct ASKTEXT
		{
			const wchar_t* ti;
			const wchar_t* as;
			wchar_t* re;
			int P;
			wstring* re2;
			int mx = 1000;
		};

		static INT_PTR CALLBACK A_DP(HWND hh, UINT mm, WPARAM ww, LPARAM ll)
		{
			static ASKTEXT* as = 0;
			switch (mm)
			{
				case WM_INITDIALOG:
				{
					as = (ASKTEXT*)ll;
					SetWindowText(hh, as->ti);
					if (as->P != 2)
					{
						SetWindowText(GetDlgItem(hh, 101), as->as);
						if (as->re)
							SetWindowText(GetDlgItem(hh, 102), as->re);
						if (as->re2)
							SetWindowText(GetDlgItem(hh, 102), as->re2->c_str());
					}
					else
						SetWindowText(GetDlgItem(hh, 701), as->as);
					if (as->P == 1)
					{
						auto w = GetWindowLongPtr(GetDlgItem(hh, 102), GWL_STYLE);
						w |= ES_PASSWORD;
						SetWindowLongPtr(GetDlgItem(hh, 102), GWL_STYLE, w);
					}
					return true;
				}
				case WM_COMMAND:
				{
					if (LOWORD(ww) == IDOK)
					{
						wchar_t re1[1000] = { 0 };
						wchar_t re2[1000] = { 0 };
	//					MessageBox(0, L"foo", 0, 0);
						if (as->P == 2)
						{
							GetWindowText(GetDlgItem(hh, 101), re1, 1000);
							GetWindowText(GetDlgItem(hh, 102), re2, 1000);
							if (wcscmp(re1, re2) != 0 || wcslen(re1) == 0)
							{
								SetWindowText(GetDlgItem(hh, 101), L"");
								SetWindowText(GetDlgItem(hh, 102), L"");
								SetFocus(GetDlgItem(hh, 101));
								return 0;
							}
							wcscpy_s(as->re, 1000, re1);
							EndDialog(hh, IDOK);
							return 0;
						}

						if (as->re2)
						{
							int lex = GetWindowTextLength(GetDlgItem(hh, 102));
							vector<wchar_t> re(lex + 100);
							GetWindowText(GetDlgItem(hh, 102), re.data(), lex + 100);
							*as->re2 = re.data();
							EndDialog(hh, IDOK);
						}
						else
						{
							GetWindowText(GetDlgItem(hh, 102), as->re, as->mx);
							EndDialog(hh, IDOK);
						}
						return 0;
					}
					if (LOWORD(ww) == IDCANCEL)
					{
						EndDialog(hh, IDCANCEL);
						return 0;
					}
				}
			}
			return 0;
		}

		bool AskText(HWND hh, const TCHAR* ti, const TCHAR* as, TCHAR* re, wstring* re2 = 0, int mxt = 1000)
		{
			const char* res = "\xC4\x0A\xCA\x90\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x6D\x01\x3E\x00\x00\x00\x00\x00\x00\x00\x0A\x00\x54\x00\x61\x00\x68\x00\x6F\x00\x6D\x00\x61\x00\x00\x00\x01\x00\x00\x50\x00\x00\x00\x00\x0A\x00\x09\x00\x1C\x01\x1A\x00\x65\x00\xFF\xFF\x82\x00\x00\x00\x00\x00\x00\x00\x80\x00\x81\x50\x00\x00\x00\x00\x0A\x00\x29\x00\x1D\x01\x0F\x00\x66\x00\xFF\xFF\x81\x00\x00\x00\x00\x00\x00\x00\x00\x03\x01\x50\x00\x00\x00\x00\x2F\x01\x16\x00\x32\x00\x0E\x00\x01\x00\xFF\xFF\x80\x00\x4F\x00\x4B\x00\x00\x00\x00\x00\x00\x00\x00\x03\x01\x50\x00\x00\x00\x00\x2F\x01\x29\x00\x32\x00\x0E\x00\x02\x00\xFF\xFF\x80\x00\x43\x00\x61\x00\x6E\x00\x63\x00\x65\x00\x6C\x00\x00\x00\x00\x00";
			ASKTEXT a = { ti,as,re,0,re2,mxt };
			return (DialogBoxIndirectParam(GetModuleHandle(0), (LPCDLGTEMPLATEW)res, hh, A_DP, (LPARAM)&a) == IDOK);
		}


		HCURSOR ArrowCursor = LoadCursor(0, IDC_ARROW);
		HCURSOR ResizeCursor = LoadCursor(0, IDC_SIZENS);
		
		vector<shared_ptr<EQCALLBACK>> cb;
		D2D1_COLOR_F bg = { 0.1f,0.1f,0.1f,1.0f };
		D2D1_COLOR_F whitecolor = { 1.0f,1.0f,1.0f,1.0f };
		D2D1_COLOR_F graycolor = { 0.5f,0.5f,0.5f,0.6f };
		D2D1_COLOR_F yellowcolor = { 0.9f,0.9f,0.1f,0.5f };
		D2D1_COLOR_F selectcolor = { 0.0f,9.0f,3.0f,1.0f };
		D2D1_COLOR_F blackcolor = { 0.0f,0.0f,0.0f,1.0f };

		CComPtr<IDWriteFactory> WriteFactory;
		CComPtr<IDWriteTextFormat> Text;
		CComPtr<ID2D1SolidColorBrush> BGBrush;
		CComPtr<ID2D1SolidColorBrush> WhiteBrush;
		CComPtr<ID2D1SolidColorBrush> GrayBrush;
		CComPtr<ID2D1SolidColorBrush> YellowBrush;
		CComPtr<ID2D1SolidColorBrush> SelectBrush;
		CComPtr<ID2D1SolidColorBrush> BlackBrush;


		template <typename T = float> bool InRect(D2D1_RECT_F & r, T x, T y)
		{
			if (x >= r.left && x <= r.right && y >= r.top && y <= r.bottom)
				return true;
			return false;
		}


		D2D1_RECT_F FromR(RECT rcc)
		{
			D2D1_RECT_F r;
			r.left = (FLOAT)rcc.left;
			r.top = (FLOAT)rcc.top;
			r.right = (FLOAT)rcc.right;
			r.bottom = (FLOAT)rcc.bottom;
			return r;
		}
		void RemoveCallbacks()
		{
			cb.clear();
		}
		void AddCallback(shared_ptr<EQCALLBACK> cx)
		{
			cb.push_back(cx);
		}
		void DestroyBrushes()
		{
			BGBrush = 0;
			WhiteBrush = 0;
			GrayBrush = 0;
			YellowBrush = 0;
			BlackBrush = 0;
			SelectBrush = 0;
			Text = 0;
			WriteFactory = 0;
		}
		CComPtr<ID2D1SolidColorBrush> GetD2SolidBrush(ID2D1RenderTarget* p, D2D1_COLOR_F cc)
		{
			CComPtr<ID2D1SolidColorBrush> b = 0;
			p->CreateSolidColorBrush(cc, &b);
			return b;
		}
		void CreateBrushes(ID2D1RenderTarget* p, bool F = false)
		{
			if (WhiteBrush && !F)
				return; // OK

			SelectBrush = GetD2SolidBrush(p, selectcolor);
			WhiteBrush = GetD2SolidBrush(p, whitecolor);
			BlackBrush = GetD2SolidBrush(p, blackcolor);
			GrayBrush = GetD2SolidBrush(p, graycolor);
			YellowBrush = GetD2SolidBrush(p, yellowcolor);
			BGBrush = GetD2SolidBrush(p, bg);
			DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(IDWriteFactory), (IUnknown * *)& WriteFactory);

			LOGFONT lf;
			GetObject(GetStockObject(DEFAULT_GUI_FONT), sizeof(lf), &lf);
			DWRITE_FONT_STYLE fst = DWRITE_FONT_STYLE_NORMAL;
			if (lf.lfItalic)
				fst = DWRITE_FONT_STYLE_ITALIC;
			DWRITE_FONT_STRETCH fsr = DWRITE_FONT_STRETCH_NORMAL;
			FLOAT fs = (FLOAT)abs(lf.lfHeight);
			WriteFactory->CreateTextFormat(lf.lfFaceName, 0, lf.lfWeight > 500 ? DWRITE_FONT_WEIGHT_BOLD : DWRITE_FONT_WEIGHT_NORMAL, fst, fsr, fs, L"", &Text);
		}

		void Redraw()
		{
			for (auto cc : cb)
				cc->RedrawRequest(this);
		}

		void Dirty(bool x)
		{
			for (auto cc : cb)
				cc->Dirty(this,x);
		}


		virtual void Paint(ID2D1Factory*fact, ID2D1RenderTarget* r, RECT rc) = 0;
		virtual void LeftDown(WPARAM ww, LPARAM ll) = 0;
		virtual void RightDown(WPARAM ww, LPARAM ll) = 0;
		virtual void LeftUp(WPARAM ww, LPARAM ll) = 0;
		virtual void MouseMove(WPARAM ww, LPARAM ll) = 0;
		virtual void MouseWheel(WPARAM ww, LPARAM ll) = 0;
		virtual void KeyDown(WPARAM ww, LPARAM ll) = 0;
		virtual void LeftDoubleClick(WPARAM ww,LPARAM ll) = 0;
		D2D1_RECT_F rc = {};
		virtual void Ser(XML3::XMLElement& e) = 0;
		virtual void Unser(XML3::XMLElement& e) = 0;

		virtual void Prepare(int SR) = 0;
		virtual void Run(int SR, float* in, int ns, float* out) = 0;
		virtual bool Run2(int SR, int nch,float** in, int ns, float** out) = 0;
		virtual void Build(int SR) = 0;

		int Max = 22000;

		int XMode = 1; // Logarithmic

		float dBFSToPercent(float dBFS = 0)
		{
			return 100.0f * (pow(10.0f, (dBFS / 20.0f)));
		}
		float PercentTodBFS(float P)
		{
			return 20.0f * log10(P / 100.0f);
		}


		float dB2V(float dB)
		{
			if (dB == 0)
				return 1.0f;
			if (dB > 0)
			{
				float p = dBFSToPercent(-dB) / 100.0f;
				float V = 1 + (1.0f - p);
				return V;
			}
			float p = dBFSToPercent(dB) / 100.0f;
			float V = p;
			return V;
		}

		float V2Y(float V)
		{
			// In 2, height
			// in V,  ?
			float he = rc.bottom - rc.top;
			return he - (he)*V / 2.0f;
		}


		float Y2V(float Y)
		{
			// In V=2, height
			// ?   , y
			return 2.0f - 2.0f * Y / (rc.bottom - rc.top);
		}

		float V2dB(float V,bool NeedNZ)
		{
			float dBX = 0;
			if (V < 1)
				dBX = PercentTodBFS((V) * 100);
			if (V > 1)
				dBX = -PercentTodBFS((2.0f - V) * 100);
			if (dBX == 0.0 && NeedNZ)
				dBX = 0.001f;
			return dBX;
		}

		wstring HZString(float z)
		{
			wchar_t a[100] = { 0 };
			if (z < 1000)
				swprintf_s(a, 100, L"%.0f Hz", z);
			else
				swprintf_s(a, 100, L"%.0f KHz", z / 1000.0f);
			return a;
		}

		wstring HZString(float z1, float z2)
		{
			wchar_t a[100] = { 0 };
			if (z1 < 1000 && z2 < 1000)
				swprintf_s(a, 100, L"%.0f - %.0f Hz", z1, z2);
			else
				if (z1 < 1000 && z2 >= 1000)
					swprintf_s(a, 100, L"%.0f Hz - %.0f KHz", z1, z2 / 1000);
				else
					swprintf_s(a, 100, L"%.0f - %.0f KHz", z1 / 1000.0f, z2 / 1000.0f);
			return a;
		}

		void ShowInDialog(HWND hh)
		{
			struct Z
			{
				EQ* c = 0;
				CComPtr<ID2D1HwndRenderTarget> d;
				CComPtr<ID2D1Factory> fa;
			};

			const char* res = "\x01\x00\xFF\xFF\x00\x00\x00\x00\x00\x00\x00\x00\xC8\x08\xCF\x80\x00\x00\x00\x00\x00\x00\xC5\x02\x95\x01\x00\x00\x00\x00\x00\x00\x08\x00\x90\x01\x00\x01\x4D\x00\x53\x00\x20\x00\x53\x00\x68\x00\x65\x00\x6C\x00\x6C\x00\x20\x00\x44\x00\x6C\x00\x67\x00\x00\x00";

			auto dp = [](HWND hh, UINT mm, WPARAM ww, LPARAM ll) -> INT_PTR
			{
				Z* z = (Z*)GetWindowLongPtr(hh, GWLP_USERDATA);
				EQ* c = 0;
				if (z)
					c = z->c;

				switch (mm)
				{
					case WM_INITDIALOG:
					{
						SetWindowLongPtr(hh, GWLP_USERDATA, ll);
						z = (Z*)GetWindowLongPtr(hh, GWLP_USERDATA);
						c = z->c;

						c->SetWindow(hh);
						auto mmd = std::make_shared<MMCB>();
						mmd->hC = hh;
						c->AddCallback(mmd);
						return true;
					}
					case WM_CLOSE:
					{
						EndDialog(hh, 0);
						return 0;
					}
					case WM_KEYDOWN:
					case WM_SYSKEYDOWN:
					{
						c->KeyDown(ww, ll);
						return 0;
					}

					case WM_MOUSEMOVE:
					{
						c->MouseMove(ww, ll);
						return 0;
					}
					case WM_MOUSEWHEEL:
					{
						c->MouseWheel(ww, ll);
						return 0;
					}
					case WM_LBUTTONDOWN:
					{
						c->LeftDown(ww, ll);
						return 0;
					}
					case WM_RBUTTONDOWN:
					{
						c->RightDown(ww, ll);
						return 0;
					}
					case WM_LBUTTONUP:
					{
						c->LeftUp(ww, ll);
						return 0;
					}
					case WM_LBUTTONDBLCLK:
					{
						c->LeftDoubleClick(ww, ll);
						return 0;
					}
					case WM_ERASEBKGND:
					{
						return 1;
					}
					case WM_PAINT:
					{
						PAINTSTRUCT ps;
						BeginPaint(hh, &ps);

						RECT rc;
						GetClientRect(hh, &rc);
						if (!z->fa)
							D2D1CreateFactory(D2D1_FACTORY_TYPE::D2D1_FACTORY_TYPE_MULTI_THREADED, &z->fa);
						if (!z->d)
						{
							//				D2D1_RENDER_TARGET_PROPERTIES p;
							D2D1_HWND_RENDER_TARGET_PROPERTIES hp;
							hp.hwnd = hh;
							hp.pixelSize.width = rc.right;
							hp.pixelSize.height = rc.bottom;
							z->d.Release();

							z->fa->CreateHwndRenderTarget(D2D1::RenderTargetProperties(), D2D1::HwndRenderTargetProperties(hh, D2D1::SizeU(rc.right - rc.left, rc.bottom - rc.top)), &z->d);
						}
						z->d->BeginDraw();
						z->c->Paint(z->fa, z->d, rc);
						[[maybe_unused]] auto hr = z->d->EndDraw();
						EndPaint(hh, &ps);
						return 0;
					}

					case WM_SIZE:
					{
						if (!z->d)
							return 0;

						RECT rc;
						GetClientRect(hh, &rc);
						D2D1_SIZE_U u;
						u.width = rc.right;
						u.height = rc.bottom;
						z->d->Resize(u);
						return 0;
					}

				}
				return 0;
			};
			Z z;
			z.c = this;
			DialogBoxIndirectParam(0, (LPCDLGTEMPLATE)res, hh ? hh : hParent, dp, (LPARAM)& z);
			DestroyBrushes();
		}

		void PaintDBLines(ID2D1RenderTarget* r)
		{
			if (true)
			{
				// Lines at 3,6,12,15 dB
				for (float dbs : { -15.0f, -12.0f, -6.0f, -3.0f, 0.0f, 3.0f, 6.0f, 12.0f, 15.0f })
				{
					float V = (float)dB2V((float)dbs);
					float y = (float)V2Y(V);

					D2D1_POINT_2F p1, p2;

					p1.y = y;
					p1.x = rc.left;
					p2.x = rc.right;
					p2.y = y;

					wchar_t t[100] = { 0 };
					swprintf_s(t, L"%.f", dbs);
					Text->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
					Text->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_CENTER);
#ifdef TURBO_PLAY
					auto rrs = AE::MeasureString(WriteFactory, Text, t, (UINT32)wcslen(t));
#else
					std::tuple<float, float> rrs = std::make_tuple(30.0f, 30.0f);
#endif

					D2D1_RECT_F ly = {};
					ly.left = p1.x + 2;
					ly.top = p1.y - 25;
					ly.bottom = ly.top + 50;
					p1.x += std::get<0>(rrs) + 3;
					ly.right = p1.x + 1;

					D2D1_RECT_F ly2 = {};
					ly2.left = p2.x - std::get<0>(rrs);
					ly2.right = p2.x;
					ly2.top = p1.y - 25;
					ly2.bottom = ly.top + 50;
					p2.x -= std::get<0>(rrs) + 3;


					if (dbs == 0)
						r->DrawLine(p1, p2, SelectBrush);
					else
						r->DrawLine(p1, p2, GrayBrush);

					r->DrawTextW(t, (UINT32)wcslen(t), Text, ly, WhiteBrush);
					r->DrawTextW(t, (UINT32)wcslen(t), Text, ly2, WhiteBrush);

				}
			}

		}

	};

	class BAND
	{
	public:
		D2D1_RECT_F r = {};
		D2D1_RECT_F r2 = {};
		sf_biquad_state_st state;
		float from, to;
		float V = 1.0f;

		virtual void Ser(XML3::XMLElement& e)
		{
			e.vv("f").SetValueFloat(from);
			e.vv("t").SetValueFloat(to);
			e.vv("v").SetValueFloat(V);
		}
		virtual void Unser(XML3::XMLElement& e)
		{
			from = e.vv("f").GetValueFloat();
			to = e.vv("t").GetValueFloat();
			V = e.vv("v").GetValueFloat();
		}

	};

	typedef Dsp::SimpleFilter<Dsp::Butterworth::LowPass<100>, 1> ButtLP;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevI::LowPass<100>, 1> Che1LP;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevII::LowPass<100>, 1> Che2LP;
	typedef Dsp::SimpleFilter<Dsp::Elliptic::LowPass<100>, 1> EllLP;

	typedef Dsp::SimpleFilter<Dsp::Butterworth::HighPass<100>, 1> ButtHP;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevI::HighPass<100>, 1> Che1HP;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevII::HighPass<100>, 1> Che2HP;
	typedef Dsp::SimpleFilter<Dsp::Elliptic::HighPass<100>, 1> EllHP;

	typedef Dsp::SimpleFilter<Dsp::Butterworth::LowShelf<100>, 1> ButtLPs;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevI::LowShelf<100>, 1> Che1LPs;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevII::LowShelf<100>, 1> Che2LPs;

	typedef Dsp::SimpleFilter<Dsp::Butterworth::HighShelf<100>, 1> ButtHPs;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevI::HighShelf<100>, 1> Che1HPs;
	typedef Dsp::SimpleFilter<Dsp::ChebyshevII::HighShelf<100>, 1> Che2HPs;

	struct PFILT
	{
		D2D1_RECT_F r = {};
		int SpecialType = 0; // 0 Biq, 1 Butt, 2 Che I , 3 Che II, 4 Elliptic
		int SpecialFilter = 0; // 1 LC, 2 HC, 0 Middle
		bool S = false;
		bool A = false;
		float fr = 0.0f;
		float Q = 1;
		float dB = 0;
		int Order = 3;
		float ripple = 0.5f;
		float rolloff = 12.0f;
		int Type = 0; // 0 LP, 1 HP , 2 LS , 3 HS, 4 Notch, 5 Peaking

		void Ser(XML3::XMLElement& e)
		{
			e.vv("a").SetValueInt(A);
			e.vv("Type").SetValueInt(Type);
			e.vv("st").SetValueInt(SpecialType);
			e.vv("sf").SetValueInt(SpecialFilter);
			e.vv("o").SetValueInt(Order);
			e.vv("fr").SetValueFloat(fr);
			e.vv("Q").SetValueFloat(Q);
			e.vv("dB").SetValueFloat(dB);
			e.vv("ripple").SetValueFloat(ripple);
			e.vv("rolloff").SetValueFloat(rolloff);

		}
		void Unser(XML3::XMLElement& e)
		{
			fr = e.vv("fr").GetValueFloat();
			Q = e.vv("Q").GetValueFloat();
			dB = e.vv("dB").GetValueFloat();
			ripple = e.vv("ripple").GetValueFloat(0.5f);
			rolloff = e.vv("rolloff").GetValueFloat(24.0f);
			if (ripple <= 0.0)
				ripple = 0.5f;
			if (rolloff <= 0.0 || rolloff >= 18.0f)
				rolloff = 1.0f;
			Type = e.vv("Type").GetValueInt();
			Order = e.vv("o").GetValueInt(3);
			if (Order == 0)
				Order = 3;
			SpecialType = e.vv("st").GetValueInt();
			SpecialFilter = e.vv("sf").GetValueInt();
			A = (bool)e.vv("a").GetValueInt();
		}

		bool operator <(const PFILT& f2)
		{
			if (fr < f2.fr)
				return true;
			return false;
		}

		sf_biquad_state_st st = { 0 };
		std::shared_ptr<void> sf = 0;


		void Build(int SR)
		{
			sf = 0;
			switch (Type)
			{
				case 0: // Low Pass
				{
					sf_lowpass(&st, SR, fr,0);
					
					// And the DSP
					if (SpecialType == 0)
					{
						// Keep Biquad
					}
					if (SpecialType == 1)
					{
						auto sf2 = std::make_shared<ButtLP>();
						sf2->setup(Order, SR, fr);
						sf = sf2;
					}
					if (SpecialType == 2)
					{
						auto sf2 = std::make_shared<Che1LP>();
						sf2->setup(Order, SR, fr,ripple);
						sf = sf2;
					}
					if (SpecialType == 3)
					{
						auto sf2 = std::make_shared<Che2LP>();
						sf2->setup(Order, SR, fr, ripple);
						sf = sf2;
					}
					if (SpecialType == 4)
					{
						auto sf2 = std::make_shared<EllLP>();
						sf2->setup(Order, SR, fr, ripple,rolloff);
						sf = sf2;
					}

					break;
				}
				case 1: // High Pass
				{
					sf_highpass(&st, SR, fr, 0);

					// And the DSP
					if (SpecialType == 0)
					{
						// Keep Biquad
					}
					if (SpecialType == 1)
					{
						auto sf2 = std::make_shared<ButtHP>();
						sf2->setup(Order, SR, fr);
						sf = sf2;
					}
					if (SpecialType == 2)
					{
						auto sf2 = std::make_shared<Che1HP>();
						sf2->setup(Order, SR, fr, ripple);
						sf = sf2;
					}
					if (SpecialType == 3)
					{
						auto sf2 = std::make_shared<Che2HP>();
						sf2->setup(Order, SR, fr, ripple);
						sf = sf2;
					}
					if (SpecialType == 4)
					{
						auto sf2 = std::make_shared<EllHP>();
						sf2->setup(Order, SR, fr, ripple, rolloff);
						sf = sf2;
					}

					break;
				}
				case 2: // Low Shelf
				{
					sf_lowshelf(&st, SR, fr, Q,dB);

					// And the DSP
					if (SpecialType == 6)
					{
						auto sf2 = std::make_shared<ButtLPs>();
						sf2->setup(Order, SR, fr,dB);
						sf = sf2;
					}
					if (SpecialType == 7)
					{
						auto sf2 = std::make_shared<Che1LPs>();
						sf2->setup(Order, SR, fr, dB,ripple);
						sf = sf2;
					}
					if (SpecialType == 8)
					{
						auto sf2 = std::make_shared<Che2LPs>();
						sf2->setup(Order, SR, fr,dB,ripple);
						sf = sf2;
					}

					break;
				}
				case 3: // High Shelf
				{
					sf_highshelf(&st, SR, fr, Q,dB);

					// And the DSP
					if (SpecialType == 6)
					{
						auto sf2 = std::make_shared<ButtHPs>();
						sf2->setup(Order, SR, fr, dB);
						sf = sf2;
					}
					if (SpecialType == 7)
					{
						auto sf2 = std::make_shared<Che1HPs>();
						sf2->setup(Order, SR, fr, dB, ripple);
						sf = sf2;
					}
					if (SpecialType == 8)
					{
						auto sf2 = std::make_shared<Che2LPs>();
						sf2->setup(Order, SR, fr, dB, ripple);
						sf = sf2;
					}

					break;
				}
				case 4: // Notch
				{
					sf_notch(&st, SR, fr, Q);
					break;
				}
				case 5: // peaking
				{
					sf_peaking(&st, SR, fr, Q,dB);
					break;
				}
			}
		}



	};

	class PARAMETRICEQ : public EQ
	{
	public:

		int rad = 5;

		PARAMETRICEQ() : EQ()
		{
			filters.resize(2);
			filters[0].A = 0;
			filters[0].fr = 200;
			filters[0].Type = 1;
			filters[1].fr = 4000;
			filters[0].SpecialFilter = 1;
			filters[0].SpecialType = 4; // Elliptic 
			filters[1].SpecialFilter = 2;
			filters[1].SpecialType = 4; // 
		}

		vector<PFILT> filters;

		void PaintTop(ID2D1RenderTarget* r, RECT rrc)
		{
			wchar_t y[1000] = {};
			D2D1_RECT_F r2 = {};
			r2.left = (FLOAT)rrc.left;
			r2.top = (FLOAT)rrc.top;
			r2.right = (FLOAT)rrc.right;
			r2.bottom = r2.top + 25;
			Text->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
			Text->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_CENTER);

			swprintf_s(y, 1000, L"Parametric Equalizer");
			for (auto& f : filters)
			{
				if (f.S)
				{
					if (f.SpecialFilter == 0) // Peaking
					{
						float f1 = f.fr * (sqrt(1.0f + (1.0f / (4.0f * f.Q * f.Q))) - 1.0f / (2.0f * f.Q));
						float f2 = f.fr * (sqrt(1.0f + (1.0f / (4.0f * f.Q * f.Q))) + 1.0f / (2.0f * f.Q));

						swprintf_s(y, 1000, L"Peak filter %.1f Hz, %.1f dB, Q = %.1f (Range %.1f - %.1f)", f.fr, f.dB, f.Q,f1,f2);
					}
					else
					{
						if (f.SpecialFilter == 1) // Low Cut
						{
							if (f.SpecialType == 0)
								swprintf_s(y, 1000, L"Biquad low cut %.1f Hz", f.fr);
							if (f.SpecialType == 1)
								swprintf_s(y, 1000, L"Butterworth low cut order %i %.1f", f.Order, f.fr);
							if (f.SpecialType == 2)
								swprintf_s(y, 1000, L"Chebyshev I low cut order %i %.1f Hz ripple %.1f dB", f.Order, f.fr,f.ripple);
							if (f.SpecialType == 3)
								swprintf_s(y, 1000, L"Chebyshev II low cut order %i %.1f Hz ripple %.1f dB", f.Order, f.fr, f.ripple);
							if (f.SpecialType == 4)
								swprintf_s(y, 1000, L"Elliptic low cut order %i %.1f Hz ripple %.1f dB rolloff %.1f", f.Order, f.fr,f.ripple,f.rolloff);


							if (f.SpecialType == 5)
								swprintf_s(y, 1000, L"Biquad low shelf order %i %.1f dB %.1f Hz", f.Order,  f.dB,f.fr);
							if (f.SpecialType == 6)
								swprintf_s(y, 1000, L"Butterworth low shelf order %i %.1f dB %.1f Hz", f.Order, f.dB, f.fr);
							if (f.SpecialType == 7)
								swprintf_s(y, 1000, L"Chebyshev I low shelf order %i %.1f dB %.1f Hz ripple %.1f dB", f.Order, f.dB, f.fr, f.ripple);
							if (f.SpecialType == 8)
								swprintf_s(y, 1000, L"Chebyshev II low shelf order %i %.1f dB %.1f Hz ripple %.1f dB", f.Order, f.dB, f.fr, f.ripple);

						}
						if (f.SpecialFilter == 2) // High Cut
						{
							if (f.SpecialType == 0)
								swprintf_s(y, 1000, L"Biquad high cut %.1f Hz", f.fr);
							if (f.SpecialType == 1)
								swprintf_s(y, 1000, L"Butterworth high cut order %i %.1f", f.Order, f.fr);
							if (f.SpecialType == 2)
								swprintf_s(y, 1000, L"Chebyshev I high cut order %i %.1f Hz ripple %.1f dB", f.Order, f.fr, f.ripple);
							if (f.SpecialType == 3)
								swprintf_s(y, 1000, L"Chebyshev II high cut order %i %.1f Hz ripple %.1f dB", f.Order, f.fr, f.ripple);
							if (f.SpecialType == 4)
								swprintf_s(y, 1000, L"Elliptic high cut order %i %.1f Hz ripple %.1f dB rolloff %.1f", f.Order, f.fr, f.ripple, f.rolloff);


							if (f.SpecialType == 5)
								swprintf_s(y, 1000, L"Biquad high shelf order %i %.1f dB %.1f Hz", f.Order, f.dB, f.fr);
							if (f.SpecialType == 6)
								swprintf_s(y, 1000, L"Butterworth high shelf order %i %.1f dB %.1f Hz", f.Order, f.dB, f.fr);
							if (f.SpecialType == 7)
								swprintf_s(y, 1000, L"Chebyshev I high shelf order %i %.1f dB %.1f Hz ripple %.1f dB", f.Order, f.dB, f.fr, f.ripple);
							if (f.SpecialType == 8)
								swprintf_s(y, 1000, L"Chebyshev II high shelf order %i %.1f dB %.1f Hz ripple %.1f dB", f.Order, f.dB, f.fr, f.ripple);

						}
					}
				}
			}
			r->DrawTextW(y, (UINT32)wcslen(y), Text, r2, WhiteBrush);
		}


		void CreateTheResponse(PFILT& f, CComPtr<ID2D1GeometrySink> pgs)
		{
			auto abs = [](Dsp::complex_t c) -> double
			{
				auto a = sqrt(c.real() * c.real() + c.imag() * c.imag());
				return 10.0f * log10(a);
			};

			auto& b = f;

			// Paint response
			vector<D2D1_POINT_2F> pts;
			vector<float> freqs;
			for (int x2 = 1; x2 <= (int)rc.right; x2++)
				freqs.push_back(X2Freqr((float)x2));

			auto fill = [&](Dsp::Cascade* fx)
			{
				for (auto& ff : freqs)
				{
					D2D1_POINT_2F pp;
					pp.x = Freq2X(ff);
					float db = (float)abs(fx->response(ff / sr));
					pp.y = V2Y(dB2V(db));
					pts.push_back(pp);
				}
			};

			if (b.SpecialType == 1)
			{
				std::shared_ptr<ButtHP> fx;
				fx = std::static_pointer_cast<ButtHP>(f.sf);
				fill(fx.get());
			}
			if (b.SpecialType == 2)
			{
				std::shared_ptr<Che1HP> fx;
				fx = std::static_pointer_cast<Che1HP>(f.sf);
				fill(fx.get());
			}
			if (b.SpecialType == 3)
			{
				std::shared_ptr<Che2HP> fx;
				fx = std::static_pointer_cast<Che2HP>(f.sf);
				fill(fx.get());
			}
			if (b.SpecialType == 4)
			{
				std::shared_ptr<EllHP> fx;
				fx = std::static_pointer_cast<EllHP>(f.sf);
				fill(fx.get());
			}
			if (b.SpecialType == 6)
			{
				std::shared_ptr<ButtLPs> fx;
				fx = std::static_pointer_cast<ButtLPs>(f.sf);
				fill(fx.get());
			}
			if (b.SpecialType == 7)
			{
				std::shared_ptr<Che1LPs> fx;
				fx = std::static_pointer_cast<Che1LPs>(f.sf);
				fill(fx.get());
			}
			if (b.SpecialType == 8)
			{
				std::shared_ptr<Che2LPs> fx;
				fx = std::static_pointer_cast<Che2LPs>(f.sf);
				fill(fx.get());
			}

			if (pts.size() > 1)
			{
				D2D1_POINT_2F fb = { pts[0].x,pts[0].y };
				pgs->BeginFigure(fb, D2D1_FIGURE_BEGIN_HOLLOW);

				pgs->AddLines(pts.data() + 1, (UINT32)(pts.size() - 1));
				pgs->EndFigure(D2D1_FIGURE_END_OPEN);
			}
		}

		virtual void Paint(ID2D1Factory* fact, ID2D1RenderTarget* r, RECT rrc)
		{
			wchar_t t[1000] = { 0 };
			CreateBrushes(r);
			rc = FromR(rrc);
			auto rr = rc;
			r->FillRectangle(rc, BGBrush);

			PaintDBLines(r);

			CComPtr<ID2D1PathGeometry> pg = 0;
			CComPtr<ID2D1GeometrySink> pgs = 0;
			fact->CreatePathGeometry(&pg);
			pg->Open(&pgs);

//			float xmiddle = (rc.right - rc.left) / 2.0f + rc.left;
			float ymiddle = (rc.bottom - rc.top) / 2.0f + rc.top;

			// Filters
			for (size_t i = 0; i < filters.size(); i++)
			{
				auto& f = filters[i];

				auto x = Freq2X(f.fr);
				auto y = V2Y(dB2V(f.dB));
				D2D1_ELLIPSE e;
				e.point.x = x;
				e.point.y = y;
				e.radiusX = (float)rad*2;
				e.radiusY = (float)rad*2;
				if (!f.A)
				{
					r->FillEllipse(e, GrayBrush);
				}
				else
				{
					if (f.S)
						r->FillEllipse(e, SelectBrush);
					else
						r->FillEllipse(e, WhiteBrush);
				}

				f.r.left = e.point.x - rad*2;
				f.r.top = e.point.y - rad*2;
				f.r.right = e.point.x + rad*2;
				f.r.bottom = e.point.y + rad*2;

				if (f.Type == 4 || f.Type == 5) // Peaking/Notch
				{
					// BW = f/Q
					float f1 = f.fr * (sqrt(1.0f + (1.0f / (4.0f * f.Q * f.Q))) - 1.0f / (2.0f * f.Q));
					float f2 = f.fr * (sqrt(1.0f + (1.0f / (4.0f * f.Q * f.Q))) + 1.0f / (2.0f * f.Q));

					float z1 = Freq2X(f1);
					float z2 = Freq2X(f2);

//					r->FillRectangle({z1,e.point.y,z2,ymiddle}, YellowBrush);
					pgs->BeginFigure({z1,ymiddle}, D2D1_FIGURE_BEGIN_HOLLOW);
					D2D1_QUADRATIC_BEZIER_SEGMENT b1;
					b1.point2 = e.point;
					b1.point1.x = e.point.x;
					b1.point1.y = ymiddle;
					pgs->AddQuadraticBezier(b1);
					D2D1_QUADRATIC_BEZIER_SEGMENT b2;
					b2.point2.x = z2;
					b2.point2.y = ymiddle;
					b2.point1.x = e.point.x;
					b2.point1.y = ymiddle;
					pgs->AddQuadraticBezier(b2);

					pgs->EndFigure(D2D1_FIGURE_END_OPEN);
				}

				if (f.sf && f.A)
					CreateTheResponse(f, pgs);

				// Line to fr
				D2D1_POINT_2F l1, l2;
				l1 = e.point;
				l1.y += rad;
				if (f.Type == 0 || f.Type == 1)
					l1.y = rc.top;
				l2.x = e.point.x;
				l2.y = rc.bottom - 35;
				r->DrawLine(l1, l2, f.S ? SelectBrush : WhiteBrush, 1.0f, 0);

				// Fr
				D2D1_RECT_F ly = {};
				ly.left = e.point.x - 20;
				ly.top = rc.bottom - 35;
				ly.right = e.point.x + 20;
				ly.bottom = rc.bottom - 10;
				swprintf_s(t, 1000, L"%s", HZString(f.fr).c_str());
				Text->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
				Text->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_FAR);
				r->DrawTextW(t, (UINT32)wcslen(t), Text, ly, WhiteBrush);
			}

			pgs->Close();
			r->DrawGeometry(pg, YellowBrush,2.5f);
			
			PaintTop(r,rrc);

		}

		float X2Freqr(float x)
		{
			// in width, Max
			// in x , ?
	//		return  Max * x / (rc.right - rc.left);
//			return  0 + exp(log(Max) * x / (rc.right - rc.left));
			return  (float)(20 + exp(log(Max) * x / (rc.right - rc.left)));
		}
		float Freq2X(float f)
		{
			// in width, Max
			// ? ,       f
//			return (rc.right - rc.left) * f / Max;
			return (float)((rc.right - rc.left)*log(f - 20)/log(Max));
		}


		int FilterHitTest(float x,float y)
		{
			for (size_t i = 0 ; i < filters.size() ; i++)
			{
				auto& f = filters[i];
				if (InRect<>(f.r, x, y))
					return (int)i;
			}
			return -1;
		}

		virtual void LeftDown(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);
			for (auto& f : filters)
				f.S = false;

			// Hit a filter?
			int h = FilterHitTest(x, y);
			if (h >= 0)
			{
				filters[h].S = true;
			}
			Redraw();
			return;
		}

		virtual void RightDown(WPARAM ww, LPARAM ll) {
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);
			int h = FilterHitTest(x, y);
			wchar_t re[1000] = { 0 };
			if (h >= 0)
			{

				HMENU hPr = CreatePopupMenu();
				if (filters[h].SpecialFilter == 1)
				{
//					AppendMenu(hPr, MF_STRING, 401, L"Biquad Low cut");
					AppendMenu(hPr, MF_STRING, 402, L"Butterworth Low cut");
					AppendMenu(hPr, MF_STRING, 403, L"Chebyshev I Low cut");
					AppendMenu(hPr, MF_STRING, 404, L"Chebyshev II Low cut");
					AppendMenu(hPr, MF_STRING, 405, L"Elliptic Low cut");

//					AppendMenu(hPr, MF_STRING, 406, L"Biquad Low shelf");
					AppendMenu(hPr, MF_STRING, 407, L"Butterworth Low shelf");
					AppendMenu(hPr, MF_STRING, 408, L"Chebyshev I Low shelf");
					AppendMenu(hPr, MF_STRING, 409, L"Chebyshev II Low shelf");

					CheckMenuItem(hPr, filters[h].SpecialType + 401, MF_CHECKED);
				}
				if (filters[h].SpecialFilter == 2)
				{
//					AppendMenu(hPr, MF_STRING, 401, L"Biquad High cut");
					AppendMenu(hPr, MF_STRING, 402, L"Butterworth High cut");
					AppendMenu(hPr, MF_STRING, 403, L"Chebyshev I High cut");
					AppendMenu(hPr, MF_STRING, 404, L"Chebyshev II High cut");
					AppendMenu(hPr, MF_STRING, 405, L"Elliptic High cut");

//					AppendMenu(hPr, MF_STRING, 406, L"Biquad High shelf");
					AppendMenu(hPr, MF_STRING, 407, L"Butterworth High shelf");
					AppendMenu(hPr, MF_STRING, 408, L"Chebyshev I High shelf");
					AppendMenu(hPr, MF_STRING, 409, L"Chebyshev II High shelf");

					CheckMenuItem(hPr, filters[h].SpecialType + 401, MF_CHECKED);
				}
				if (filters[h].SpecialFilter == 0)
				{
//					AppendMenu(hPr, MF_STRING, 102, L"Low cut");
//					AppendMenu(hPr, MF_STRING, 101, L"High cut");
//					AppendMenu(hPr, MF_STRING, 103, L"Low shelf");
//					AppendMenu(hPr, MF_STRING, 104, L"High shelf");
/*					AppendMenu(hPr, MF_STRING, 105, L"Notch");
					AppendMenu(hPr, MF_STRING, 106, L"Peaking");
					if (filters[h].Type == -1)
						CheckMenuItem(hPr, 100, MF_CHECKED);
					else
						CheckMenuItem(hPr, 101 + filters[h].Type, MF_CHECKED);
*/
				}

				if (filters[h].Type != -1)
				{
					AppendMenu(hPr, MF_SEPARATOR, 0, L"");
					AppendMenu(hPr, MF_STRING, 201, L"Frequency...");
					if (filters[h].SpecialFilter == 0)
					{
						AppendMenu(hPr, MF_STRING, 202, L"dB...");
					}
					else
					{
						if (filters[h].SpecialFilter != 0 && filters[h].SpecialType != 0)
						{
							AppendMenu(hPr, MF_STRING, 251, L"Order...\tUp/Down");
							AppendMenu(hPr, MF_STRING, 252, L"Ripple...\tCtrl+Up/Down");
							if (filters[h].SpecialType == 4)
								AppendMenu(hPr, MF_STRING, 253, L"Roll off...\tAlt+Up/Down");
						}

						// Order, Ripple, ...
					}
					if (filters[h].Type == 4 || filters[h].Type == 5)
						AppendMenu(hPr, MF_STRING, 203, L"Q...");
					AppendMenu(hPr, MF_SEPARATOR, 0, L"");
					AppendMenu(hPr, MF_STRING, 207, L"Active");
					if (filters[h].A)
						CheckMenuItem(hPr, 207, MF_CHECKED);
				}
				if (filters[h].SpecialFilter == 0)
				{
					AppendMenu(hPr, MF_SEPARATOR, 0, L"");
					AppendMenu(hPr, MF_STRING, 299, L"Delete");
				}

				POINT po;
				GetCursorPos(&po);
				int tcmd = TrackPopupMenu(hPr, TPM_CENTERALIGN | TPM_RETURNCMD, po.x, po.y, 0, hParent, 0);
				DestroyMenu(hPr);
				if (tcmd == 0)
					return;

				if (tcmd == 251)
				{
					swprintf_s(re, 1000, L"%i", filters[h].Order);
					if (!AskText(hParent, L"Order", L"Enter filter order:", re))
						return;
					filters[h].Order = (int)_wtoi(re);
					if (filters[h].Order <= 1)
						filters[h].Order = 1;
					Dirty(true);
					Redraw();
				}

				if (tcmd == 252)
				{
					swprintf_s(re, 1000, L"%.1f", filters[h].ripple);
					if (!AskText(hParent, L"Ripple", L"Enter filter ripple (dB):", re))
						return;
					filters[h].ripple = (float)_wtof(re);
					if (filters[h].ripple <= 0.0f)
						filters[h].ripple = 0.5f;
					Dirty(true);
					Redraw();
				}

				if (tcmd == 253)
				{
					swprintf_s(re, 1000, L"%.1f", filters[h].rolloff);
					if (!AskText(hParent, L"Roll off", L"Enter elliptic filter roll off(dB):", re))
						return;
					filters[h].rolloff = (float)_wtof(re);
					if (filters[h].rolloff <= 1)
						filters[h].rolloff = 1;
					Dirty(true);
					Redraw();
				}

				if (tcmd == 207)
				{
					filters[h].A = !filters[h].A;
					Dirty(true);
					Redraw();
					return;
				}
				if (tcmd >= 401 && tcmd <= 410)
				{
					auto& fh = filters[h];
					fh.SpecialType = tcmd - 401;
					if (fh.SpecialType > 4)
					{
						// Shelf
						if (fh.Type == 0 || fh.Type == 3)
							fh.Type = 3;
						if (fh.Type == 1 || fh.Type == 2)
							fh.Type = 2;
					}
					else
					{
						// Shelf
						if (fh.Type == 0 || fh.Type == 3)
							fh.Type = 0;
						if (fh.Type == 1 || fh.Type == 2)
							fh.Type = 1;
						fh.dB = 0;

					}
					Dirty(true);
					Redraw();
					return;
				}

				if (tcmd == 201)
				{
					swprintf_s(re, 1000, L"%.1f", filters[h].fr);
					if (!AskText(hParent, L"Frequency", L"Enter center frequency:", re))
						return;
					filters[h].fr = (float)_wtof(re);
					Dirty(true);
					Redraw();
					return;
				}
				if (tcmd == 202)
				{
					swprintf_s(re, 1000, L"%.1f", filters[h].dB);
					if (!AskText(hParent, L"dB", L"Enter dB:", re))
						return;
					filters[h].dB = (float)_wtof(re);
					Dirty(true);
					Redraw();
					return;
				}
				if (tcmd == 203)
				{
					swprintf_s(re, 1000, L"%.1f", filters[h].Q);
					if (!AskText(hParent, L"Q", L"Enter Q:", re))
						return;
					filters[h].Q = (float)_wtof(re);
					Dirty(true);
					Redraw();
					return;
				}

				if (tcmd == 299)
				{
					filters.erase(filters.begin() + h);
					Dirty(true);
					Redraw();
					return;
				}

				if (tcmd == 101)
				{
					filters[h].Type =  0;
					filters[h].dB = 0;
				}
				if (tcmd == 102)
				{
					filters[h].Type = 1;
					filters[h].dB = 0;
				}
				if (tcmd == 103)
				{
					filters[h].Type = 2;
				}
				if (tcmd == 104)
				{
					filters[h].Type = 3;
				}
				if (tcmd == 105)
				{
					filters[h].Type = 4;
					filters[h].dB = 0;
				}
				if (tcmd == 106)
				{
					filters[h].Type = 5;
				}

				Dirty(true);
				Redraw();
			}
		}
		virtual void LeftUp(WPARAM ww, LPARAM ll)
		{

		}

		virtual void KeyDown(WPARAM ww, LPARAM ll)
		{
			bool Shift = ((GetAsyncKeyState(VK_SHIFT) & 0x8000) != 0);
			bool Control = ((GetAsyncKeyState(VK_CONTROL) & 0x8000) != 0);
			bool Alt = ((GetAsyncKeyState(VK_MENU) & 0x8000) != 0);

			// Change Parameters
			for (auto& f : filters)
			{
				if (f.S && f.sf)
				{
					if (ww == 'A')
					{
						f.A = !f.A;
						f.Build(sr);
						Redraw();
					}

					if (ww == VK_DOWN && !Shift && !Control && !Alt) // Change Order
					{
						if (f.Order > 1)
							f.Order--;
						f.Build(sr);
						Redraw();
					}
					if (ww == VK_UP && !Shift && !Control && !Alt) // Change Order
					{
						f.Order++;
						f.Build(sr);
						Redraw();
					}
					if (ww == VK_DOWN && Control && !Alt) // Change Ripple
					{
						if (Shift)
						{
							if (f.ripple > 0.1)
								f.ripple -= 0.1f;
						}
						else
						{
							if (f.ripple > 1)
								f.ripple -= 1;

						}
						f.Build(sr);
						Redraw();
					}
					if (ww == VK_UP && Control && !Alt) // Change ripple
					{
						if (Shift)
						{
							f.ripple += 0.1f;
						}
						else
						{
							f.ripple += 1;
						}
						f.Build(sr);
						Redraw();
					}
					if (ww == VK_DOWN && !Control && Alt) // Change Rolloff
					{
						if (Shift)
						{
							if (f.rolloff > 0.1)
								f.rolloff -= 0.1f;
						}
						else
						{
							if (f.rolloff > 1)
								f.rolloff -= 1;

						}
						f.Build(sr);
						Redraw();
					}
					if (ww == VK_UP && !Control && Alt) // Change rolloff
					{
						if (Shift)
						{
							if (f.rolloff <= 17.8)
								f.rolloff += 0.1f;
						}
						else
						{
							if (f.rolloff < 17)
								f.rolloff += 1;
						}
						f.Build(sr);
						Redraw();
					}
				}
			}
		}


		virtual void MouseWheel(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			POINT pt;
			GetCursorPos(&pt);
			ScreenToClient(hParent, &pt);
			x = (FLOAT)pt.x;
			y = (FLOAT)pt.y;
			signed short HW = HIWORD(ww);

			// Change Q
			for (auto& f : filters)
			{
				if (f.S && (f.Type == 5 || f.Type == 4))
				{
					if (HW < 0)
						f.Q += 0.5;
					else
						f.Q -= 0.5;
					if (f.Q < 1.0f)
						f.Q = 1.0f;
					if (f.Q > 20)
						f.Q = 20;
					Dirty(true);
					Redraw();
				}
			}

		}


		virtual void MouseMove(WPARAM ww, LPARAM ll)
		{
			bool LeftClick = ((GetAsyncKeyState(VK_LBUTTON) & 0x8000) != 0);

			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);

			for (auto& f : filters)
			{
				if (f.S && LeftClick)
				{
					auto fr2 = X2Freqr(x);
					int Found = 0;
					bool OK = true;
					for (size_t ifi = 0; ifi < filters.size(); ifi++)
					{
						if (&filters[ifi] == &f)
						{
							Found = 1;
							continue;
						}
						if (!Found)
						{
							if (filters[ifi].fr >= fr2)
							{
								OK = false;
								break;
							}
						}
						if (Found)
						{
							if (filters[ifi].fr <= fr2)
							{
								OK = false;
								break;
							}
						}
					}
					if (OK)
					{
						f.fr = fr2;
						if (f.Type > 1 && f.Type != -1 && f.Type != 4)
							f.dB = V2dB(Y2V(y),f.SpecialFilter ? true : false);
						Dirty(true);
						Redraw();
					}
					return;
				}
			}
		}
		virtual void LeftDoubleClick(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);
			int h = FilterHitTest(x, y);
			if (h >= 0)
			{
				if (filters[h].SpecialFilter == 0)
				{
					filters.erase(filters.begin() + h);
					Dirty(true);
					Redraw();
				}
				return;
			}



			// Insert
			PFILT f;
			f.Type = 5;
			f.S = true;
			f.A = true;
			f.fr = X2Freqr(x);
			filters.push_back(f);
			std::sort(filters.begin(), filters.end());
			Dirty(true);
			Redraw();

		}

		virtual void Ser(XML3::XMLElement& e)
		{
			e.RemoveAllElements();
			auto& ee = e["filters"];
			for (auto& p : filters)
				p.Ser(ee.AddElement("f"));
		}
		virtual void Unser(XML3::XMLElement& e)
		{
			auto& ee = e["filters"];
			if (ee.GetChildrenNum() > 0)
			{
				filters.clear();

				for (auto& eee : ee)
				{
					PFILT p;
					p.Unser(eee);
					filters.push_back(p);
				}
			}
		}

		std::vector<sf_sample_st> inb;
		std::vector<sf_sample_st> outb;
		int sr = 0;

		virtual void Build(int SR)
		{
			sr = SR;
			for (auto& b : filters)
				b.Build(SR);
		}

		virtual void Prepare(int SR)
		{
			Build(SR);
		}
		
		// Multi channel
		virtual bool Run2(int SR, int nch, float** in, int ns, float** out)
		{
			// TODO: Stereo Processing
			for (int i = 0; i < nch; i++)
				Run(SR, in[i], ns, out[i]);

/*			if (filters.empty())
				return true;
			for (auto& b : filters)
			{
				if (b.st.b0 == 0 && b.st.b1 == 0 && b.st.b2 == 0 && b.sf == 0)
					b.Build(SR);
			}

			inb.resize(ns);
			outb.resize(ns);
			while(nch > 1)
			{
				for (size_t i = 0; i < ns; i++)
				{
					inb[i].L = in[nch - 2][i];
					inb[i].R = in[nch - 1][i];
				}
				for (auto& b : filters)
				{
					if (!b.A)
						continue;
					if (b.sf)
					{
						if (b.SpecialFilter == 1) // Low cut
						{
						}
					}
					else
					{
						sf_biquad_process(&b.st, (int)ns, inb.data(), outb.data());
					}
				}
				for (size_t i = 0; i < ns; i++)
				{
					out[nch - 2][i] = outb[i].L;
					out[nch - 1][i] = outb[i].R;
				}

				nch -= 2;
			}
			if (nch == 1)
				Run(SR, in[0], ns, out[0]);
*/
			return true;
		}

		// Single channel
		virtual void Run(int SR, float* in, int ons, float* outd)
		{
			if (filters.empty())
				return;
			for (auto& b : filters)
			{
				if (b.st.b0 == 0 && b.st.b1 == 0 && b.st.b2 == 0)
					b.Build(SR);
			}
		
			//auto nns = ons;
			inb.resize(ons);
			outb.resize(ons);

			bool One = false;
			for (auto& b : filters)
			{
				One = true;
				if (b.sf)
				{
					if (b.SpecialFilter == 1) // Low cut - this is the first filter
					{
						// Copy clear to out
						memcpy(outd, in, ons * sizeof(float));

						if (!b.A)
							continue;

						if (b.SpecialType == 1)
						{
							std::shared_ptr<ButtHP> fx;
							fx = std::static_pointer_cast<ButtHP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 2)
						{
							std::shared_ptr<Che1HP> fx;
							fx = std::static_pointer_cast<Che1HP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 3)
						{
							std::shared_ptr<Che2HP> fx;
							fx = std::static_pointer_cast<Che2HP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 4)
						{
							std::shared_ptr<EllHP> fx;
							fx = std::static_pointer_cast<EllHP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 6)
						{
							std::shared_ptr<ButtLPs> fx;
							fx = std::static_pointer_cast<ButtLPs>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 7)
						{
							std::shared_ptr<Che1LPs> fx;
							fx = std::static_pointer_cast<Che1LPs>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 8)
						{
							std::shared_ptr<Che2LPs> fx;
							fx = std::static_pointer_cast<Che2LPs>(b.sf);
							fx->process(ons, &outd);
						}
					}
					if (b.SpecialFilter == 2) // High cut - this is the last filter
					{

						if (!b.A)
							continue;

						if (b.SpecialType == 1)
						{
							std::shared_ptr<ButtLP> fx;
							fx = std::static_pointer_cast<ButtLP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 2)
						{
							std::shared_ptr<Che1LP> fx;
							fx = std::static_pointer_cast<Che1LP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 3)
						{
							std::shared_ptr<Che2LP> fx;
							fx = std::static_pointer_cast<Che2LP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 4)
						{
							std::shared_ptr<EllLP> fx;
							fx = std::static_pointer_cast<EllLP>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 6)
						{
							std::shared_ptr<ButtHPs> fx;
							fx = std::static_pointer_cast<ButtHPs>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 7)
						{
							std::shared_ptr<Che1HPs> fx;
							fx = std::static_pointer_cast<Che1HPs>(b.sf);
							fx->process(ons, &outd);
						}
						if (b.SpecialType == 8)
						{
							std::shared_ptr<Che2HPs> fx;
							fx = std::static_pointer_cast<Che2HPs>(b.sf);
							fx->process(ons, &outd);
						}
					}
				}
				else
				{
					if (!b.A)
						continue;

					for (size_t i = 0; i < ons; i++)
					{
						inb[i].L = outd[i];
						inb[i].R = outd[i];
					}

					sf_biquad_process(&b.st, (int)ons, inb.data(), outb.data());
					for (size_t i = 0; i < ons; i++)
					{
						outd[i] = outb[i].L;
					}
				}
			}
			if (!One)
			{
				memcpy(outd, in, ons * sizeof(float));
			}
			return;
		}
	};

	class GRAPHICEQ : public EQ
	{
	public:
		BAND* ResizingBand = 0;
		vector<BAND> bands;

		virtual void Ser(XML3::XMLElement& e)
		{
			e.RemoveAllElements();
			auto& ee = e["bands"];
			for (auto& p : bands)
				p.Ser(ee.AddElement("b"));
		}
		virtual void Unser(XML3::XMLElement& e)
		{
			auto& ee = e["bands"];
			bands.clear();
			for (auto& eee : ee)
			{
				BAND p;
				p.Unser(eee);
				bands.push_back(p);
			}
		}

		virtual void LeftDoubleClick(WPARAM ww, LPARAM ll)
		{
			if (true)
				return;
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);

			for (size_t i = 0; i < bands.size(); i++)
			{
				auto& b = bands[i];
				if (InRect<>(b.r, x, y))
				{
					// X pos frequency split
					// In full w, Max
					// in x   ?
					float F = Max * x / (rc.right - rc.left);


					float ot = b.to;
					b.to = F;


					BAND b2;
					b2.from = F;
					b2.to = ot;
					b2.V = b.V;
					bands.insert(bands.begin() + i + 1, b2);



					Redraw();
					break;
				}
			}

		}

		virtual void KeyDown(WPARAM ww, LPARAM ll)
		{

		}


		virtual void MouseWheel(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			POINT pt;
			GetCursorPos(&pt);
			ScreenToClient(hParent, &pt);
			x = (FLOAT)pt.x;
			y = (FLOAT)pt.y;
			signed short HW = HIWORD(ww);

			// Change Q
			for (auto& f : bands)
			{
				if (InRect<>(f.r,x,y))
				{
					bool Shift = ((GetAsyncKeyState(VK_SHIFT) & 0x8000) != 0);
					float dB = V2dB(f.V,false);

					if (Shift)
					{
						if (HW < 0)
							dB -= 0.1f;
						else
							dB += 0.1f;
					}
					else
					{
						if (HW < 0)
							dB -= 1.0f;
						else
							dB += 1.0f;
					}

					f.V = dB2V(dB);
					Dirty(true);
					Redraw();
				}
			}

		}

		virtual void MouseMove(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);

			bool Left = ((GetAsyncKeyState(VK_LBUTTON) & 0x8000) != 0);
			if (!Left)
				ResizingBand = 0;


			if (!ResizingBand)
			{
				for (auto& b : bands)
				{
					if (InRect<>(b.r2, x, y))
					{
						SetCursor(ResizeCursor);
						return;
					}
				}
				SetCursor(ArrowCursor);
			}
			else
			{
				SetCursor(ResizeCursor);

				auto& b = *ResizingBand;
				// In full, 2
				// In y, ?
				b.V = (2 * y) / (rc.bottom - rc.top);
				b.V = 2 - b.V;
				
				for (float dbs : { -15.0f, -12.0f, -6.0f, -3.0f, 0.0f, 3.0f, 6.0f, 12.0f, 15.0f })
				{
					float dB = V2dB(b.V,false);
					if (fabs(dB - dbs) < 0.3f)
					{
						b.V = dB2V(dbs);
						break;
					}
				}

				Dirty(true);
				Redraw();
			}
			/*

			bool Left = ((GetAsyncKeyState(VK_LBUTTON) & 0x8000) != 0);
			if (Left)
				LeftDown(ww, ll);
				*/
//			Redraw();
		}

		virtual void RightDown(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);

			HMENU hPr = CreatePopupMenu();
			AppendMenu(hPr, MF_STRING, 101, L"10 bands");
			AppendMenu(hPr, MF_STRING, 102, L"20 bands");
			AppendMenu(hPr, MF_STRING, 103, L"31 bands");
			AppendMenu(hPr, MF_STRING | MF_SEPARATOR, 0, L"");
			//wchar_t t[1000] = { 0 };

			HMENU hPr2 = CreatePopupMenu();
			AppendMenu(hPr, MF_STRING | MF_POPUP, (UINT_PTR)hPr2, L"FFT Size");
			AppendMenu(hPr2, FFTSize == 1024 ? MF_STRING | MF_CHECKED : MF_STRING, 111, L"1024");
			AppendMenu(hPr2, FFTSize == 2048 ? MF_STRING | MF_CHECKED : MF_STRING, 112, L"2048");
			AppendMenu(hPr2, FFTSize == 4096 ? MF_STRING | MF_CHECKED : MF_STRING, 113, L"4096");
			AppendMenu(hPr, MF_STRING, 121, L"Post Gain...");
			AppendMenu(hPr, MF_STRING | MF_SEPARATOR, 0, L"");
			AppendMenu(hPr, MF_STRING, 201, L"Melody boost");
			AppendMenu(hPr, MF_STRING, 202, L"Vocal boost");
			AppendMenu(hPr, MF_STRING, 203, L"Bass lift");
			AppendMenu(hPr, MF_STRING, 204, L"Bass cut");
			AppendMenu(hPr, MF_STRING, 205, L"High lift");
			AppendMenu(hPr, MF_STRING, 206, L"High cut");
			AppendMenu(hPr, MF_STRING, 207, L"Low pass");
			AppendMenu(hPr, MF_STRING, 208, L"High pass");

			POINT po;
			GetCursorPos(&po);
			int tcmd = TrackPopupMenu(hPr, TPM_CENTERALIGN | TPM_RETURNCMD, po.x, po.y, 0, hParent, 0);
			DestroyMenu(hPr);
			if (tcmd == 0)
				return;

			if (tcmd == 111) FFTSize = 1024;
			if (tcmd == 112) FFTSize = 2048;
			if (tcmd == 113) FFTSize = 4096;

			if (tcmd == 101) FixBands(10);
			if (tcmd == 102) FixBands(20);
			if (tcmd == 103) FixBands(31);

			if (tcmd == 208)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i <= 4)
						bands[i].V = dB2V(-48);
				}
			}
			if (tcmd == 207)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i >= 5)
						bands[i].V = dB2V(-48);
				}
			}
			if (tcmd == 206)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i >= 8)
						bands[i].V = dB2V(-5);
					else
						if (i >= 5)
							bands[i].V = dB2V(-3);
				}
			}
			if (tcmd == 205)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i >= 8)
						bands[i].V = dB2V(5);
					else
						if (i >= 5)
							bands[i].V = dB2V(3);
				}
			}
			if (tcmd == 204)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i <= 1)
						bands[i].V = dB2V(-5);
					else
						if (i <= 4)
							bands[i].V = dB2V(-3);
				}
			}
			if (tcmd == 203)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i <= 1)
						bands[i].V = dB2V(5);
					else
					if (i <= 4)
						bands[i].V = dB2V(3);
				}
			}
			if (tcmd == 202)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i == 6) // 2Khhz
						bands[i].V = dB2V(7);
					else
						bands[i].V = dB2V(0);
				}
			}
			if (tcmd == 201)
			{
				FixBands(10);
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (i == 6) // 2Khhz
						bands[i].V = dB2V(5);
					else
						if (i == 5) // 2Khhz
							bands[i].V = dB2V(3);
						else
							if (i == 7) // 2Khhz
								bands[i].V = dB2V(3);
							else
							bands[i].V = dB2V(0);
				}
			}

			Dirty(true);
			Redraw();
		}

		virtual void LeftDown(WPARAM ww, LPARAM ll)
		{
			float x = 0, y = 0;
			x = (FLOAT)LOWORD(ll);
			y = (FLOAT)HIWORD(ll);

			for (auto& b : bands)
			{
				if (InRect<>(b.r2, x, y))
				{
					ResizingBand = &b;
					Redraw();
					break;
				}
			}
		}
		virtual void LeftUp(WPARAM ww, LPARAM ll)
		{
			ResizingBand = 0;
		}



		void FixBands(int n) // 10,20,31
		{
			XMode = n;
			if (n == 10)
			{
				bands.clear();
				bands.resize(10);
				float st = 32.0f;
				bands[0].from = 0;
				for (int i = 0; i < 10; i++)
				{
					bands[i].to = st;
					st *= 2;
					if (i > 0)
						bands[i].from = bands[i - 1].to;
				}
			}
			if (n == 20)
			{
				bands.clear();
				bands.resize(n);
				float st = 10.0f;
				bands[0].from = 0;
				for (int i = 0; i < n; i++)
				{
					bands[i].to = st;
					st *= 3.0f;
					st /= 2.0f;
					if (i > 0)
						bands[i].from = bands[i - 1].to;
				}
			}
			if (n == 31)
			{
				bands.clear();
				bands.resize(n);
//				float st = 20.0f;
				bands[0].from = 0;

				float b31[] = { 20,25,31.5,40,50 , 63 , 80 , 100 , 125 , 160 , 200 , 250 , 315 , 400 , 500 , 630 , 800 , 1000 , 1250 , 1600 , 2000 , 2500 , 3150 , 4000 , 5000 , 6300 , 8000 , 10000 , 12500 , 16000 , 20000 };

				for (int i = 0; i < n; i++)
				{
					bands[i].to = b31[i];
//					st += st/3.0f;
					if (i > 0)
						bands[i].from = bands[i - 1].to;
				}
			}

			if (!bands.empty())
			bands[bands.size() - 1].to = (float)Max*2;
		}

		float Freq2X(float fr)
		{
			float width = rc.right - rc.left;
			if (XMode > 0 && !bands.empty())
			{
				size_t iSection = 0;
				for (size_t i = 0; i < bands.size(); i++)
				{
					if (bands[i].to < fr)
						continue;
					iSection = i;
					break;
				}

				// In bsize , width
				// in iSection, ? 

				float w0 = ((iSection + 0) * width) / bands.size();
				float w1 = ((iSection + 1) * width) / bands.size();
				float wi = w1 - w0;
				float fri = bands[iSection].to - bands[iSection].from;

				// linear from w0 to w1
				// in wi, fri
				// ?    , fr
				fr -= bands[iSection].from;
				float Y = (wi * fr) / fri;
				Y += w0;
				return Y;
			}



			// else Linear
			return width * fr / Max;
		}


		void Paint(ID2D1Factory* fact, ID2D1RenderTarget* r, RECT rrc)
		{
			wchar_t t[1000] = { 0 };
			CreateBrushes(r);
			rc = FromR(rrc);
			auto rr = rc;
			r->FillRectangle(rc, BGBrush);
			if (bands.empty())
			{
				FixBands(10);
			}

			PaintDBLines(r);


			for (size_t i = 0; i < bands.size(); i++)
			{

				// in max, width
				// in from  ?
				D2D1_POINT_2F p1, p2;
				p1.x = Freq2X(bands[i].from);
				p2.x = Freq2X(bands[i].to);

				bands[i].r.left = p1.x;
				bands[i].r.right = p2.x;
				bands[i].r2.left = p1.x;
				bands[i].r2.right = p2.x;

				// in 2.0, Max height
				// in V    ?
				float H = bands[i].V * (rc.bottom - rc.top) / 2.0f;
				H = (rc.bottom - rc.top) - H;
				bands[i].r.top = H;
				bands[i].r.bottom = rc.bottom;

				float j = 2;
				p1.y = rc.top + H;
				p2.y = rc.top + H;
				p1.y -= j;
				p2.y -= j;
				D2D1_RECT_F re = { 0 };
				re.left = p1.x;
				re.top = p2.y;
				re.right = p2.x;
//				r->DrawLine(p1, p2, WhiteBrush);
				p1.y += j * 2;
				p2.y += j * 2;
//				r->DrawLine(p1, p2, WhiteBrush);
				re.bottom = p2.y;

				r->FillRectangle(bands[i].r, GrayBrush);


				r->FillRectangle(re, WhiteBrush);



				p1.y -= j;
				p2.y -= j;
				bands[i].r2.top = p1.y - j * 2;
				bands[i].r2.bottom = p1.y + j * 2;


				float rad = 12.0f;
				// must take 10% of width
				rad = bands[i].r2.right - bands[i].r2.left;
				rad /= 10.0f;



				if (true)
				{
					D2D1_ELLIPSE el;
					el.point = p2;
					el.radiusX = el.radiusY = rad;
					r->FillEllipse(el, WhiteBrush);
				}

				D2D1_POINT_2F p3, p4;
				if (true)
				{
					p3.x = p1.x;
					p4.x = p1.x;
					p3.y = p1.y;
					p4.y = rc.bottom;
					r->DrawLine(p3, p4, WhiteBrush);

					D2D1_ELLIPSE el;
					el.point = p3;
					el.radiusX = el.radiusY = rad;
					r->FillEllipse(el, WhiteBrush);
				}

				D2D1_RECT_F ly;
				if (true)
				{
					p3.x = p2.x;
					p4.x = p2.x;
					p3.y = p2.y;
					p4.y = rc.bottom;
					r->DrawLine(p3, p4, WhiteBrush);

					D2D1_ELLIPSE el;
					el.point = p3;
					el.radiusX = el.radiusY = rad;
					r->FillEllipse(el, WhiteBrush);


					Text->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
					Text->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_FAR);
					ly.left = p1.x;
					ly.top = bands[i].r2.top;
					ly.right = p2.x;
					ly.bottom = rc.bottom;
				}

				ly.bottom -= 20;
//				swprintf_s(t, 1000, L"%s", HZString(bands[i].from, bands[i].to).c_str());
				swprintf_s(t, 1000, L"%s", HZString(bands[i].to).c_str());
				r->DrawTextW(t, (UINT32)wcslen(t), Text, ly, WhiteBrush);
				Text->SetParagraphAlignment(DWRITE_PARAGRAPH_ALIGNMENT_NEAR);
				ly.top += 20;
				swprintf_s(t, 1000, L"%.0f dB", V2dB(bands[i].V,false));
				r->DrawTextW(t, (UINT32)wcslen(t), Text, ly, WhiteBrush);
			}


			if (true)
			{
				POINT p;
				GetCursorPos(&p);
				ScreenToClient(hParent, &p);
				D2D1_POINT_2F p1, p2;

				p1.y = (FLOAT)p.y;
				p1.x = rc.left;
				p2.x = rc.right;
				p2.y = (FLOAT)p.y;
//				r->DrawLine(p1, p2, SelectBrush);

				// In height , 2.0f
				// in y , ?
/*				float V = Y2V(p.y);
				swprintf_s(t,1000,L"%.0f dB",V2dB(V));

				D2D1_RECT_F ly;
				ly.left = p.x + 35;
				ly.top = p.y - 35;
				ly.right = ly.left + 100;
				ly.bottom = ly.top + 100;
				Text->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
				r->DrawTextW(t, wcslen(t), Text, ly,WhiteBrush);
				ly.top = p.y + 35;
				ly.bottom = ly.top + 100;
				swprintf_s(t, 1000, L"%s", HZString(X2Freq(p.x)).c_str());
				r->DrawTextW(t, wcslen(t), Text, ly, WhiteBrush);
				*/
				p1.x = (FLOAT)p.x;
				p1.y = rc.top;
				p2.y = rc.bottom;
				p2.x = (FLOAT)p1.x;
//				r->DrawLine(p1, p2, SelectBrush);
			}

		}
		int FFTSize = 4096;



		virtual void Build(int SR)
		{
			// BW = f0/Q; - f0 center frequency
			for (auto& b : bands)
			{
				float BW = b.to - b.from;
				float f0 = b.from + BW / 2.0f;
				float Q = f0 / BW;
				sf_biquad_state_st& s1 = b.state;
				sf_peaking(&s1, SR, f0, Q, V2dB(b.V,false));
			}
			
		}


		std::vector<sf_sample_st> inb;
		std::vector<sf_sample_st> outb;

		bool UseBiquad = true;

		virtual void Prepare(int SR)
		{
			if (UseBiquad)
				Build(SR);
		}
		virtual bool Run2(int SR, int nch, float** in, int ns, float** out)
		{
			// TODO: Stereo Processing
			for (int i = 0; i < nch; i++)
				Run(SR, in[i], ns, out[i]);
			return true;
		}

		virtual void Run(int SR, float* in, int ons, float* outd)
		{
			if (bands.empty())
				return;
			if (UseBiquad)
			{
//				auto nns = ons;
				if (bands[0].state.b0 == 0 && bands[0].state.b1 == 0 && bands[0].state.b2 == 0)
					Build(SR);
				inb.resize(ons);
				outb.resize(ons);

				for (size_t i = 0; i < ons; i++)
				{
					inb[i].L = in[i];
					inb[i].R = in[i];
				}

				outb = inb;
				for (auto& b : bands)
				{
					if (b.V != 1.0f)
					{
						sf_biquad_process(&b.state, (int)ons, inb.data(), outb.data());
						inb = outb;
					}
				}

				for (size_t i = 0; i < ons; i++)
				{
					outd[i] = outb[i].L;
				}

				if (true)
					return;
				return;
			}

// FFT based, obsolete
/*			float rep[4096] = { 0 };
//			float mags[4096] = { 0 };

			// Make sure we have power of 2
			int jo = 0;
			for (;;)
			{
				if (ons <= 2)
					break;
				int ns = min(ons, FFTSize);
				while (ns > 0 && (ns & (ns - 1)) != 0)
					ns--;
				if (ns <= 2)
					break;

				ffft::FFTReal <float> fft(ns);

				// Transform
				fft.do_fft(rep,in + jo);

				auto GetF = [&](int sr, int i)
				{
					float x = (float)(sr * i);
					return (float)(x / (float)FFTSize);
				};



				for (int j = 0; j < ns / 2; j++)
				{
					float fr = GetF(SR, j);

					bool U = false;
					// Find band
					for (auto& b : bands)
					{
						if (b.from <= fr && b.to >= fr && b.V != 1.0f)
						{
							float aa = rep[j];
							float bb = rep[ns/2 + j];
							float mg = sqrt(aa * aa + bb * bb);
							mg /= FFTSize;
							float phase = atan2(bb, aa);

							float MAG_dB = 20 * log10(mg);

							float dB = V2dB(b.V);
							MAG_dB += dB;

							// log10mg = (db + 20log10mx)/20
							mg = MAG_dB / 20.0f;
							mg = pow(10.0f, mg);

							mg *= FFTSize;

							rep[j] = ((float)(mg * cos(phase)));
							rep[ns/2 + j] = ((float)(mg * sin(phase)));

							aa = rep[j];
							bb = rep[ns / 2 + j];
							//float mgn = sqrt(aa * aa + bb * bb);

							U = true;
							break;
						}
					}
				}

				fft.do_ifft(rep, outd + jo);
				fft.rescale(outd + jo);

				jo += ns;
				ons -= ns;
			}


		// Normalize to 
	//	AUDIOFX::Normalize<float>(1.0f,outd, nns);
*/


		}
	


	};





}