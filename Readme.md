# An equalizer design for Windows with Direct2D

CodeProject article: https://www.codeproject.com/Articles/5253995/The-Equalizer-Design-Graphic-and-Parametric-Win32


![PEQ](https://www.codeproject.com/KB/miscctrl/5253995/eq.jpg)

![GEQ](https://www.codeproject.com/KB/miscctrl/5253995/geq.jpg)


Uses:

* My XML for serialization: https://github.com/WindowsNT/xml
* My syncs for playback: https://github.com/WindowsNT/mt
* DSPFilters for Low/High cut of the parametric equalizer: https://github.com/vinniefalco/DSPFilters
* SNDFilter for biquad peaking filters in parametric/graphic equalizer: https://github.com/voidqk/sndfilter

Includes a VS 2019 solution to test with a simple MP3 player. Run the app, press Space and load an MP3.
Use EQ::PARAMETRICEQ for parametric equalizer and EQ::GRAPHICEQ for graphic.

Features:

* Graphic equalizer for 10,20,31 bands with some presets, using Peaking biquad filters.
* Parametric equalizer with a low/high cut or low/high shelf filter (Choose from Butterworth, Chebyshev I,II, Elliptic) and unlimited peaking filters with Q.
* Draw the response with Direct2D
* Keyboard shortcuts for activation, order, ripple 
* Mouse movement for dB, frequency, and wheel for Q


```C++
class EQ
{
	virtual void Paint(ID2D1Factory*fact, ID2D1RenderTarget* r, RECT rc) = 0;
	virtual void LeftDown(WPARAM ww, LPARAM ll) = 0;
	virtual void RightDown(WPARAM ww, LPARAM ll) = 0;
	virtual void LeftUp(WPARAM ww, LPARAM ll) = 0;
	virtual void MouseMove(WPARAM ww, LPARAM ll) = 0;
	virtual void MouseWheel(WPARAM ww, LPARAM ll) = 0;
	virtual void KeyDown(WPARAM ww, LPARAM ll) = 0;
	virtual void LeftDoubleClick(WPARAM ww,LPARAM ll) = 0;
	virtual void Ser(XML3::XMLElement& e) = 0;
	virtual void Unser(XML3::XMLElement& e) = 0;

	virtual void Prepare(int SR) = 0;
	virtual void Run(int SR, float* in, int ns, float* out) = 0; // Single channel
	virtual bool Run2(int SR, int nch,float** in, int ns, float** out) = 0; // Multiple channel (ToDo)
	virtual void Build(int SR) = 0; // Create the filters
};

// Callbacks
class EQCALLBACK
{
public:

	virtual void RedrawRequest(EQ* pr) = 0; // Called when the equalizer needs to redraw
	virtual void Dirty(EQ* e,bool) = 0; // Called when the equalizer has changed
};


```

