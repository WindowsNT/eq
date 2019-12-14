# An equalizer design for Windows with Direct2D

CodeProject article: 


Uses:

My XML for serialization: https://github.com/WindowsNT/xml
My syncs for playback: https://github.com/WindowsNT/mt
DSPFilters for Low/High cut of the parametric equalizer: https://github.com/vinniefalco/DSPFilters
SNDFilter for biquad peaking filters in parametric/graphic equalizer: https://github.com/voidqk/sndfilter



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
![TT](https://www.codeproject.com/KB/miscctrl/1279856/1.jpg)

![Test](https://www.codeproject.com/KB/miscctrl/1279856/2.jpg)


