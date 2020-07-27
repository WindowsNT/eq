class WHDR
{
public:

	WAVEHDR wH;
	vector<char> d;

	WHDR(DWORD b = 0, const char* data = 0);
	void r(DWORD b, const char* data = 0);
	operator WAVEHDR*() { return &wH; }
};

class WIN
{
private:

	MMRESULT mm = 0;
	HANDLE hEv = 0;
	HWAVEIN wIN = 0;
	vector<WHDR> buffers;
	bool Started = 0;
	int nbuf = 10;
	function<void(const char*d, DWORD sz)> ff;


public:

	WIN(function<void(const char*d, DWORD sz)> ffx, int bufnum = 10);
	bool Open(int wid = WAVE_MAPPER, DWORD SR = 22050, WORD BPS = 16,int CH = 1);
	void thrx();
	bool Start();
	bool Stop();
	~WIN();
};



class WOUT
{
private:

	MMRESULT mm = 0;
	HWAVEOUT wOUT = 0;
	HANDLE hEv = 0;
	HANDLE hDied = 0;
	HANDLE hCanSend = 0;
	bool Dying = false;
	volatile LONGLONG Count = 0;
	tlock<list<WHDR*>> Pending;
	unsigned int Lat = 0;
	bool First = false;

public:

	WOUT(int lt = 5);
	bool Open(int wid = WAVE_MAPPER, DWORD SR = 22050, WORD BPS = 16,int CH = 1);
	void thrx();
	bool Write(const char* d, size_t sz);
	~WOUT();
};


class CONV
{

	vector<float> SSEHelpBufferA1;
	float* SSEHelpBufferPointerAligned16A1;
	vector<float> SSEHelpBufferA2;
	float* SSEHelpBufferPointerAligned16A2;
	vector<float> SSEHelpBufferB1;
	float* SSEHelpBufferPointerAligned16B1;
	vector<float> SSEHelpBufferB2;
	float* SSEHelpBufferPointerAligned16B2;

public:

	void InitSSE();
	void PrepareSSEHelpBuffer(long frames);
	bool SSEFrom32To16(float *source, void* target, long frames);
	bool SSEFrom16To32(const short*source, float* target, long frames);
	float Saturate(float input, float);
	void f16t32(const short*src, long frm, float*target);
	void f32t16(float *src, long frm, short *dst);
	void f8t32(const signed char*src, long frames, float *target);
	void f32t8(float *source, long frames, signed char *dst);

};


WHDR::WHDR(DWORD b, const char* data)
{
	r(b, data);
}

void WHDR::r(DWORD b, const char* data)
{
	memset(&wH, 0, sizeof(wH));
	d.resize(b);
	if (data)
		memcpy(d.data(), data, b);

	wH.dwBufferLength = b;
	wH.lpData = d.data();
}






WIN::WIN(function<void(const char*d, DWORD sz)> ffx, int bufnum)
{
	nbuf = bufnum;
	ff = ffx;
	hEv = CreateEvent(0, 0, 0, 0);
}

bool WIN::Open(int wid, DWORD SR, WORD BPS,int CH)
{
	WAVEFORMATEX wEx = { 0 };
	wEx.wFormatTag = WAVE_FORMAT_PCM;
	wEx.nChannels = (WORD)CH;
	wEx.nSamplesPerSec = SR;
	wEx.wBitsPerSample = BPS;
	wEx.nBlockAlign = (wEx.wBitsPerSample / 8) * wEx.nChannels;
	wEx.nAvgBytesPerSec = wEx.nBlockAlign * wEx.nSamplesPerSec;
	mm = waveInOpen(&wIN, wid, &wEx, (DWORD_PTR)hEv, 0, CALLBACK_EVENT);
	if (mm)
		return false;

	buffers.resize(nbuf);
	for (auto& b : buffers)
	{
		b.r(wEx.nAvgBytesPerSec / nbuf);
		mm = waveInPrepareHeader(wIN, b, sizeof(WAVEHDR));
		mm = waveInAddBuffer(wIN, b, sizeof(WAVEHDR));
	}


	return true;
}

void WIN::thrx()
{
	__try
	{
		for (;;)
		{
			WaitForSingleObject(hEv, INFINITE);
			if (!Started)
				break;
			if (Ending)
				break;

			for (auto& b : buffers)
			{
				if (b.wH.dwFlags & WHDR_DONE)
				{
					ff(b.wH.lpData, b.wH.dwBytesRecorded);
					if (Ending)
						break;
					mm = waveInAddBuffer(wIN, b, sizeof(WAVEHDR));
					break;
				}
			}

		}
	}
	__except (true)
	{

	}
}

bool WIN::Start()
{
	if (Started)
		return true;
	ResetEvent(hEv);
	Started = true;
	thread t(&WIN::thrx, this);
	t.detach();
	return (waveInStart(wIN) == 0);
}

bool WIN::Stop()
{
	if (!Started)
		return true;
	waveInStop(wIN);
	waveInReset(wIN);
	Started = false;
	SetEvent(hEv);
	return true;
}

WIN::~WIN()
{
	if (wIN)
	{
		waveInStop(wIN);
		waveInReset(wIN);
		waveInClose(wIN);
	}
	wIN = 0;
	if (hEv)
		CloseHandle(hEv);
	hEv = 0;
}




WOUT::WOUT(int lt)
{
	Lat = lt;
	hEv = CreateEvent(0, 0, 0, 0);
	hDied = CreateEvent(0, 0, 0, 0);
	hCanSend = CreateEvent(0, TRUE, TRUE, 0);
}

bool WOUT::Open(int wid, DWORD SR, WORD BPS,int CH)
{
	WAVEFORMATEX wEx = { 0 };
	wEx.wFormatTag = WAVE_FORMAT_PCM;
	wEx.nChannels = (WORD)CH;
	wEx.nSamplesPerSec = SR;
	wEx.wBitsPerSample = BPS;
	wEx.nBlockAlign = (wEx.wBitsPerSample / 8) * wEx.nChannels;
	wEx.nAvgBytesPerSec = wEx.nBlockAlign * wEx.nSamplesPerSec;
	mm = waveOutOpen(&wOUT, wid, &wEx, (DWORD_PTR)hEv, 0, CALLBACK_EVENT);
	if (mm)
		return false;

	thread t(&WOUT::thrx, this);
	t.detach();
	return true;
}

void WOUT::thrx()
{
	for (;;)
	{
		WaitForSingleObject(hEv, INFINITE);
		if (!wOUT || Dying)
		{
			SetEvent(hDied);
			break;
		}

		Pending.writelock([&](list<WHDR*>& p) {

			auto it = p.begin();
			while (it != p.end())
			{
				if ((*it)->wH.dwFlags & WHDR_DONE)
				{
					mm = waveOutUnprepareHeader(wOUT, (*it)->operator WAVEHDR *(), sizeof(WAVEHDR));
					p.erase(it++);
				}
				else
					it++;
			}

			if (p.size() < Lat)
				SetEvent(hCanSend);
			});
	}
}

bool WOUT::Write(const char* d, size_t sz)
{

	auto e = WaitForSingleObject(hCanSend, 10000);
	Pending.writelock([&](list<WHDR*>& p) {

		WHDR* w = new WHDR((DWORD)sz, d);
		p.push_back(w);
		if (p.size() < Lat)
			return;

		if (!First)
		{
			for (auto& a : p)
			{
				mm = waveOutPrepareHeader(wOUT, a->operator WAVEHDR *(), sizeof(WAVEHDR));
				mm = waveOutWrite(wOUT, a->operator WAVEHDR *(), sizeof(WAVEHDR));
			}
			First = true;
		}
		else
		{
			if (e == WAIT_TIMEOUT)
			{
				MessageBox(0, 0, 0, 0);
			}
			mm = waveOutPrepareHeader(wOUT, w->operator WAVEHDR *(), sizeof(WAVEHDR));
			mm = waveOutWrite(wOUT, w->operator WAVEHDR *(), sizeof(WAVEHDR));
			if (p.size() >= Lat)
				ResetEvent(hCanSend);
		}
		});

	return true;

}

WOUT::~WOUT()
{
	Dying = true;
	if (wOUT)
	{
		waveOutReset(wOUT);
		waveOutClose(wOUT);
	}
	wOUT = 0;
	if (hEv)
	{
		SetEvent(hEv);
		if (hDied)
		{
			WaitForSingleObject(hDied, 1000);
			CloseHandle(hDied);
		}
		CloseHandle(hEv);
	}
	hEv = 0;
	if (hCanSend)
		CloseHandle(hCanSend);
	hCanSend = 0;
}

void CONV::InitSSE()
{
	SSEHelpBufferA1.resize(1);
	SSEHelpBufferPointerAligned16A1 = 0;
	SSEHelpBufferA2.resize(1);
	SSEHelpBufferPointerAligned16A2 = 0;
	SSEHelpBufferB1.resize(1);
	SSEHelpBufferPointerAligned16B1 = 0;
	SSEHelpBufferB2.resize(1);
	SSEHelpBufferPointerAligned16B2 = 0;
}

void CONV::PrepareSSEHelpBuffer(long frames)
{
	if (SSEHelpBufferA1.size() <= (unsigned long)frames)
	{
		// Create the buffer and align it
		SSEHelpBufferA1.resize(frames * 2 + 100);
		SSEHelpBufferPointerAligned16A1 = SSEHelpBufferA1.data();
		unsigned long long a = (unsigned long long)SSEHelpBufferPointerAligned16A1;
		while (a % 16)
			a++;
		SSEHelpBufferPointerAligned16A1 = (float*)a;
	}
	if (SSEHelpBufferA2.size() <= (unsigned long)frames)
	{
		// Create the buffer and align it
		SSEHelpBufferA2.resize(frames * 2 + 100);
		SSEHelpBufferPointerAligned16A2 = SSEHelpBufferA2.data();
		unsigned long long a = (unsigned long long)SSEHelpBufferPointerAligned16A2;
		while (a % 16)
			a++;
		SSEHelpBufferPointerAligned16A2 = (float*)a;
	}
	if (SSEHelpBufferB1.size() <= (unsigned long)frames)
	{
		// Create the buffer and align it
		SSEHelpBufferB1.resize(frames * 2 + 100);
		SSEHelpBufferPointerAligned16B1 = SSEHelpBufferB1.data();
		unsigned long long a = (unsigned long long)SSEHelpBufferPointerAligned16B1;
		while (a % 16)
			a++;
		SSEHelpBufferPointerAligned16B1 = (float*)a;
	}
	if (SSEHelpBufferB2.size() <= (unsigned long)frames)
	{
		// Create the buffer and align it
		SSEHelpBufferB2.resize(frames * 2 + 100);
		SSEHelpBufferPointerAligned16B2 = SSEHelpBufferB2.data();
		unsigned long long a = (unsigned long long)SSEHelpBufferPointerAligned16B2;
		while (a % 16)
			a++;
		SSEHelpBufferPointerAligned16B2 = (float*)a;
	}
}

bool CONV::SSEFrom32To16(float *source, void* target, long frames)
{
#ifdef _WIN64
	return false;
#else
	PrepareSSEHelpBuffer(frames);
	unsigned long long jX1 = (unsigned long long)SSEHelpBufferPointerAligned16A1;
	if (jX1 % 16 || jX1 == 0)
		return false;
	unsigned long long jX2 = (unsigned long long)SSEHelpBufferPointerAligned16A2;
	if (jX2 % 16 || jX2 == 0)
		return false;
	memcpy(SSEHelpBufferPointerAligned16A1, source, frames * sizeof(float));

	static const float L = 32767.0f;//(1 << 15);
	static const __m128 LI = { (float)L, (float)L, (float)L, (float)L };

	__m128* floats = (__m128*)SSEHelpBufferPointerAligned16A1;
	__m64* ints = (__m64*)SSEHelpBufferPointerAligned16A2;

	for (int i = 0; i < frames / 4; i++)
		ints[i] = _mm_cvtps_pi16(_mm_mul_ps(floats[i], LI));
	memcpy(target, SSEHelpBufferPointerAligned16A2, frames * 2);
	_mm_empty();
	return true;
#endif
}

bool CONV::SSEFrom16To32(const short*source, float* target, long frames)
{
#ifdef _WIN64
	return false;
#else
	PrepareSSEHelpBuffer(frames);
	unsigned long long jX1 = (unsigned long long)SSEHelpBufferPointerAligned16B1;
	if (jX1 % 16 || jX1 == 0)
		return false;
	unsigned long long jX2 = (unsigned long long)SSEHelpBufferPointerAligned16B2;
	if (jX2 % 16 || jX2 == 0)
		return false;
	memcpy(SSEHelpBufferPointerAligned16B1, source, frames * 2);

	static const float L = 32767.0f;//(1 << 15);
	static const __m128 LI = { (float)L, (float)L, (float)L, (float)L };

	__m64* ints = (__m64*)SSEHelpBufferPointerAligned16B1;
	__m128* floats = (__m128*)SSEHelpBufferPointerAligned16B2;

	for (int i = 0; i < frames / 4; i++)
		floats[i] = _mm_div_ps(_mm_cvtpi16_ps(ints[i]), LI);
	_mm_empty();
	memcpy(target, SSEHelpBufferPointerAligned16B2, frames * sizeof(float));
	return true;
#endif
}

float CONV::Saturate(float input, float)
{
	return input;
	/*
	auto fMax = 1.0f;
	static const float fGrdDiv = 0.5f;
	float x1 = fabs(input + fMax);
	float x2 = fabs(input - fMax);
	return fGrdDiv * (x1 - x2);
*/

}

void CONV::f16t32(const short*src, long frm, float*target)
{
#ifdef _WIN64
	while (--frm >= 0)
		*target++ = ((*src++) + .5f) * (1.0f / 32767.503f);
#else
	SSEFrom16To32(src, target, frm);
#endif
}

void CONV::f32t16(float *src, long frm, short *dst)
{
#ifdef _WIN64
	float finter;
	while (--frm >= 0)
	{
		finter = Saturate(*src++, 1.f);
		*dst++ = (short)((finter * 32767.505f) - .5f);
	}
#else
	SSEFrom32To16(src, dst, frm);
#endif


}

void CONV::f8t32(const signed char*src, long frames, float *target)
{
	while (--frames >= 0)
		*target++ = ((*src++) + .5f) * (1.0f / 256.0f);
}

void CONV::f32t8(float *source, long frames, signed char *dst)
{
	float finter = 0;
	while (--frames >= 0)
	{
		finter = Saturate(*source++, 1.f);
		*dst++ = (signed char)(((finter * 128.0f)) + 0);
	}
}
