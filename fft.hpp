// FFT
/*****************************************************************************
*                                                                            *
*       DIGITAL SIGNAL PROCESSING TOOLS                                      *
*       Version 1.03, 2001/06/15                                             *
*       (c) 1999 - Laurent de Soras                                          *
*                                                                            *
*       FFTReal.h                                                            *
*       Fourier transformation of real number arrays.                        *
*       Portable ISO C++                                                     *
*                                                                            *
* Tab = 3                                                                    *
*****************************************************************************/
#ifndef _FFT_HPP
#define _FFT_HPP


#if defined (FFTReal_CURRENT_HEADER)
#error Recursive inclusion of FFTReal header file.
#endif
#define	FFTReal_CURRENT_HEADER

#if ! defined (FFTReal_HEADER_INCLUDED)
#define	FFTReal_HEADER_INCLUDED



#if defined (_MSC_VER)
#pragma pack (push, 8)
#endif	// _MSC_VER



class	FFTReal
{

/*\\\ PUBLIC \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

public:

	/* Change this typedef to use a different floating point type in your FFTs
		(i.e. float, double or long double). */
	typedef float	flt_t;

	explicit			FFTReal(const long length)
			: _length(length)
			, _nbr_bits(int(floor(log((float)length) / log(2.0f) + 0.5f)))
			, _bit_rev_lut(int(floor(log((float)length) / log(2.0f) + 0.5f)))
			, _trigo_lut(int(floor(log((float)length) / log(2.0f) + 0.5f)))
			, _sqrt2_2(flt_t(sqrt(2.0f) * 0.5f))
		{
		//	ASSERT((1L << _nbr_bits) == length);

			_buffer_ptr = 0;
			if (_nbr_bits > 2)
			{
				_buffer_ptr = new flt_t[_length];
			}
		}

	~FFTReal(void)
	{
		delete[] _buffer_ptr;
		_buffer_ptr = 0;
	}



	void do_fft(flt_t f[], const flt_t x[]) const
	{

	/*______________________________________________
	 *
	 * General case
	 *______________________________________________
	 */

		if (_nbr_bits > 2)
		{
			flt_t* sf;
			flt_t* df;

			if (_nbr_bits & 1)
			{
				df = _buffer_ptr;
				sf = f;
			}
			else
			{
				df = f;
				sf = _buffer_ptr;
			}

			/* Do the transformation in several pass */
			{
				int		pass;
				long		nbr_coef;
				long		h_nbr_coef;
				long		d_nbr_coef;
				long		coef_index;

				/* First and second pass at once */
				{
					const long* const	bit_rev_lut_ptr = _bit_rev_lut.get_ptr();
					coef_index = 0;
					do
					{
						const long		rev_index_0 = bit_rev_lut_ptr[coef_index];
						const long		rev_index_1 = bit_rev_lut_ptr[coef_index + 1];
						const long		rev_index_2 = bit_rev_lut_ptr[coef_index + 2];
						const long		rev_index_3 = bit_rev_lut_ptr[coef_index + 3];

						flt_t* const	df2 = df + coef_index;
						df2[1] = x[rev_index_0] - x[rev_index_1];
						df2[3] = x[rev_index_2] - x[rev_index_3];

						const flt_t		sf_0 = x[rev_index_0] + x[rev_index_1];
						const flt_t		sf_2 = x[rev_index_2] + x[rev_index_3];

						df2[0] = sf_0 + sf_2;
						df2[2] = sf_0 - sf_2;

						coef_index += 4;
					} while (coef_index < _length);
				}

				/* Third pass */
				{
					coef_index = 0;
					const flt_t		sqrt2_2 = _sqrt2_2;
					do
					{
						flt_t				v;

						sf[coef_index] = df[coef_index] + df[coef_index + 4];
						sf[coef_index + 4] = df[coef_index] - df[coef_index + 4];
						sf[coef_index + 2] = df[coef_index + 2];
						sf[coef_index + 6] = df[coef_index + 6];

						v = (df[coef_index + 5] - df[coef_index + 7]) * sqrt2_2;
						sf[coef_index + 1] = df[coef_index + 1] + v;
						sf[coef_index + 3] = df[coef_index + 1] - v;

						v = (df[coef_index + 5] + df[coef_index + 7]) * sqrt2_2;
						sf[coef_index + 5] = v + df[coef_index + 3];
						sf[coef_index + 7] = v - df[coef_index + 3];

						coef_index += 8;
					} while (coef_index < _length);
				}

				/* Next pass */
				for (pass = 3; pass < _nbr_bits; ++pass)
				{
					coef_index = 0;
					nbr_coef = 1 << pass;
					h_nbr_coef = nbr_coef >> 1;
					d_nbr_coef = nbr_coef << 1;
					const flt_t* const	cos_ptr = _trigo_lut.get_ptr(pass);
					do
					{
						long				i;
						const flt_t* const sf1r = sf + coef_index;
						const flt_t* const sf2r = sf1r + nbr_coef;
						flt_t* const dfr = df + coef_index;
						flt_t* const dfi = dfr + nbr_coef;

						/* Extreme coefficients are always real */
						dfr[0] = sf1r[0] + sf2r[0];
						dfi[0] = sf1r[0] - sf2r[0];	// dfr [nbr_coef] =
						dfr[h_nbr_coef] = sf1r[h_nbr_coef];
						dfi[h_nbr_coef] = sf2r[h_nbr_coef];

						/* Others are conjugate complex numbers */
						const flt_t* const	sf1i = sf1r + h_nbr_coef;
						const flt_t* const	sf2i = sf1i + nbr_coef;
						for (i = 1; i < h_nbr_coef; ++i)
						{
							const flt_t		c = cos_ptr[i];					// cos (i*PI/nbr_coef);
							const flt_t		s = cos_ptr[h_nbr_coef - i];	// sin (i*PI/nbr_coef);
							flt_t				v;

							v = sf2r[i] * c - sf2i[i] * s;
							dfr[i] = sf1r[i] + v;
							dfi[-i] = sf1r[i] - v;	// dfr [nbr_coef - i] =

							v = sf2r[i] * s + sf2i[i] * c;
							dfi[i] = v + sf1i[i];
							dfi[nbr_coef - i] = v - sf1i[i];
						}

						coef_index += d_nbr_coef;
					} while (coef_index < _length);

					/* Prepare to the next pass */
					{
						flt_t* const		temp_ptr = df;
						df = sf;
						sf = temp_ptr;
					}
				}
			}
		}

	/*______________________________________________
	 *
	 * Special cases
	 *______________________________________________
	 */

		/* 4-point FFT */
		else if (_nbr_bits == 2)
		{
			f[1] = x[0] - x[2];
			f[3] = x[1] - x[3];

			const flt_t			b_0 = x[0] + x[2];
			const flt_t			b_2 = x[1] + x[3];

			f[0] = b_0 + b_2;
			f[2] = b_0 - b_2;
		}

		/* 2-point FFT */
		else if (_nbr_bits == 1)
		{
			f[0] = x[0] + x[1];
			f[1] = x[0] - x[1];
		}

		/* 1-point FFT */
		else
		{
			f[0] = x[0];
		}
	}


	void	do_ifft(const flt_t f[], flt_t x[]) const
	{

	/*______________________________________________
	 *
	 * General case
	 *______________________________________________
	 */

		if (_nbr_bits > 2)
		{
			flt_t* sf = const_cast <flt_t*> (f);
			flt_t* df;
			flt_t* df_temp;

			if (_nbr_bits & 1)
			{
				df = _buffer_ptr;
				df_temp = x;
			}
			else
			{
				df = x;
				df_temp = _buffer_ptr;
			}

			/* Do the transformation in several pass */
			{
				int			pass;
				long			nbr_coef;
				long			h_nbr_coef;
				long			d_nbr_coef;
				long			coef_index;

				/* First pass */
				for (pass = _nbr_bits - 1; pass >= 3; --pass)
				{
					coef_index = 0;
					nbr_coef = 1 << pass;
					h_nbr_coef = nbr_coef >> 1;
					d_nbr_coef = nbr_coef << 1;
					const flt_t* const cos_ptr = _trigo_lut.get_ptr(pass);
					do
					{
						long				i;
						const flt_t* const sfr = sf + coef_index;
						const flt_t* const sfi = sfr + nbr_coef;
						flt_t* const df1r = df + coef_index;
						flt_t* const df2r = df1r + nbr_coef;

						/* Extreme coefficients are always real */
						df1r[0] = sfr[0] + sfi[0];		// + sfr [nbr_coef]
						df2r[0] = sfr[0] - sfi[0];		// - sfr [nbr_coef]
						df1r[h_nbr_coef] = sfr[h_nbr_coef] * 2;
						df2r[h_nbr_coef] = sfi[h_nbr_coef] * 2;

						/* Others are conjugate complex numbers */
						flt_t* const	df1i = df1r + h_nbr_coef;
						flt_t* const	df2i = df1i + nbr_coef;
						for (i = 1; i < h_nbr_coef; ++i)
						{
							df1r[i] = sfr[i] + sfi[-i];		// + sfr [nbr_coef - i]
							df1i[i] = sfi[i] - sfi[nbr_coef - i];

							const flt_t		c = cos_ptr[i];					// cos (i*PI/nbr_coef);
							const flt_t		s = cos_ptr[h_nbr_coef - i];	// sin (i*PI/nbr_coef);
							const flt_t		vr = sfr[i] - sfi[-i];		// - sfr [nbr_coef - i]
							const flt_t		vi = sfi[i] + sfi[nbr_coef - i];

							df2r[i] = vr * c + vi * s;
							df2i[i] = vi * c - vr * s;
						}

						coef_index += d_nbr_coef;
					} while (coef_index < _length);

					/* Prepare to the next pass */
					if (pass < _nbr_bits - 1)
					{
						flt_t* const	temp_ptr = df;
						df = sf;
						sf = temp_ptr;
					}
					else
					{
						sf = df;
						df = df_temp;
					}
				}

				/* Antepenultimate pass */
				{
					const flt_t		sqrt2_2 = _sqrt2_2;
					coef_index = 0;
					do
					{
						df[coef_index] = sf[coef_index] + sf[coef_index + 4];
						df[coef_index + 4] = sf[coef_index] - sf[coef_index + 4];
						df[coef_index + 2] = sf[coef_index + 2] * 2;
						df[coef_index + 6] = sf[coef_index + 6] * 2;

						df[coef_index + 1] = sf[coef_index + 1] + sf[coef_index + 3];
						df[coef_index + 3] = sf[coef_index + 5] - sf[coef_index + 7];

						const flt_t		vr = sf[coef_index + 1] - sf[coef_index + 3];
						const flt_t		vi = sf[coef_index + 5] + sf[coef_index + 7];

						df[coef_index + 5] = (vr + vi) * sqrt2_2;
						df[coef_index + 7] = (vi - vr) * sqrt2_2;

						coef_index += 8;
					} while (coef_index < _length);
				}

				/* Penultimate and last pass at once */
				{
					coef_index = 0;
					const long* bit_rev_lut_ptr = _bit_rev_lut.get_ptr();
					const flt_t* sf2 = df;
					do
					{
						{
							const flt_t		b_0 = sf2[0] + sf2[2];
							const flt_t		b_2 = sf2[0] - sf2[2];
							const flt_t		b_1 = sf2[1] * 2;
							const flt_t		b_3 = sf2[3] * 2;

							x[bit_rev_lut_ptr[0]] = b_0 + b_1;
							x[bit_rev_lut_ptr[1]] = b_0 - b_1;
							x[bit_rev_lut_ptr[2]] = b_2 + b_3;
							x[bit_rev_lut_ptr[3]] = b_2 - b_3;
						}
						{
							const flt_t		b_0 = sf2[4] + sf2[6];
							const flt_t		b_2 = sf2[4] - sf2[6];
							const flt_t		b_1 = sf2[5] * 2;
							const flt_t		b_3 = sf2[7] * 2;

							x[bit_rev_lut_ptr[4]] = b_0 + b_1;
							x[bit_rev_lut_ptr[5]] = b_0 - b_1;
							x[bit_rev_lut_ptr[6]] = b_2 + b_3;
							x[bit_rev_lut_ptr[7]] = b_2 - b_3;
						}

						sf2 += 8;
						coef_index += 8;
						bit_rev_lut_ptr += 8;
					} while (coef_index < _length);
				}
			}
		}

	/*______________________________________________
	 *
	 * Special cases
	 *______________________________________________
	 */

		/* 4-point IFFT */
		else if (_nbr_bits == 2)
		{
			const flt_t		b_0 = f[0] + f[2];
			const flt_t		b_2 = f[0] - f[2];

			x[0] = b_0 + f[1] * 2;
			x[2] = b_0 - f[1] * 2;
			x[1] = b_2 + f[3] * 2;
			x[3] = b_2 - f[3] * 2;
		}

		/* 2-point IFFT */
		else if (_nbr_bits == 1)
		{
			x[0] = f[0] + f[1];
			x[1] = f[0] - f[1];
		}

		/* 1-point IFFT */
		else
		{
			x[0] = f[0];
		}
	}



	void	rescale(flt_t x[]) const
	{
		const flt_t		mul = flt_t(1.0 / _length);
		long				i = _length - 1;

		do
		{
			x[i] *= mul;
			--i;
		} while (i >= 0);
	}




/*\\\ PRIVATE \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

private:

	/* Bit-reversed look-up table nested class */
	class BitReversedLUT
	{
	public:
		
		BitReversedLUT(const int nbr_bits)
		{
			long				length;
			long				cnt;
			long				br_index;
			long				bit;

			length = 1L << nbr_bits;
			_ptr = new long[length];

			br_index = 0;
			_ptr[0] = 0;
			for (cnt = 1; cnt < length; ++cnt)
			{
				/* ++br_index (bit reversed) */
				bit = length >> 1;
				while (((br_index ^= bit) & bit) == 0)
				{
					bit >>= 1;
				}

				_ptr[cnt] = br_index;
			}
		}



		~BitReversedLUT(void)
		{
			delete[] _ptr;
			_ptr = 0;
		}

		const long* get_ptr() const
		{
			return (_ptr);
		}
	private:
		long* _ptr;
	};

	/* Trigonometric look-up table nested class */
	class	TrigoLUT
	{
	public:
		TrigoLUT(const int nbr_bits)
		{
			long		total_len;

			_ptr = 0;
			if (nbr_bits > 3)
			{
				total_len = (1L << (nbr_bits - 1)) - 4;
				_ptr = new flt_t[total_len];
				const double	PIX = atan(1.0f) * 4;
				for (int level = 3; level < nbr_bits; ++level)
				{
					const long		level_len = 1L << (level - 1);
					flt_t* const	level_ptr = const_cast<flt_t*> (get_ptr(level));
					const double	mul = PIX / (level_len << 1);

					for (long i = 0; i < level_len; ++i)
					{
						level_ptr[i] = (flt_t)cos(i * mul);
					}
				}
			}
		}


		~TrigoLUT(void)
		{
			delete[] _ptr;
			_ptr = 0;
		}

		const flt_t* get_ptr(const int level) const
		{
#ifdef _WIN64
			return (_ptr + (1i64 << (level - 1)) - 4);
#else
			return (_ptr + (1L << (level - 1)) - 4);
#endif
		};
	private:
		flt_t* _ptr;
	};

	const BitReversedLUT	_bit_rev_lut;
	const TrigoLUT	_trigo_lut;
	const flt_t		_sqrt2_2;
	const long		_length;
	const int		_nbr_bits;
	flt_t* _buffer_ptr;



/*\\\ FORBIDDEN MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

private:

	FFTReal(const FFTReal& other);
	const FFTReal& operator = (const FFTReal& other);
	int				operator == (const FFTReal& other);
	int				operator != (const FFTReal& other);
};



#if defined (_MSC_VER)
#pragma pack (pop)
#endif	// _MSC_VER



#endif	// FFTReal_HEADER_INCLUDED

#undef FFTReal_CURRENT_HEADER



/*****************************************************************************

	LEGAL

	Source code may be freely used for any purpose, including commercial
	applications. Programs must display in their "About" dialog-box (or
	documentation) a text telling they use these routines by Laurent de Soras.
	Modified source code can be distributed, but modifications must be clearly
	indicated.

	CONTACT

	Laurent de Soras
	92 avenue Albert 1er
	92500 Rueil-Malmaison
	France

	ldesoras@club-internet.fr

*****************************************************************************/



/*\\\ EOF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/

/*****************************************************************************
*                                                                            *
*       DIGITAL SIGNAL PROCESSING TOOLS                                      *
*       Version 1.03, 2001/06/15                                             *
*       (c) 1999 - Laurent de Soras                                          *
*                                                                            *
*       FFTReal.cpp                                                          *
*       Fourier transformation of real number arrays.                        *
*       Portable ISO C++                                                     *
*                                                                            *
* Tab = 3                                                                    *
*****************************************************************************/



/*\\\ INCLUDE FILES \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/





#if defined (_MSC_VER)
#pragma pack (push, 8)
#endif	// _MSC_VER



/*\\\ PUBLIC MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/



/*==========================================================================*/
/*      Name: Constructor                                                   */
/*      Input parameters:                                                   */
/*        - length: length of the array on which we want to do a FFT.       */
/*                  Range: power of 2 only, > 0.                            */
/*      Throws: std::bad_alloc, anything                                    */
/*==========================================================================*/

/*==========================================================================*/
/*      Name: Destructor                                                    */
/*==========================================================================*/


/*==========================================================================*/
/*      Name: do_fft                                                        */
/*      Description: Compute the FFT of the array.                          */
/*      Input parameters:                                                   */
/*        - x: pointer on the source array (time).                          */
/*      Output parameters:                                                  */
/*        - f: pointer on the destination array (frequencies).              */
/*             f [0...length(x)/2] = real values,                           */
/*             f [length(x)/2+1...length(x)-1] = imaginary values of        */
/*               coefficents 1...length(x)/2-1.                             */
/*      Throws: Nothing                                                     */
/*==========================================================================*/




/*==========================================================================*/
/*      Name: do_ifft                                                       */
/*      Description: Compute the inverse FFT of the array. Notice that      */
/*                   IFFT (FFT (x)) = x * length (x). Data must be          */
/*                   post-scaled.                                           */
/*      Input parameters:                                                   */
/*        - f: pointer on the source array (frequencies).                   */
/*             f [0...length(x)/2] = real values,                           */
/*             f [length(x)/2+1...length(x)] = imaginary values of          */
/*               coefficents 1...length(x)-1.                               */
/*      Output parameters:                                                  */
/*        - x: pointer on the destination array (time).                     */
/*      Throws: Nothing                                                     */
/*==========================================================================*/


/*==========================================================================*/
/*      Name: rescale                                                       */
/*      Description: Scale an array by divide each element by its length.   */
/*                   This function should be called after FFT + IFFT.       */
/*      Input/Output parameters:                                            */
/*        - x: pointer on array to rescale (time or frequency).             */
/*      Throws: Nothing                                                     */
/*==========================================================================*/




/*\\\ NESTED CLASS MEMBER FUNCTIONS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/



/*==========================================================================*/
/*      Name: Constructor                                                   */
/*      Input parameters:                                                   */
/*        - nbr_bits: number of bits of the array on which we want to do a  */
/*                    FFT. Range: > 0                                       */
/*      Throws: std::bad_alloc                                              */
/*==========================================================================*/



/*==========================================================================*/
/*      Name: Destructor                                                    */
/*==========================================================================*/




/*==========================================================================*/
/*      Name: Constructor                                                   */
/*      Input parameters:                                                   */
/*        - nbr_bits: number of bits of the array on which we want to do a  */
/*                    FFT. Range: > 0                                       */
/*      Throws: std::bad_alloc, anything                                    */
/*==========================================================================*/



/*==========================================================================*/
/*      Name: Destructor                                                    */
/*==========================================================================*/




#if defined (_MSC_VER)
#pragma pack (pop)
#endif	// _MSC_VER



/*****************************************************************************

	LEGAL

	Source code may be freely used for any purpose, including commercial
	applications. Programs must display in their "About" dialog-box (or
	documentation) a text telling they use these routines by Laurent de Soras.
	Modified source code can be distributed, but modifications must be clearly
	indicated.

	CONTACT

	Laurent de Soras
	92 avenue Albert 1er
	92500 Rueil-Malmaison
	France

	ldesoras@club-internet.fr

*****************************************************************************/



/*\\\ EOF \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


template <typename T = float> class QuickFFT3
{
public:
	std::shared_ptr<FFTReal> fx = 0;
	int PrevS = 0;
	std::vector<T> b1;
	std::vector<T> b2;
	std::vector<T> mags;

	void Reset()
	{
		fx = 0;
		PrevS = 0;
	}

	void Prepare(T* in,size_t S)
	{
		if (PrevS != S)
		{
			fx = std::make_shared<FFTReal>(S);
			PrevS = S;
		}
		b1.clear();
		b2.clear();
		mags.clear();
		b1.resize(S);
		b2.resize(S);
		mags.resize(S);
		memcpy(b1.data(), in, S * sizeof(T));
	}

	void PrepareI(T* re, T* im, size_t S)
	{
		if (PrevS != S)
		{
			fx = std::make_shared<FFTReal>(S);
			PrevS = S;
		}
		b1.clear();
		b2.clear();
		mags.clear();
		b1.resize(S);
		b2.resize(S);
		mags.resize(S);

		size_t hs = S / 2;
		memcpy(b1.data(), re, S / 2);
		memcpy(b1.data() + hs, im, S / 2);
	}

	void ITransform()
	{
		b1.clear();
		b1.resize(b2.size());
		fx->do_ifft(b2.data(), b1.data());
	}

	void Transform(float *re = 0,float*im = 0)
	{
		size_t hs = b1.size() / 2;

		fx->do_fft(b2.data(), b1.data());
		if (re && im)
		{
			memcpy(re, b2.data(), hs * sizeof(T));
			memcpy(im, b2.data() + hs, hs * sizeof(T));
		}
	}

	void Mags()
	{
		int hs = b2.size() / 2;
		float* r = mags.data();
		for (int i = 0; i < b2.size() / 2; i++)
		{
			r[i] = sqrt(b2[i] * b2[i] + b2[i + hs] * b2[i + hs]);
			r[hs - i - 1] = r[i];
		}
	}

};


template <typename T = float> class QuickFFT2
{
public:

	std::shared_ptr<FFTReal> fx = 0;
	int PrevS = 0;

	void Reset()
	{
		fx = 0;
		PrevS = 0;
	}

	int Prepare(T* in, int S)
	{
		while (S > 0 && (S & (S - 1)) != 0)
		{
			S--;
		}
		if (PrevS != S)
		{
			fx = std::make_shared<FFTReal>(S);
			PrevS = S;
		}
		b1.clear();
		b1.resize(S);
		memcpy(b1.data(), in, S * sizeof(T));
		b2.clear();
		b2.resize(S);
		return S;
	}

	void Mags(float* r)
	{
		int hs = b2.size() / 2;
		for (int i = 0; i < b2.size() / 2; i++)
		{
			r[i] = sqrt(b2[i] * b2[i] + b2[i + hs] * b2[i + hs]);
			r[hs - i - 1] = r[i];
		}
	}

	T* Transform()
	{
		fx->do_fft(b2.data(), b1.data());
		return b2.data();
	}
	T* ITransform()
	{
		fx->do_ifft(b2.data(), b1.data());
		return b1.data();
	}

	T GetMaxIntensity()
	{
		int iM = -1;
		T vM = 0.0f;
		for (int i = 0; i < b2.size() / 2; i++)
		{
			T in = b2[i];
			if (in > vM)
			{
				iM = i;
				vM = in;
			}
		}
		return vM;
	}

	T GetF(int sr,int i)
	{
		T x = (T)(sr * i);
		return (T)(x / (T)b2.size());
	}


	std::vector<T> b1;
	std::vector<T> b2;
};

template <typename T = float> class QuickFFT
{
public:

	QuickFFT()
	{
	}

	QuickFFT(int SR, T* d, int s)
	{
		sr = SR;
		dd = d;
		ss = s;
	}

	void ReInput(int SR, T* data, int s)
	{
		sr = SR;
		dd = data;
		ss = s;
	}

	~QuickFFT()
	{
	}

	void Transform()
	{
		nP = ss;
		while (!IsPOfTwo(nP))
			nP++;

		AllIntensities.Resize(nP);
		ReIn.Resize(nP);
		ReOut.Resize(nP);
		ImIn.Resize(nP);
		ImOut.Resize(nP);

		AllIntensities.clear();
		memcpy(ReIn, dd, ss * sizeof(T));
		ImIn.clear();
		ReOut.clear();
		ImOut.clear();

		_fft<T>(nP, 0, ReIn, ImIn, ReOut, ImOut);
	}

	void ITransform()
	{
		_fft<T>(nP, 1, ReOut, ImOut, ReIn, ImIn);
	}

	T GetFR()
	{
		int savei = 0;
		T saveI = 0;
		for (int i = 0; i < nP / 2; i++)
		{
			T I = GetI(i);
			if (saveI < I)
			{
				saveI = I; // save maximum intensity
				savei = i; // and its index
				maxI = I;
				maxIi = i;
			}
		}
		return GetF(savei);
	}

	T GetFRD(T fr, int dom = 70)
	{
	// Dominant

	/*
	Formula by www.cipoo.net
	Find the first maximum that is at least dom% of the MaxInt
	*and* at least 6 commas before it

	*/
		T PerInt = (dom * maxI) / 100;
		int DomInt = -1;
		for (int i = 0; i < nP / 2; i++)
		{
			T I = GetI(i);
			T intFr = GetF(i);
			if (I < PerInt)
				continue; // intensity too low

			// intensity high! frequency 
			//int S = 72*(_F_)(log10((_F_)((_F_)fr / (_F_)intFr))/log10(2.0f));
			if (intFr == 0)
				continue;

			int S = abs(HzDiff(fr, intFr));
			if (S < 50) // semitone ?
			{
				break;
			}

			DomInt = i;
			break;
		}

	// If maximum is not found, thats all.
		if (DomInt == -1)
			return fr;

		// Found, get closer to max
		while (GetI(DomInt) < GetI(DomInt + 1))
			DomInt++;


		T DomFr = GetF(DomInt);
		return DomFr;
	}

	T GetDominantCorrect()
	{
	// Get the array of intensities, calculate peaks 
		int iM = GetMaxIntensity();
		T vM = GetI(iM);

		// And go back from here, finding tops
		vector<int> IntensityPeaks;
		IntensityPeaks.push_back(iM);
		int GoingMode = 0; // 0 lowering, 1 highering
		T PrevV = vM;
		T PrevF = GetF(iM);

		for (int j = (iM - 1); j >= 0; j--)
		{
			T inte = GetI(j);
			if (inte <= PrevV && GoingMode == 0)
			{
				PrevV = inte;
				continue;
			}
			if (inte >= PrevV && GoingMode == 1)
			{
				PrevV = inte;
				continue;
			}

			int S = HzDiffAbs(GetF(j), PrevF);
			if (S > 500)
			{
				IntensityPeaks.push_back(j);
				PrevF = GetF(j);
			}
			if (GoingMode == 0)
				GoingMode = 1;
			else
				GoingMode = 0;
		}
		return 0;
	}

	int GetNP()
	{
		return nP;
	}

	int GetMaxIntensity()
	{
		int iM = -1;
		T vM = 0.0f;
		for (int i = 0; i < nP / 2; i++)
		{
			T in = GetI(i);
			if (in > vM)
			{
				iM = i;
				vM = in;
			}
		}
		return iM;
	}

	T GetI(int i)
	{
		if (AllIntensities[i])
			return AllIntensities[i];

		T y = 0;
		if (ImOut[i] == 0)
			return ReOut[i];

		y = sqrt((ReOut[i] * ReOut[i]) + (ImOut[i] * ImOut[i]));
		y /= ((T)sqrt((T)nP));

		AllIntensities[i] = y;
		return y;
	}

	T GetNormI(int i)
	{
		int m = GetMaxIntensity();
		return GetI(i) / GetI(m);
	}

	T GetIdB(int i)
	{
		T inte = GetI(i);
		return 20 * log10(inte);
	}

	T GetNormIdB(int i)
	{
		T inte = GetNormI(i);
		return 20 * log10(inte);
	}

	T GetF(int i)
	{
		T x = sr * i;
		return (T)(x / (T)nP);
	}

/*	void GetAllFrequenciesSortedByIntensity(std::vector<PEAK>& p)
	{
	// p.p is the frequency, p.v is the intensity
		p.resize(nP / 2);
		for (int i = 0; i < nP / 2; i++)
		{
			p[i].p = GetF(i);
			if (i == 40)
				nop();
			T inte = GetI(i);
			//		inte = 20*log10(inte);
			p[i].v = inte;
			nop();
		}
		std::sort(p.begin(), p.end());
		std::reverse(p.begin(), p.end());
	}
*/
	std::vector<T>& GetReIn() { return ReIn; }
	std::vector<T>& GetImIn() { return ImIn; }
	std::vector<T>& GetReOut() { return ReOut; }
	std::vector<T>& GetImOut() { return ImOut; }



private:


	T* dd;
	int ss;
	int sr;
	int nP;
	std::vector<T> AllIntensities;
	std::vector<T> ReIn;
	std::vector<T> ReOut;
	std::vector<T> ImIn;
	std::vector<T> ImOut;
	T maxI;
	int maxIi;


};


#endif