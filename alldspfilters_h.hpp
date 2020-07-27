/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_COMMON_H
#define DSPFILTERS_COMMON_H

//
// This must be the first file included in every DspFilters header and source
//

#ifdef _MSC_VER
#  pragma warning (disable: 4100)
#endif

//#include <assert.h>
#include <stdlib.h>

#include <cassert>
#include <cfloat>
#include <cmath>
#include <complex>
#include <cstring>
#include <string>
#include <limits>
#include <vector>
#include <algorithm>

namespace tr1 = std;


#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_MATHSUPPLEMENT_H
#define DSPFILTERS_MATHSUPPLEMENT_H


namespace Dsp {

    const double doublePi = 3.1415926535897932384626433832795028841971;
    const double doublePi_2 = 1.5707963267948966192313216916397514420986;
    const double doubleLn2 = 0.69314718055994530941723212145818;//?????
    const double doubleLn10 = 2.3025850929940456840179914546844;//??????

    typedef std::complex<double> complex_t;
    typedef std::pair<complex_t, complex_t> complex_pair_t;

    template<typename Real>
    inline std::complex<Real> solve_quadratic_1(Real a, Real b, Real c)
    {
        return (-b + sqrt(std::complex<Real>(b * b - 4 * a * c, 0))) / (2. * a);
    }

    template<typename Real>
    inline std::complex<Real> solve_quadratic_2(Real a, Real b, Real c)
    {
        return (-b - sqrt(std::complex<Real>(b * b - 4 * a * c, 0))) / (2. * a);
    }

    inline const complex_t infinity()
    {
        return complex_t(std::numeric_limits<double>::infinity());
    }

    inline const complex_t adjust_imag(const complex_t& c)
    {
        if (fabs(c.imag()) < 1e-30)
            return complex_t(c.real(), 0);
        else
            return c;
    }

    template <typename Ty, typename To>
    inline std::complex<Ty> addmul(const std::complex<Ty>& c,
        Ty v,
        const std::complex<To>& c1)
    {
        return std::complex <Ty>(
            c.real() + v * c1.real(), c.imag() + v * c1.imag());
    }

    template <typename Ty>
    inline std::complex<Ty> recip(const std::complex<Ty>& c)
    {
        Ty n = 1.0 / std::norm(c);

        return std::complex<Ty>(n * c.real(), n * c.imag());
    }

    template <typename Ty>
    inline Ty asinh(Ty x)
    {
        return log(x + std::sqrt(x * x + 1));
    }

    template <typename Ty>
    inline Ty acosh(Ty x)
    {
        return log(x + std::sqrt(x * x - 1));
    }

    template <typename Ty>
    inline bool is_nan(Ty v)
    {
        return !(v == v);
    }

    template <>
    inline bool is_nan<complex_t>(complex_t v)
    {
        return Dsp::is_nan(v.real()) || Dsp::is_nan(v.imag());
    }

    //------------------------------------------------------------------------------

    /*
     * Hack to prevent denormals
     *
     */

     //const double anti_denormal_vsa = 1e-16; // doesn't prevent denormals
     //const double anti_denormal_vsa = 0;
    const double anti_denormal_vsa = 1e-8;

    class DenormalPrevention
    {
    public:
        DenormalPrevention()
            : m_v(anti_denormal_vsa)
        {
        }

        // small alternating current
        inline double ac()
        {
            return m_v = -m_v;
        }

        // small direct current
        static inline double dc()
        {
            return anti_denormal_vsa;
        }

    private:
        double m_v;
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_TYPES_H
#define DSPFILTERS_TYPES_H


namespace Dsp {

    // A conjugate or real pair
    struct ComplexPair : complex_pair_t
    {
        ComplexPair()
        {
        }

        explicit ComplexPair(const complex_t& c1)
            : complex_pair_t(c1, 0.)
        {
            assert(isReal());
        }

        ComplexPair(const complex_t& c1,
            const complex_t& c2)
            : complex_pair_t(c1, c2)
        {
        }

        bool isConjugate() const
        {
            return second == std::conj(first);
        }

        bool isReal() const
        {
            return first.imag() == 0 && second.imag() == 0;
        }

        // Returns true if this is either a conjugate pair,
        // or a pair of reals where neither is zero.
        bool isMatchedPair() const
        {
            if (first.imag() != 0)
                return second == std::conj(first);
            else
                return second.imag() == 0 &&
                second.real() != 0 &&
                first.real() != 0;
        }

        bool is_nan() const
        {
            return Dsp::is_nan(first) || Dsp::is_nan(second);
        }
    };

    // A pair of pole/zeros. This fits in a biquad (but is missing the gain)
    struct PoleZeroPair
    {
        ComplexPair poles;
        ComplexPair zeros;

        PoleZeroPair() { }

        // single pole/zero
        PoleZeroPair(const complex_t& p, const complex_t& z)
            : poles(p), zeros(z)
        {
        }

        // pole/zero pair
        PoleZeroPair(const complex_t& p1, const complex_t& z1,
            const complex_t& p2, const complex_t& z2)
            : poles(p1, p2)
            , zeros(z1, z2)
        {
        }

        bool isSinglePole() const
        {
            return poles.second == 0. && zeros.second == 0.;
        }

        bool is_nan() const
        {
            return poles.is_nan() || zeros.is_nan();
        }
    };

    // Identifies the general class of filter
    enum Kind
    {
        kindLowPass,
        kindHighPass,
        kindBandPass,
        kindBandStop,
        kindLowShelf,
        kindHighShelf,
        kindBandShelf,
        kindOther
    };

}

#endif

/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_BIQUAD_H
#define DSPFILTERS_BIQUAD_H


namespace Dsp {

    struct BiquadPoleState;

    /*
     * Holds coefficients for a second order Infinite Impulse Response
     * digital filter. This is the building block for all IIR filters.
     *
     */

     // Factored interface to prevent outsiders from fiddling
    class BiquadBase
    {
    public:
        template <class StateType>
        struct State : StateType, private DenormalPrevention
        {
            template <typename Sample>
            inline Sample process(const Sample in, const BiquadBase& b)
            {
                return static_cast<Sample> (StateType::process1(in, b, ac()));
            }
        };

    public:
        // Calculate filter response at the given normalized frequency.
        complex_t response(double normalizedFrequency) const;

        std::vector<PoleZeroPair> getPoleZeros() const;

        double getA0() const { return m_a0; }
        double getA1() const { return m_a1 * m_a0; }
        double getA2() const { return m_a2 * m_a0; }
        double getB0() const { return m_b0 * m_a0; }
        double getB1() const { return m_b1 * m_a0; }
        double getB2() const { return m_b2 * m_a0; }

        // Process a block of samples in the given form
        template <class StateType, typename Sample>
        void process(int numSamples, Sample* dest, StateType& state) const
        {
            while (--numSamples >= 0) {
                *dest = state.process(*dest, *this);
                dest++;
            }
        }

    protected:
        //
        // These are protected so you can't mess with RBJ biquads
        //

        void setCoefficients(double a0, double a1, double a2,
            double b0, double b1, double b2);

        void setOnePole(complex_t pole, complex_t zero);

        void setTwoPole(complex_t pole1, complex_t zero1,
            complex_t pole2, complex_t zero2);

        void setPoleZeroPair(const PoleZeroPair& pair)
        {
            if (pair.isSinglePole())
                setOnePole(pair.poles.first, pair.zeros.first);
            else
                setTwoPole(pair.poles.first, pair.zeros.first,
                    pair.poles.second, pair.zeros.second);
        }

        void setPoleZeroForm(const BiquadPoleState& bps);

        void setIdentity();

        void applyScale(double scale);

    public:
        double m_a0;
        double m_a1;
        double m_a2;
        double m_b1;
        double m_b2;
        double m_b0;
    };

    //------------------------------------------------------------------------------

    // Expresses a biquad as a pair of pole/zeros, with gain
    // values so that the coefficients can be reconstructed precisely.
    struct BiquadPoleState : PoleZeroPair
    {
        BiquadPoleState() { }

        explicit BiquadPoleState(const BiquadBase& s);

        double gain;
    };

    // More permissive interface for fooling around
    class Biquad : public BiquadBase
    {
    public:
        Biquad();

        explicit Biquad(const BiquadPoleState& bps);

    public:
        // Process a block of samples, interpolating from the old section's coefficients
        // to this section's coefficients, over numSamples. This implements smooth
        // parameter changes.

        template <class StateType, typename Sample>
        void smoothProcess1(int numSamples,
            Sample* dest,
            StateType& state,
            Biquad sectionPrev) const
        {
            double t = 1. / numSamples;
            double da1 = (m_a1 - sectionPrev.m_a1) * t;
            double da2 = (m_a2 - sectionPrev.m_a2) * t;
            double db0 = (m_b0 - sectionPrev.m_b0) * t;
            double db1 = (m_b1 - sectionPrev.m_b1) * t;
            double db2 = (m_b2 - sectionPrev.m_b2) * t;

            while (--numSamples >= 0)
            {
                sectionPrev.m_a1 += da1;
                sectionPrev.m_a2 += da2;
                sectionPrev.m_b0 += db0;
                sectionPrev.m_b1 += db1;
                sectionPrev.m_b2 += db2;

                *dest = state.process(*dest, sectionPrev);
                dest++;
            }
        }

        // Process a block of samples, interpolating from the old section's pole/zeros
        // to this section's pole/zeros, over numSamples. The interpolation is done
        // in the z-plane using polar coordinates.
        template <class StateType, typename Sample>
        void smoothProcess2(int numSamples,
            Sample* dest,
            StateType& state,
            BiquadPoleState zPrev) const
        {
            BiquadPoleState z(*this);
            double t = 1. / numSamples;
            complex_t dp0 = (z.poles.first - zPrev.poles.first) * t;
            complex_t dp1 = (z.poles.second - zPrev.poles.second) * t;
            complex_t dz0 = (z.zeros.first - zPrev.zeros.first) * t;
            complex_t dz1 = (z.zeros.second - zPrev.zeros.second) * t;
            double dg = (z.gain - zPrev.gain) * t;

            while (--numSamples >= 0)
            {
                zPrev.poles.first += dp0;
                zPrev.poles.second += dp1;
                zPrev.zeros.first += dz0;
                zPrev.zeros.second += dz1;
                zPrev.gain += dg;

                *dest = state.process(*dest, Biquad(zPrev));
                dest++;
            }
        }

    public:
        // Export these as public

        void setOnePole(complex_t pole, complex_t zero)
        {
            BiquadBase::setOnePole(pole, zero);
        }

        void setTwoPole(complex_t pole1, complex_t zero1,
            complex_t pole2, complex_t zero2)
        {
            BiquadBase::setTwoPole(pole1, zero1, pole2, zero2);
        }

        void setPoleZeroPair(const PoleZeroPair& pair)
        {
            BiquadBase::setPoleZeroPair(pair);
        }

        void applyScale(double scale)
        {
            BiquadBase::applyScale(scale);
        }
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_PARAMS_H
#define DSPFILTERS_PARAMS_H

namespace Dsp {

    /*
     * System for abstracting parameterizable filter specifications.
     *
     * This provides a "GUI-friendly" layer to the filters. Note that
     * it is not necessary to use this layer, it is possible to instantiate
     * the filters and their associated processing state directly,
     * and bypass the overhead for this API if it is not needed.
     *
     */

     // Unique IDs to help identify parameters
    enum ParamID
    {
        idSampleRate,
        idFrequency,
        idQ,
        idBandwidth,
        idBandwidthHz,
        idGain,
        idSlope,
        idOrder,
        idRippleDb,
        idStopDb,
        idRolloff,

        idPoleRho,
        idPoleTheta,
        idZeroRho,
        idZeroTheta,

        idPoleReal,
        idZeroReal
    };

    enum
    {
        maxParameters = 8
    };

    struct Params
    {
        void clear()
        {
            for (int i = 0; i < maxParameters; ++i)
                value[i] = 0;
        }

        double& operator[] (int index)
        {
            return value[index];
        }

        const double& operator[] (int index) const
        {
            return value[index];
        }

        double value[maxParameters];
    };

    //
    // Provides meta-information about a filter parameter
    // to achieve run-time introspection.
    //
    class ParamInfo
    {
    public:
        typedef double (ParamInfo::* toControlValue_t) (double) const;
        typedef double (ParamInfo::* toNativeValue_t) (double) const;
        typedef std::string(ParamInfo::* toString_t) (double) const;

        // dont use this one
        ParamInfo(); // throws std::logic_error

        ParamInfo(ParamID id,
            const char* szLabel,
            const char* szName,
            double arg1,
            double arg2,
            double defaultNativeValue,
            toControlValue_t toControlValue_proc,
            toNativeValue_t toNativeValue_proc,
            toString_t toString_proc)
            : m_id(id)
            , m_szLabel(szLabel)
            , m_szName(szName)
            , m_arg1(arg1)
            , m_arg2(arg2)
            , m_defaultNativeValue(defaultNativeValue)
            , m_toControlValue(toControlValue_proc)
            , m_toNativeValue(toNativeValue_proc)
            , m_toString(toString_proc)
        {
        }

        // Used to identify well-known parameters (like cutoff frequency)
        ParamID getId() const
        {
            return m_id;
        }

        // Returns a short label suitable for placement on a control
        const char* getLabel() const
        {
            return m_szLabel;
        }

        // Returns the full name
        const char* getName() const
        {
            return m_szName;
        }

        double getDefaultValue() const
        {
            return m_defaultNativeValue;
        }

        //
        // Control value is always in the range [0..1]
        //
        double toControlValue(double nativeValue) const
        {
            return (this->*m_toControlValue) (nativeValue);
        }

        //
        // Native value is in filter-specific units. For example,
        // cutoff frequency would probably be in Hertz.
        //
        double toNativeValue(double controlValue) const
        {
            return (this->*m_toNativeValue) (controlValue);
        }

        std::string toString(double nativeValue) const
        {
            return (this->*m_toString) (nativeValue);
        }

        double clamp(double nativeValue) const;

        //
        // These routines are used as function pointers when
        // constructing the various ParamInfo used by filters
        //

        double Int_toControlValue(double nativeValue) const;
        double Int_toNativeValue(double controlValue) const;

        double Real_toControlValue(double nativeValue) const;
        double Real_toNativeValue(double controlValue) const;

        double Log_toControlValue(double nativeValue) const;
        double Log_toNativeValue(double controlValue) const;

        double Pow2_toControlValue(double nativeValue) const;
        double Pow2_toNativeValue(double controlValue) const;

        std::string Int_toString(double nativeValue) const;
        std::string Hz_toString(double nativeValue) const;
        std::string Real_toString(double nativeValue) const;
        std::string Db_toString(double nativeValue) const;

        //
        // Creates the specified ParamInfo
        //

        static ParamInfo defaultSampleRateParam();
        static ParamInfo defaultCutoffFrequencyParam();
        static ParamInfo defaultCenterFrequencyParam();
        static ParamInfo defaultQParam();
        static ParamInfo defaultBandwidthParam();
        static ParamInfo defaultBandwidthHzParam();
        static ParamInfo defaultGainParam();
        static ParamInfo defaultSlopeParam();
        static ParamInfo defaultRippleDbParam();
        static ParamInfo defaultStopDbParam();
        static ParamInfo defaultRolloffParam();
        static ParamInfo defaultPoleRhoParam();
        static ParamInfo defaultPoleThetaParam();
        static ParamInfo defaultZeroRhoParam();
        static ParamInfo defaultZeroThetaParam();
        static ParamInfo defaultPoleRealParam();
        static ParamInfo defaultZeroRealParam();

    private:
        ParamID m_id;
        const char* m_szLabel;
        const char* m_szName;
        double m_arg1;
        double m_arg2;
        double m_defaultNativeValue;
        toControlValue_t m_toControlValue;
        toNativeValue_t m_toNativeValue;
        toString_t m_toString;
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_STATE_H
#define DSPFILTERS_STATE_H

#include <stdexcept>

namespace Dsp {

    /*
     * Various forms of state information required to
     * process channels of actual sample data.
     *
     */

     //------------------------------------------------------------------------------

     /*
      * State for applying a second order section to a sample using Direct Form I
      *
      * Difference equation:
      *
      *  y[n] = (b0/a0)*x[n] + (b1/a0)*x[n-1] + (b2/a0)*x[n-2]
      *                      - (a1/a0)*y[n-1] - (a2/a0)*y[n-2]
      */
    class DirectFormI
    {
    public:
        DirectFormI()
        {
            reset();
        }

        void reset()
        {
            m_x1 = 0;
            m_x2 = 0;
            m_y1 = 0;
            m_y2 = 0;
        }

        template <typename Sample>
        inline Sample process1(const Sample in,
            const BiquadBase& s,
            const double vsa) // very small amount
        {
            double out = s.m_b0 * in + s.m_b1 * m_x1 + s.m_b2 * m_x2
                - s.m_a1 * m_y1 - s.m_a2 * m_y2
                + vsa;
            m_x2 = m_x1;
            m_y2 = m_y1;
            m_x1 = in;
            m_y1 = out;

            return static_cast<Sample> (out);
        }

    protected:
        double m_x2; // x[n-2]
        double m_y2; // y[n-2]
        double m_x1; // x[n-1]
        double m_y1; // y[n-1]
    };

    //------------------------------------------------------------------------------

    /*
     * State for applying a second order section to a sample using Direct Form II
     *
     * Difference equation:
     *
     *  v[n] =         x[n] - (a1/a0)*v[n-1] - (a2/a0)*v[n-2]
     *  y(n) = (b0/a0)*v[n] + (b1/a0)*v[n-1] + (b2/a0)*v[n-2]
     *
     */
    class DirectFormII
    {
    public:
        DirectFormII()
        {
            reset();
        }

        void reset()
        {
            m_v1 = 0;
            m_v2 = 0;
        }

        template <typename Sample>
        Sample process1(const Sample in,
            const BiquadBase& s,
            const double vsa)
        {
            double w = in - s.m_a1 * m_v1 - s.m_a2 * m_v2 + vsa;
            double out = s.m_b0 * w + s.m_b1 * m_v1 + s.m_b2 * m_v2;

            m_v2 = m_v1;
            m_v1 = w;

            return static_cast<Sample> (out);
        }

    private:
        double m_v1; // v[-1]
        double m_v2; // v[-2]
    };

    //------------------------------------------------------------------------------

    /*
     * Transposed Direct Form I and II
     * by lubomir i. ivanov (neolit123 [at] gmail)
     *
     * Reference:
     * http://www.kvraudio.com/forum/viewtopic.php?p=4430351
     *
     */

     // I think this one is broken
    class TransposedDirectFormI
    {
    public:
        TransposedDirectFormI()
        {
            reset();
        }

        void reset()
        {
            m_v = 0;
            m_s1 = 0;
            m_s1_1 = 0;
            m_s2 = 0;
            m_s2_1 = 0;
            m_s3 = 0;
            m_s3_1 = 0;
            m_s4 = 0;
            m_s4_1 = 0;
        }

        template <typename Sample>
        inline Sample process1(const Sample in,
            const BiquadBase& s,
            const double vsa)
        {
            double out;

            // can be: in += m_s1_1;
            m_v = in + m_s1_1;
            out = s.m_b0 * m_v + m_s3_1;
            m_s1 = m_s2_1 - s.m_a1 * m_v;
            m_s2 = -s.m_a2 * m_v;
            m_s3 = s.m_b1 * m_v + m_s4_1;
            m_s4 = s.m_b2 * m_v;

            m_s4_1 = m_s4;
            m_s3_1 = m_s3;
            m_s2_1 = m_s2;
            m_s1_1 = m_s1;

            return static_cast<Sample> (out);
        }

    private:
        double m_v;
        double m_s1;
        double m_s1_1;
        double m_s2;
        double m_s2_1;
        double m_s3;
        double m_s3_1;
        double m_s4;
        double m_s4_1;
    };

    //------------------------------------------------------------------------------

    class TransposedDirectFormII
    {
    public:
        TransposedDirectFormII()
        {
            reset();
        }

        void reset()
        {
            m_s1 = 0;
            m_s1_1 = 0;
            m_s2 = 0;
            m_s2_1 = 0;
        }

        template <typename Sample>
        inline Sample process1(const Sample in,
            const BiquadBase& s,
            const double vsa)
        {
            double out;

            out = m_s1_1 + s.m_b0 * in + vsa;
            m_s1 = m_s2_1 + s.m_b1 * in - s.m_a1 * out;
            m_s2 = s.m_b2 * in - s.m_a2 * out;
            m_s1_1 = m_s1;
            m_s2_1 = m_s2;

            return static_cast<Sample> (out);
        }

    private:
        double m_s1;
        double m_s1_1;
        double m_s2;
        double m_s2_1;
    };

    //------------------------------------------------------------------------------

    // Holds an array of states suitable for multi-channel processing
    template <int ChannelsX, class StateType>
    class ChannelsState
    {
        int Channels = 1;
    public:
        ChannelsState()
        {
            SetChannelNumber(ChannelsX);
        }
        void SetChannelNumber(int nch)
        {
            Channels = nch;
            m_state.resize(Channels);
        }

        const int getNumChannels() const
        {
            return Channels;
        }

        void reset()
        {
            for (int i = 0; i < Channels; ++i)
                m_state[i].reset();
        }

        StateType& operator[] (int index)
        {
            assert(index >= 0 && index < Channels);
            return m_state[index];
        }

        template <class Filter, typename Sample>
        void process(int numSamples,
            Sample* const* arrayOfChannels,
            Filter& filter)
        {
            for (int i = 0; i < Channels; ++i)
                filter.process(numSamples, arrayOfChannels[i], m_state[i]);
        }

    private:
        std::vector<StateType> m_state;
    };

    // Empty state, can't process anything
    template <class StateType>
    class ChannelsState <0, StateType>
    {
    public:
        const int getNumChannels() const
        {
            return 0;
        }

        void reset()
        {
            throw std::logic_error("attempt to reset empty ChannelState");
        }

        template <class FilterDesign, typename Sample>
        void process(int numSamples,
            Sample* const* arrayOfChannels,
            FilterDesign& filter)
        {
            throw std::logic_error("attempt to process empty ChannelState");
        }
    };

    //------------------------------------------------------------------------------

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_FILTER_H
#define DSPFILTERS_FILTER_H

namespace Dsp {

    /*
     * Filter
     *
     * Full abstraction of a digital IIR filter.
     * Supports run-time introspection and modulation of filter
     * parameters.
     *
     */
    class Filter
    {
    public:
        virtual ~Filter();

        virtual Kind getKind() const = 0;

        virtual const std::string getName() const = 0;

        virtual int getNumParams() const = 0;

        virtual ParamInfo getParamInfo(int index) const = 0;

        Params getDefaultParams() const;

        const Params& getParams() const
        {
            return m_params;
        }

        double getParam(int paramIndex) const
        {
            assert(paramIndex >= 0 && paramIndex <= getNumParams());
            return m_params[paramIndex];
        }

        void setParam(int paramIndex, double nativeValue)
        {
            assert(paramIndex >= 0 && paramIndex <= getNumParams());
            m_params[paramIndex] = nativeValue;
            doSetParams(m_params);
        }

        int findParamId(int paramId);

        void setParamById(int paramId, double nativeValue);

        void setParams(const Params& parameters)
        {
            m_params = parameters;
            doSetParams(parameters);
        }

        // This makes a best-effort to pick up the values
        // of matching parameters from another set. It uses
        // the ParamID information to make the match.
        void copyParamsFrom(Dsp::Filter const* other);

        virtual std::vector<PoleZeroPair> getPoleZeros() const = 0;

        virtual complex_t response(double normalizedFrequency) const = 0;

        virtual int getNumChannels() = 0;
        virtual void reset() = 0;
        virtual void process(int numSamples, float* const* arrayOfChannels) = 0;
        virtual void process(int numSamples, double* const* arrayOfChannels) = 0;

    protected:
        virtual void doSetParams(const Params& parameters) = 0;

    private:
        Params m_params;
    };

    //------------------------------------------------------------------------------

    /*
     * FilterDesign
     *
     * This container holds a filter Design (Gui-friendly layer) and
     * optionally combines it with the necessary state information to
     * process channel data.
     *
     */

     // Factored to reduce template instantiations
    template <class DesignClass>
    class FilterDesignBase : public Filter
    {
    public:
        Kind getKind() const
        {
            return m_design.getKind();
        }

        const std::string getName() const
        {
            return m_design.getName();
        }

        int getNumParams() const
        {
            return DesignClass::NumParams;
        }

        Params getDefaultParams() const
        {
            return m_design.getDefaultParams();
        }

        ParamInfo getParamInfo(int index) const
        {
            switch (index)
            {
            case 0: return m_design.getParamInfo_0();
            case 1: return m_design.getParamInfo_1();
            case 2: return m_design.getParamInfo_2();
            case 3: return m_design.getParamInfo_3();
            case 4: return m_design.getParamInfo_4();
            case 5: return m_design.getParamInfo_5();
            case 6: return m_design.getParamInfo_6();
            case 7: return m_design.getParamInfo_7();
            };

            return ParamInfo();
        }

        std::vector<PoleZeroPair> getPoleZeros() const
        {
            return m_design.getPoleZeros();
        }

        complex_t response(double normalizedFrequency) const
        {
            return m_design.response(normalizedFrequency);
        }

    protected:
        void doSetParams(const Params& parameters)
        {
            m_design.setParams(parameters);
        }

    protected:
        DesignClass m_design;
    };



    template <class DesignClass,
        int Channels = 0,
        class StateType = DirectFormII>
        class FilterDesign : public FilterDesignBase <DesignClass>
    {
    public:
        FilterDesign()
        {
        }

        int getNumChannels()
        {
            return Channels;
        }

        void reset()
        {
            m_state.reset();
        }

        void process(int numSamples, float* const* arrayOfChannels)
        {
            m_state.process(numSamples, arrayOfChannels,
                FilterDesignBase<DesignClass>::m_design);
        }

        void process(int numSamples, double* const* arrayOfChannels)
        {
            m_state.process(numSamples, arrayOfChannels,
                FilterDesignBase<DesignClass>::m_design);
        }

    protected:
        ChannelsState <Channels,
            typename DesignClass::template State <StateType> > m_state;
    };

    //------------------------------------------------------------------------------

    /*
     * This container combines a raw filter with state information
     * so it can process channels. In order to set up the filter you
     * must call a setup function directly. Smooth changes are
     * not supported, but this class has a smaller footprint.
     *
     */
    template <class FilterClass,
        int ChannelsX = 0,
        class StateType = DirectFormII>
        class SimpleFilter : public FilterClass
    {
        int Channels;
    public:

        SimpleFilter() : FilterClass()
        {
            Channels = ChannelsX;
        }
        int getNumChannels()
        {
            return Channels;
        }
        void setNumChannels(int nch)
        {
            Channels = nch;
            m_state.SetChannelNumber(nch);
        }

        void reset()
        {
            m_state.reset();
        }

        template <typename Sample>
        void process(int numSamples, Sample* const* arrayOfChannels)
        {
            m_state.process(numSamples, arrayOfChannels, *((FilterClass*)this));
        }

    protected:
        ChannelsState <ChannelsX,
            typename FilterClass::template State <StateType> > m_state;
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_LAYOUT_H
#define DSPFILTERS_LAYOUT_H


namespace Dsp {

    //
    // Describes a filter as a collection of poles and zeros along with
    // normalization information to achieve a specified gain at a specified
    // frequency. The poles and zeros may lie either in the s or the z plane.
    //

    // Base uses pointers to reduce template instantiations
    class LayoutBase
    {
    public:
        LayoutBase()
            : m_numPoles(0)
            , m_maxPoles(0)
        {
        }

        LayoutBase(int maxPoles, PoleZeroPair* pairs)
            : m_numPoles(0)
            , m_maxPoles(maxPoles)
            , m_pair(pairs)
        {
        }

        void setStorage(const LayoutBase& other)
        {
            m_numPoles = 0;
            m_maxPoles = other.m_maxPoles;
            m_pair = other.m_pair;
        }

        void reset()
        {
            m_numPoles = 0;
        }

        int getNumPoles() const
        {
            return m_numPoles;
        }

        int getMaxPoles() const
        {
            return m_maxPoles;
        }

        void add(const complex_t& pole, const complex_t& zero)
        {
            assert(!(m_numPoles & 1)); // single comes last
            assert(!Dsp::is_nan(pole));
            m_pair[m_numPoles / 2] = PoleZeroPair(pole, zero);
            ++m_numPoles;
        }

        void addPoleZeroConjugatePairs(const complex_t pole,
            const complex_t zero)
        {
            assert(!(m_numPoles & 1)); // single comes last
            assert(!Dsp::is_nan(pole));
            m_pair[m_numPoles / 2] = PoleZeroPair(
                pole, zero, std::conj(pole), std::conj(zero));
            m_numPoles += 2;
        }

        void add(const ComplexPair& poles, const ComplexPair& zeros)
        {
            assert(!(m_numPoles & 1)); // single comes last
            assert(poles.isMatchedPair());
            assert(zeros.isMatchedPair());
            m_pair[m_numPoles / 2] = PoleZeroPair(poles.first, zeros.first,
                poles.second, zeros.second);
            m_numPoles += 2;
        }

        const PoleZeroPair& getPair(int pairIndex) const
        {
            assert(pairIndex >= 0 && pairIndex < (m_numPoles + 1) / 2);
            return m_pair[pairIndex];
        }

        const PoleZeroPair& operator[] (int pairIndex) const
        {
            return getPair(pairIndex);
        }

        double getNormalW() const
        {
            return m_normalW;
        }

        double getNormalGain() const
        {
            return m_normalGain;
        }

        void setNormal(double w, double g)
        {
            m_normalW = w;
            m_normalGain = g;
        }

    private:
        int m_numPoles;
        int m_maxPoles;
        PoleZeroPair* m_pair;
        double m_normalW;
        double m_normalGain;
    };

    //------------------------------------------------------------------------------

    // Storage for Layout
    template <int MaxPoles>
    class Layout
    {
    public:
        operator LayoutBase ()
        {
            return LayoutBase(MaxPoles, m_pairs);
        }

    private:
        PoleZeroPair m_pairs[(MaxPoles + 1) / 2];
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_ROOTFINDER_H
#define DSPFILTERS_ROOTFINDER_H


namespace Dsp {

    //
    // Finds the complex roots of the given polynomial with
    // complex-valued coefficients using a numerical method.
    //

    class RootFinderBase
    {
    public:
        struct Array
        {
            Array(int max, complex_t* values)
                // : m_max (max)
                // , m_values (values)
            {
            }

            //complex_t& operator[] (int index)
            //{
            //};
        };

        //
        // Find roots of polynomial f(x)=a[0]+a[1]*x+a[2]*x^2...+a[degree]*x^degree
        // The input coefficients are set using coef()[].
        // The solutions are placed in roots.
        //
        void solve(int degree,
            bool polish = true,
            bool doSort = true);

        // Evaluates the polynomial at x
        complex_t eval(int degree,
            const complex_t& x);

        // Direct access to the input coefficient array of size degree+1.
        complex_t* coef()
        {
            return m_a;
        }

        // Direct access to the resulting roots array of size degree
        complex_t* root()
        {
            return m_root;
        }

        // sort the roots by descending imaginary part
        void sort(int degree);

    private:
        // Improves x as a root using Laguerre's method.
        // The input coefficient array has degree+1 elements.
        void laguerre(int degree,
            complex_t a[],
            complex_t& x,
            int& its);

    protected:
        int m_maxdegree;
        complex_t* m_a;		// input coefficients (m_maxdegree+1 elements)
        complex_t* m_ad;	// copy of deflating coefficients
        complex_t* m_root; // array of roots (maxdegree elements)
    };

    //------------------------------------------------------------------------------

    template<int maxdegree>
    struct RootFinder : RootFinderBase
    {
        RootFinder()
        {
            m_maxdegree = maxdegree;
            m_a = m_a0;
            m_ad = m_ad0;
            m_root = m_r;
        }

    private:
        complex_t m_a0[maxdegree + 1];
        complex_t m_ad0[maxdegree + 1];
        complex_t m_r[maxdegree];
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_UTILITIES_H
#define DSPFILTERS_UTILITIES_H


namespace Dsp {

    /*
     * Utilities
     *
     * These routines are handy for manipulating buffers of samples.
     *
     */

     //------------------------------------------------------------------------------

     // Add src samples to dest, without clip or overflow checking.
    template <class Td,
        class Ts>
        void add(int samples,
            Td* dest,
            Ts const* src,
            int destSkip = 0,
            int srcSkip = 0)
    {
        if (srcSkip != 0 || destSkip != 0)
        {
            ++srcSkip;
            ++destSkip;
            while (--samples >= 0)
            {
                *dest += static_cast<Td>(*src);
                dest += destSkip;
                src += srcSkip;
            }
        }
        else
        {
            while (--samples >= 0)
                *dest++ += static_cast<Td>(*src++);
        }
    }

    // Multichannel add
    template <typename Td,
        typename Ts>
        void add(int channels,
            int samples,
            Td* const* dest,
            Ts const* const* src)
    {
        for (int i = channels; --i >= 0;)
            add(samples, dest[i], src[i]);
    }

    //--------------------------------------------------------------------------

    // Copy samples from src to dest, which may not overlap. Performs an implicit
    // type conversion if Ts and Td are different (for example, float to double).
    template <typename Td,
        typename Ts>
        void copy(int samples,
            Td* dest,
            Ts const* src,
            int destSkip = 0,
            int srcSkip = 0)
    {
        if (srcSkip != 0)
        {
            if (destSkip != 0)
            {
                ++srcSkip;
                ++destSkip;
                while (--samples >= 0)
                {
                    *dest++ = *src++;
                    dest += destSkip;
                    src += srcSkip;
                }
            }
            else
            {
                ++srcSkip;
                while (--samples >= 0)
                {
                    *dest++ = *src++;
                    src += srcSkip;
                }
            }
        }
        else if (destSkip != 0)
        {
            ++destSkip;
            while (--samples >= 0)
            {
                *dest = *src++;
                dest += destSkip;
            }
        }
        else
        {
            while (--samples >= 0)
                *dest++ = *src++;
        }
    }

    // Wrapper that uses memcpy if there is no skip and the types are the same
    template <typename Ty>
    void copy(int samples,
        Ty* dest,
        Ty const* src,
        int destSkip = 0,
        int srcSkip = 0)
    {
        if (destSkip != 0 || srcSkip != 0)
            copy<Ty, Ty>(samples, dest, src, destSkip, srcSkip);
        else
            ::memcpy(dest, src, samples * sizeof(src[0]));
    }

    // Copy a set of channels from src to dest, with implicit type conversion.
    template <typename Td,
        typename Ts>
        void copy(int channels,
            int samples,
            Td* const* dest,
            Ts const* const* src,
            int destSkip = 0,
            int srcSkip = 0)
    {
        for (int i = channels; --i >= 0;)
            copy(samples, dest[i], src[i], destSkip, srcSkip);
    }

    //--------------------------------------------------------------------------

    // Deinterleave channels. Performs implicit type conversion.
    template <typename Td, typename Ts>
    void deinterleave(int channels,
        int samples,
        Td* const* dest,
        Ts const* src)
    {
        assert(channels > 1);

        switch (channels)
        {
        case 2:
        {
            Td* l = dest[0];
            Td* r = dest[1];
            int n = (samples + 7) / 8;
            switch (samples % 8)
            {
            case 0: do
            {
                *l++ = *src++; *r++ = *src++;
            case 7:		*l++ = *src++; *r++ = *src++;
            case 6:		*l++ = *src++; *r++ = *src++;
            case 5:		*l++ = *src++; *r++ = *src++;
            case 4:		*l++ = *src++; *r++ = *src++;
            case 3:		*l++ = *src++; *r++ = *src++;
            case 2:		*l++ = *src++; *r++ = *src++;
            case 1:		*l++ = *src++; *r++ = *src++;
            } while (--n > 0);
            }
        }
        break;

        default:
        {
            for (int i = channels; --i >= 0;)
                copy(samples, dest[i], src + i, 0, channels - 1);
        }
        break;
        };
    }

    // Convenience for a stereo pair of channels
    template <typename Td,
        typename Ts>
        void deinterleave(int samples,
            Td* left,
            Td* right,
            Ts const* src)
    {
        Td* dest[2];
        dest[0] = left;
        dest[1] = right;
        deinterleave(2, samples, dest, src);
    }

    //--------------------------------------------------------------------------

    // Fade dest
    template <typename Td,
        typename Ty>
        void fade(int samples,
            Td* dest,
            Ty start = 0,
            Ty end = 1)
    {
        Ty t = start;
        Ty dt = (end - start) / samples;

        while (--samples >= 0)
        {
            *dest++ *= t;
            t += dt;
        }
    }

    // Fade dest cannels
    template <typename Td,
        typename Ty>
        void fade(int channels,
            int samples,
            Td* const* dest,
            Ty start = 0,
            Ty end = 1)
    {
        for (int i = channels; --i >= 0;)
            fade(samples, dest[i], start, end);
    }

    // Fade src into dest
    template <typename Td,
        typename Ts,
        typename Ty>
        void fade(int samples,
            Td* dest,
            Ts const* src,
            Ty start = 0,
            Ty end = 1)
    {
        Ty t = start;
        Ty dt = (end - start) / samples;

        while (--samples >= 0)
        {
            *dest = static_cast<Td>(*dest + t * (*src++ - *dest));
            dest++;
            t += dt;
        }
    }

    // Fade src channels into dest channels
    template <typename Td,
        typename Ts,
        typename Ty>
        void fade(int channels,
            int samples,
            Td* const* dest,
            Ts const* const* src,
            Ty start = 0,
            Ty end = 1)
    {
        for (int i = channels; --i >= 0;)
            fade(samples, dest[i], src[i], start, end);
    }

    //--------------------------------------------------------------------------

    // Interleave separate channels from source pointers to destination
    // (Destination requires channels*frames samples of storage). Performs
    // implicit type conversion.
    template <typename Td,
        typename Ts>
        void interleave(int channels,
            size_t samples,
            Td* dest,
            Ts const* const* src)
    {
        assert(channels > 1);

        if (samples == 0)
            return;

        switch (channels)
        {
        case 2:
        {
            const Ts* l = src[0];
            const Ts* r = src[1];

            // note that Duff's Device only works when samples>0
            int n = (samples + 7) / 8;
            switch (samples % 8)
            {
            case 0: do
            {
                *dest++ = *l++; *dest++ = *r++;
            case 7:		*dest++ = *l++; *dest++ = *r++;
            case 6:		*dest++ = *l++; *dest++ = *r++;
            case 5:		*dest++ = *l++; *dest++ = *r++;
            case 4:		*dest++ = *l++; *dest++ = *r++;
            case 3:		*dest++ = *l++; *dest++ = *r++;
            case 2:		*dest++ = *l++; *dest++ = *r++;
            case 1:		*dest++ = *l++; *dest++ = *r++;
            } while (--n > 0);
            }
        }
        break;

        default:
        {
            for (int i = channels; --i >= 0;)
                copy(samples, dest + i, src[i], channels - 1, 0);
        }
        break;
        };
    }

    //--------------------------------------------------------------------------

    // Convenience for a stereo channel pair
    template <typename Td,
        typename Ts>
        void interleave(int samples,
            Td* dest,
            Ts const* left,
            Ts const* right)
    {
        const Ts* src[2];
        src[0] = left;
        src[1] = right;
        interleave(2, samples, dest, src);
    }

    //--------------------------------------------------------------------------

    // Multiply samples by a constant, without clip or overflow checking.
    template <typename Td,
        typename Ty>
        void multiply(int samples,
            Td* dest,
            Ty factor,
            int destSkip = 0)
    {
        if (destSkip != 0)
        {
            ++destSkip;
            while (--samples >= 0)
            {
                *dest = static_cast<Td>(*dest * factor);
                dest += destSkip;
            }
        }
        else
        {
            while (--samples >= 0) {
                *dest = static_cast<Td>(*dest * factor);
                dest++;
            }
        }
    }

    // Multiply a set of channels by a constant.
    template <typename Td,
        typename Ty>
        void multiply(int channels,
            int samples,
            Td* const* dest,
            Ty factor,
            int destSkip = 0)
    {
        for (int i = channels; --i >= 0;)
            multiply(samples, dest[i], factor, destSkip);
    }

    //--------------------------------------------------------------------------

    // Copy samples from src to dest in reversed order. Performs implicit
    // type conversion. src and dest may not overlap.
    template <typename Td,
        typename Ts>
        void reverse(int samples,
            Td* dest,
            Ts const* src,
            int destSkip = 0,
            int srcSkip = 0)
    {
        src += (srcSkip + 1) * samples;

        if (srcSkip != 0 || destSkip == 0)
        {
            ++srcSkip;
            ++destSkip;
            while (--samples >= 0)
            {
                src -= srcSkip;
                *dest = *src;
                dest += destSkip;
            }
        }
        else
        {
            while (--samples >= 0)
                *dest++ = *--src;
        }
    }

    template <typename Td, typename Ts>
    void reverse(int channels, size_t frames, Td* const* dest, const Ts* const* src)
    {
        for (int i = channels; --i >= 0;)
            reverse(frames, dest[i], src[i]);
    }

    //--------------------------------------------------------------------------

    template <typename Tn>
    void to_mono(int samples, Tn* dest, Tn const* left, Tn const* right)
    {
#if 1
        while (samples-- > 0)
            *dest++ = (*left++ + *right++) * Tn(0.70710678118654752440084436210485);
#else
        while (samples-- > 0)
            *dest++ = (*left++ + *right++) * Tn(0.5);
#endif
    }

    //--------------------------------------------------------------------------

    template <typename T>
    void validate(int numChannels, int numSamples, T const* const* src)
    {
        for (int i = 0; i < numChannels; ++i)
        {
            T const* p = src[i];
            for (int j = numSamples; j > 0; --j)
            {
                T v = *p++;
                assert(v < 2 && v > -2);
            }
        }
    }

    //--------------------------------------------------------------------------

#if 0
/*
 * this stuff all depends on is_pod which is not always available
 *
 */
    namespace detail {

        template <typename Ty,
            bool isPod>
            struct zero
        {
            static void process(int samples,
                Ty* dest,
                int destSkip)
            {
                if (destSkip != 0)
                {
                    ++destSkip;
                    while (--samples >= 0)
                    {
                        *dest = Ty();
                        dest += destSkip;
                    }
                }
                else
                {
                    std::fill(dest, dest + samples, Ty());
                }
            }
        };

        template <typename Ty>
        struct zero<Ty, true>
        {
            static void process(int samples,
                Ty* dest,
                int destSkip)
            {
                if (destSkip != 0)
                    zero<Ty, false>::process(samples, dest, destSkip);
                else
                    ::memset(dest, 0, samples * sizeof(dest[0]));
            }
        };

    }

    // Fill a channel with zeros. This works even if Ty is not a basic type.
    template <typename Ty>
    void zero(int samples,
        Ty* dest,
        int destSkip = 0)
    {
        detail::zero<Ty, tr1::is_pod<Ty>::value>::process(samples, dest, destSkip);
    }

#else
// Fill a channel with zeros. This works even if Ty is not a basic type.
    template <typename Ty>
    void zero(int samples,
        Ty* dest,
        int destSkip = 0)
    {
        if (destSkip != 0)
        {
            ++destSkip;
            while (--samples >= 0)
            {
                *dest = Ty();
                dest += destSkip;
            }
        }
        else
        {
            std::fill(dest, dest + samples, Ty());
        }
    }

#endif

    // Fill a set of channels with zero.
    template <typename Ty>
    void zero(int channels,
        int samples,
        Ty* const* dest,
        int destSkip = 0)
    {
        for (int i = channels; --i >= 0;)
            zero(samples, dest[i], destSkip);
    }

    //------------------------------------------------------------------------------

    // Implementation of Brent's Method provided by
    // John D. Cook (http://www.johndcook.com/)
    // The return value of Minimize is the minimum of the function f.
    // The location where f takes its minimum is returned in the variable minLoc.
    // Notation and implementation based on Chapter 5 of Richard Brent's book
    // "Algorithms for Minimization Without Derivatives".
    //
    // Reference:
    // http://www.codeproject.com/KB/recipes/one_variable_optimize.aspx?msg=2779038

    template <class TFunction>
    double BrentMinimize
    (
        TFunction& f,	// [in] objective function to minimize
        double leftEnd,	// [in] smaller value of bracketing interval
        double rightEnd,	// [in] larger value of bracketing interval
        double epsilon,	// [in] stopping tolerance
        double& minLoc	// [out] location of minimum
    )
    {
        double d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
        static const double c = 0.5 * (3.0 - ::std::sqrt(5.0));
        static const double SQRT_DBL_EPSILON = ::std::sqrt(DBL_EPSILON);

        double& a = leftEnd;
        double& b = rightEnd;
        double& x = minLoc;

        v = w = x = a + c * (b - a);
        d = e = 0.0;
        fv = fw = fx = f(x);
        int counter = 0;
    loop:
        counter++;
        m = 0.5 * (a + b);
        tol = SQRT_DBL_EPSILON * ::fabs(x) + epsilon;
        t2 = 2.0 * tol;
        // Check stopping criteria
        if (::fabs(x - m) > t2 - 0.5 * (b - a))
        {
            p = q = r = 0.0;
            if (::fabs(e) > tol)
            {
                // fit parabola
                r = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * r;
                q = 2.0 * (q - r);
                (q > 0.0) ? p = -p : q = -q;
                r = e;
                e = d;
            }
            if (::fabs(p) < ::fabs(0.5 * q * r) && p < q * (a - x) && p < q * (b - x))
            {
                // A parabolic interpolation step
                d = p / q;
                u = x + d;
                // f must not be evaluated too close to a or b
                if (u - a < t2 || b - u < t2)
                    d = (x < m) ? tol : -tol;
            }
            else
            {
                // A golden section step
                e = (x < m) ? b : a;
                e -= x;
                d = c * e;
            }
            // f must not be evaluated too close to x
            if (::fabs(d) >= tol)
                u = x + d;
            else if (d > 0.0)
                u = x + tol;
            else
                u = x - tol;
            fu = f(u);
            // Update a, b, v, w, and x
            if (fu <= fx)
            {
                (u < x) ? b = x : a = x;
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            }
            else
            {
                (u < x) ? a = u : b = u;
                if (fu <= fw || w == x)
                {
                    v = w;
                    fv = fw;
                    w = u;
                    fw = fu;
                }
                else if (fu <= fv || v == x || v == w)
                {
                    v = u;
                    fv = fu;
                }
            }
            goto loop;  // Yes, the dreaded goto statement. But the code
            // here is faithful to Brent's orginal pseudocode.
        }
        return  fx;
    }

    //------------------------------------------------------------------------------

    // Tracks the peaks in the signal stream using the attack and release parameters
    template <int Channels = 2, typename Value = float>
    class EnvelopeFollower
    {
    public:
        EnvelopeFollower()
        {
            for (int i = 0; i < Channels; i++)
                m_env[i] = 0;
        }

        Value operator[] (int channel) const
        {
            return m_env[channel];
        }

        void Setup(int sampleRate, double attackMs, double releaseMs)
        {
            m_a = pow(0.01, 1.0 / (attackMs * sampleRate * 0.001));
            m_r = pow(0.01, 1.0 / (releaseMs * sampleRate * 0.001));
        }

        void Process(size_t samples, const Value** src)
        {
            for (int i = 0; i < Channels; i++)
            {
                const Value* cur = src[i];

                double e = m_env[i];
                for (int n = samples; n; n--)
                {
                    double v = std::abs(*cur++);
                    if (v > e)
                        e = m_a * (e - v) + v;
                    else
                        e = m_r * (e - v) + v;
                }
                m_env[i] = e;
            }
        }

        double m_env[Channels];

    protected:
        double m_a;
        double m_r;
    };

    //------------------------------------------------------------------------------

    // Helpful for discovering discontinuities in buffers
    template <int Channels = 2, typename Value = float>
    class SlopeDetector
    {
    public:
        SlopeDetector() : m_firstTime(true)
        {
            for (int i = 0; i < Channels; ++i)
                m_slope[i] = 0;
        }

        Value getSlope(int channel) const
        {
            return m_slope[channel];
        }

        void process(size_t numSamples, const Value** input)
        {
            for (int i = 0; i < Channels; ++i)
            {
                const Value* src = input[i];
                int n = numSamples;

                if (m_firstTime)
                {
                    m_prev[i] = *src++;
                    --n;
                }

                while (n > 0)
                {
                    n--;
                    Value cur = *src++;
                    Value diff = std::abs(cur - m_prev[i]);
                    m_slope[i] = std::max(diff, m_slope[i]);
                    m_prev[i] = cur;
                }
            }

            m_firstTime = false;
        }

    private:
        bool m_firstTime;
        Value m_slope[Channels];
        Value m_prev[Channels];
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_CASCADE_H
#define DSPFILTERS_CASCADE_H


namespace Dsp {

    /*
     * Holds coefficients for a cascade of second order sections.
     *
     */

     // Factored implementation to reduce template instantiations
    class Cascade
    {
    public:
        template <class StateType>
        class StateBase : private DenormalPrevention
        {
        public:
            template <typename Sample>
            inline Sample process(const Sample in, const Cascade& c)
            {
                double out = in;
                StateType* state = m_stateArray;
                Biquad const* stage = c.m_stageArray;
                const double vsa = ac();
                int i = c.m_numStages - 1;
                out = (state++)->process1(out, *stage++, vsa);
                for (; --i >= 0;)
                    out = (state++)->process1(out, *stage++, 0);
                //for (int i = c.m_numStages; --i >= 0; ++state, ++stage)
                //  out = state->process1 (out, *stage, vsa);
                return static_cast<Sample> (out);
            }

        protected:
            StateBase(StateType* stateArray)
                : m_stateArray(stateArray)
            {
            }

        protected:
            StateType* m_stateArray;
        };

        struct Stage : Biquad
        {
        };

        struct Storage
        {
            Storage(int maxStages_, Stage* stageArray_)
                : maxStages(maxStages_)
                , stageArray(stageArray_)
            {
            }

            int maxStages;
            Stage* stageArray;
        };

        int getNumStages() const
        {
            return m_numStages;
        }

        const Stage& operator[] (int index)
        {
            assert(index >= 0 && index <= m_numStages);
            return m_stageArray[index];
        }

    public:
        // Calculate filter response at the given normalized frequency.
        complex_t response(double normalizedFrequency) const;

        std::vector<PoleZeroPair> getPoleZeros() const;

        // Process a block of samples in the given form
        template <class StateType, typename Sample>
        void process(int numSamples, Sample* dest, StateType& state) const
        {
            while (--numSamples >= 0) {
                *dest = state.process(*dest, *this);
                dest++;
            }
        }

    protected:
        Cascade();

        void setCascadeStorage(const Storage& storage);

        void applyScale(double scale);
        void setLayout(const LayoutBase& proto);

    private:
        int m_numStages;
        int m_maxStages;
        Stage* m_stageArray;
    };

    //------------------------------------------------------------------------------

    // Storage for Cascade
    template <int MaxStages>
    class CascadeStages
    {
    public:
        template <class StateType>
        class State : public Cascade::StateBase <StateType>
        {
        public:
            State() : Cascade::StateBase <StateType>(m_states)
            {
                Cascade::StateBase <StateType>::m_stateArray = m_states;
                reset();
            }

            void reset()
            {
                StateType* state = m_states;
                for (int i = MaxStages; --i >= 0; ++state)
                    state->reset();
            }

        private:
            StateType m_states[MaxStages];
        };

        /*@Internal*/
        Cascade::Storage getCascadeStorage()
        {
            return Cascade::Storage(MaxStages, m_stages);
        }

    private:
        Cascade::Stage m_stages[MaxStages];
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_POLEFILTER_H
#define DSPFILTERS_POLEFILTER_H

namespace Dsp {

    /*
     * Base for filters designed via algorithmic placement of poles and zeros.
     *
     * Typically, the filter is first designed as a half-band low pass or
     * low shelf analog filter (s-plane). Then, using a transformation such
     * as the ones from Constantinides, the poles and zeros of the analog filter
     * are calculated in the z-plane.
     *
     */

     // Factored implementations to reduce template instantiations

    class PoleFilterBase2 : public Cascade
    {
    public:
        // This gets the poles/zeros directly from the digital
        // prototype. It is used to double check the correctness
        // of the recovery of pole/zeros from biquad coefficients.
        //
        // It can also be used to accelerate the interpolation
        // of pole/zeros for parameter modulation, since a pole
        // filter already has them calculated

#if 1
  // Commenting this out will pass the call to the Cascade,
  // which tries to compute the poles and zeros from the biquad
  // coefficients.
        std::vector<PoleZeroPair> getPoleZeros() const
        {
            std::vector<PoleZeroPair> vpz;
            const int pairs = (m_digitalProto.getNumPoles() + 1) / 2;
            for (int i = 0; i < pairs; ++i)
                vpz.push_back(m_digitalProto[i]);
            return vpz;
        }
#endif

    protected:
        LayoutBase m_digitalProto;
    };

    // Serves a container to hold the analog prototype
    // and the digital pole/zero layout.
    template <class AnalogPrototype>
    class PoleFilterBase : public PoleFilterBase2
    {
    protected:
        void setPrototypeStorage(const LayoutBase& analogStorage,
            const LayoutBase& digitalStorage)
        {
            m_analogProto.setStorage(analogStorage);
            m_digitalProto = digitalStorage;
        }

    protected:
        AnalogPrototype m_analogProto;
    };

    //------------------------------------------------------------------------------

    // Storage for pole filters
    template <class BaseClass,
        int MaxAnalogPoles,
        int MaxDigitalPoles = MaxAnalogPoles>
        struct PoleFilter : BaseClass
        , CascadeStages <(MaxDigitalPoles + 1) / 2>
    {
        PoleFilter()
        {
            // This glues together the factored base classes
            // with the templatized storage classes.
            BaseClass::setCascadeStorage(this->getCascadeStorage());
            BaseClass::setPrototypeStorage(m_analogStorage, m_digitalStorage);
        }

    private:
        Layout <MaxAnalogPoles> m_analogStorage;
        Layout <MaxDigitalPoles> m_digitalStorage;
    };

    //------------------------------------------------------------------------------

    /*
     * s-plane to z-plane transforms
     *
     * For pole filters, an analog prototype is created via placement of
     * poles and zeros in the s-plane. The analog prototype is either
     * a halfband low pass or a halfband low shelf. The poles, zeros,
     * and normalization parameters are transformed into the z-plane
     * using variants of the bilinear transformation.
     *
     */

     // low pass to low pass 
    class LowPassTransform
    {
    public:
        LowPassTransform(double fc,
            LayoutBase& digital,
            LayoutBase const& analog);

    private:
        complex_t transform(complex_t c);

        double f;
    };

    //------------------------------------------------------------------------------

    // low pass to high pass
    class HighPassTransform
    {
    public:
        HighPassTransform(double fc,
            LayoutBase& digital,
            LayoutBase const& analog);

    private:
        complex_t transform(complex_t c);

        double f;
    };

    //------------------------------------------------------------------------------

    // low pass to band pass transform
    class BandPassTransform
    {

    public:
        BandPassTransform(double fc,
            double fw,
            LayoutBase& digital,
            LayoutBase const& analog);

    private:
        ComplexPair transform(complex_t c);

        double wc;
        double wc2;
        double a;
        double b;
        double a2;
        double b2;
        double ab;
        double ab_2;
    };

    //------------------------------------------------------------------------------

    // low pass to band stop transform
    class BandStopTransform
    {
    public:
        BandStopTransform(double fc,
            double fw,
            LayoutBase& digital,
            LayoutBase const& analog);

    private:
        ComplexPair transform(complex_t c);

        double wc;
        double wc2;
        double a;
        double b;
        double a2;
        double b2;
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_DESIGN_H
#define DSPFILTERS_DESIGN_H

namespace Dsp {

    struct DesignBase
    {
        // Sampling rate is the first param for every Design filter
        static const ParamInfo getParamInfo_0()
        {
            return ParamInfo::defaultSampleRateParam();
        }

        // These should never get called
        static const ParamInfo getParamInfo_1() { return ParamInfo(); }
        static const ParamInfo getParamInfo_2() { return ParamInfo(); }
        static const ParamInfo getParamInfo_3() { return ParamInfo(); }
        static const ParamInfo getParamInfo_4() { return ParamInfo(); }
        static const ParamInfo getParamInfo_5() { return ParamInfo(); }
        static const ParamInfo getParamInfo_6() { return ParamInfo(); }
        static const ParamInfo getParamInfo_7() { return ParamInfo(); }
    };

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_CHEBYSHEVI_H
#define DSPFILTERS_CHEBYSHEVI_H


namespace Dsp {

    /*
     * Filters with Chebyshev response characteristics
     *
     */

    namespace ChebyshevI {

        // Half-band analog prototypes (s-plane)

        class AnalogLowPass : public LayoutBase
        {
        public:
            AnalogLowPass();

            void design(const int numPoles,
                double rippleDb);

        private:
            int m_numPoles;
            double m_rippleDb;
        };

        //------------------------------------------------------------------------------

        class AnalogLowShelf : public LayoutBase
        {
        public:
            AnalogLowShelf();

            void design(int numPoles,
                double gainDb,
                double rippleDb);

        private:
            int m_numPoles;
            double m_rippleDb;
            double m_gainDb;
        };

        //------------------------------------------------------------------------------

        // Factored implementations to reduce template instantiations

        struct LowPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double rippleDb);
        };

        struct HighPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double rippleDb);
        };

        struct BandPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double rippleDb);
        };

        struct BandStopBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double rippleDb);
        };

        struct LowShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb,
                double rippleDb);
        };

        struct HighShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb,
                double rippleDb);
        };

        struct BandShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double gainDb,
                double rippleDb);
        };

        //------------------------------------------------------------------------------

        //
        // Raw filters
        //

        template <int MaxOrder>
        struct LowPass : PoleFilter <LowPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighPass : PoleFilter <HighPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct LowShelf : PoleFilter <LowShelfBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighShelf : PoleFilter <HighShelfBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandShelf : PoleFilter <BandShelfBase, MaxOrder, MaxOrder * 2>
        {
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultRippleDbParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultRippleDbParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            struct TypeIIIBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultRippleDbParam();
                }
            };

            template <class FilterClass>
            struct TypeIII : TypeIIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            struct TypeIVBase : DesignBase
            {
                enum
                {
                    NumParams = 6
                };

                static int getNumParams()
                {
                    return 6;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_5()
                {
                    return ParamInfo::defaultRippleDbParam();
                }
            };

            template <class FilterClass>
            struct TypeIV : TypeIVBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4], params[5]);
                }
            };

            // Factored kind and name

            struct LowPassDescription
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "Chebyshev I Low Pass"; }
            };

            struct HighPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Chebyshev I High Pass"; }
            };

            struct BandPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Chebyshev I Band Pass"; }
            };

            struct BandStopDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Chebyshev I Band Stop"; }
            };

            struct LowShelfDescription
            {
                static Kind getKind() { return kindLowShelf; }
                static const char* getName() { return "Chebyshev I Low Shelf"; }
            };

            struct HighShelfDescription
            {
                static Kind getKind() { return kindHighShelf; }
                static const char* getName() { return "Chebyshev I High Shelf"; }
            };

            struct BandShelfDescription
            {
                static Kind getKind() { return kindBandShelf; }
                static const char* getName() { return "Chebyshev I Band Shelf"; }
            };

            // This glues on the Order parameter
            template <int MaxOrder,
                template <class> class TypeClass,
                template <int> class FilterClass>
            struct OrderBase : TypeClass <FilterClass <MaxOrder> >
            {
                const ParamInfo getParamInfo_1() const
                {
                    return ParamInfo(idOrder, "Order", "Order",
                        1, MaxOrder, 2,
                        &ParamInfo::Int_toControlValue,
                        &ParamInfo::Int_toNativeValue,
                        &ParamInfo::Int_toString);

                }
            };

            //------------------------------------------------------------------------------

            //
            // Design filters
            //

            template <int MaxOrder>
            struct LowPass : OrderBase <MaxOrder, TypeI, ChebyshevI::LowPass>,
                LowPassDescription
            {
            };

            template <int MaxOrder>
            struct HighPass : OrderBase <MaxOrder, TypeI, ChebyshevI::HighPass>,
                HighPassDescription
            {
            };

            template <int MaxOrder>
            struct BandPass : OrderBase <MaxOrder, TypeII, ChebyshevI::BandPass>,
                BandPassDescription
            {
            };

            template <int MaxOrder>
            struct BandStop : OrderBase <MaxOrder, TypeII, ChebyshevI::BandStop>,
                BandStopDescription
            {
            };

            template <int MaxOrder>
            struct LowShelf : OrderBase <MaxOrder, TypeIII, ChebyshevI::LowShelf>,
                LowShelfDescription
            {
            };

            template <int MaxOrder>
            struct HighShelf : OrderBase <MaxOrder, TypeIII, ChebyshevI::HighShelf>,
                HighShelfDescription
            {
            };

            template <int MaxOrder>
            struct BandShelf : OrderBase <MaxOrder, TypeIV, ChebyshevI::BandShelf>,
                BandShelfDescription
            {
            };

        }

    }

}

#endif

/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_CHEBYSHEVII_H
#define DSPFILTERS_CHEBYSHEVII_H


namespace Dsp {

    /*
     * Filters with Inverse Chebyshev response characteristics
     *
     */

    namespace ChebyshevII {

        // Half-band analog prototypes (s-plane)

        class AnalogLowPass : public LayoutBase
        {
        public:
            AnalogLowPass();

            void design(const int numPoles,
                double stopBandDb);

        private:
            int m_numPoles;
            double m_stopBandDb;
        };

        //------------------------------------------------------------------------------

        class AnalogLowShelf : public LayoutBase
        {
        public:
            AnalogLowShelf();

            void design(int numPoles,
                double gainDb,
                double stopBandDb);

        private:
            int m_numPoles;
            double m_stopBandDb;
            double m_gainDb;
        };

        //------------------------------------------------------------------------------

        // Factored implementations to reduce template instantiations

        struct LowPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double stopBandDb);
        };

        struct HighPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double stopBandDb);
        };

        struct BandPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double stopBandDb);
        };

        struct BandStopBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double stopBandDb);
        };

        struct LowShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb,
                double stopBandDb);
        };

        struct HighShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb,
                double stopBandDb);
        };

        struct BandShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double gainDb,
                double stopBandDb);
        };

        //------------------------------------------------------------------------------

        //
        // Raw filters
        //

        template <int MaxOrder>
        struct LowPass : PoleFilter <LowPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighPass : PoleFilter <HighPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct LowShelf : PoleFilter <LowShelfBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighShelf : PoleFilter <HighShelfBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandShelf : PoleFilter <BandShelfBase, MaxOrder, MaxOrder * 2>
        {
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultStopDbParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultStopDbParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            struct TypeIIIBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultStopDbParam();
                }
            };

            template <class FilterClass>
            struct TypeIII : TypeIIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            struct TypeIVBase : DesignBase
            {
                enum
                {
                    NumParams = 6
                };

                static int getNumParams()
                {
                    return 6;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_5()
                {
                    return ParamInfo::defaultStopDbParam();
                }
            };

            template <class FilterClass>
            struct TypeIV : TypeIVBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4], params[5]);
                }
            };

            // Factored kind and name

            struct LowPassDescription
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "Chebyshev II Low Pass"; }
            };

            struct HighPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Chebyshev II High Pass"; }
            };

            struct BandPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Chebyshev II Band Pass"; }
            };

            struct BandStopDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Chebyshev II Band Stop"; }
            };

            struct LowShelfDescription
            {
                static Kind getKind() { return kindLowShelf; }
                static const char* getName() { return "Chebyshev II Low Shelf"; }
            };

            struct HighShelfDescription
            {
                static Kind getKind() { return kindHighShelf; }
                static const char* getName() { return "Chebyshev II High Shelf"; }
            };

            struct BandShelfDescription
            {
                static Kind getKind() { return kindBandShelf; }
                static const char* getName() { return "Chebyshev II Band Shelf"; }
            };

            // This glues on the Order parameter
            template <int MaxOrder,
                template <class> class TypeClass,
                template <int> class FilterClass>
            struct OrderBase : TypeClass <FilterClass <MaxOrder> >
            {
                const ParamInfo getParamInfo_1() const
                {
                    return ParamInfo(idOrder, "Order", "Order",
                        1, MaxOrder, 2,
                        &ParamInfo::Int_toControlValue,
                        &ParamInfo::Int_toNativeValue,
                        &ParamInfo::Int_toString);

                }
            };

            //------------------------------------------------------------------------------

            //
            // Design Filters
            //

            template <int MaxOrder>
            struct LowPass : OrderBase <MaxOrder, TypeI, ChebyshevII::LowPass>,
                LowPassDescription
            {
            };

            template <int MaxOrder>
            struct HighPass : OrderBase <MaxOrder, TypeI, ChebyshevII::HighPass>,
                HighPassDescription
            {
            };

            template <int MaxOrder>
            struct BandPass : OrderBase <MaxOrder, TypeII, ChebyshevII::BandPass>,
                BandPassDescription
            {
            };

            template <int MaxOrder>
            struct BandStop : OrderBase <MaxOrder, TypeII, ChebyshevII::BandStop>,
                BandStopDescription
            {
            };

            template <int MaxOrder>
            struct LowShelf : OrderBase <MaxOrder, TypeIII, ChebyshevII::LowShelf>,
                LowShelfDescription
            {
            };

            template <int MaxOrder>
            struct HighShelf : OrderBase <MaxOrder, TypeIII, ChebyshevII::HighShelf>,
                HighShelfDescription
            {
            };

            template <int MaxOrder>
            struct BandShelf : OrderBase <MaxOrder, TypeIV, ChebyshevII::BandShelf>,
                BandShelfDescription
            {
            };


        }

    }

}

#endif

/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_BUTTERWORTH_H
#define DSPFILTERS_BUTTERWORTH_H


namespace Dsp {

    /*
     * Filters with Butterworth response characteristics
     *
     */

    namespace Butterworth {

        // Half-band analog prototypes (s-plane)

        class AnalogLowPass : public LayoutBase
        {
        public:
            AnalogLowPass();

            void design(const int numPoles);

        private:
            int m_numPoles;
        };

        //------------------------------------------------------------------------------

        class AnalogLowShelf : public LayoutBase
        {
        public:
            AnalogLowShelf();

            void design(int numPoles, double gainDb);

        private:
            int m_numPoles;
            double m_gainDb;
        };

        //------------------------------------------------------------------------------

        // Factored implementations to reduce template instantiations

        struct LowPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency);
        };

        struct HighPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency);
        };

        struct BandPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency);
        };

        struct BandStopBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency);
        };

        struct LowShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb);
        };

        struct HighShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb);
        };

        struct BandShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double gainDb);
        };

        //------------------------------------------------------------------------------

        //
        // Raw filters
        //

        template <int MaxOrder>
        struct LowPass : PoleFilter <LowPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighPass : PoleFilter <HighPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct LowShelf : PoleFilter <LowShelfBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighShelf : PoleFilter <HighShelfBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandShelf : PoleFilter <BandShelfBase, MaxOrder, MaxOrder * 2>
        {
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 3
                };

                static int getNumParams()
                {
                    return 3;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
                }
            };

            struct TypeIIIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultGainParam();
                }
            };

            template <class FilterClass>
            struct TypeIII : TypeIIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]),
                        params[0],
                        params[2],
                        params[3]);
                }
            };

            struct TypeIVBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultGainParam();
                }
            };

            template <class FilterClass>
            struct TypeIV : TypeIVBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            // Factored kind and name

            struct LowPassDescription
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "Butterworth Low Pass"; }
            };

            struct HighPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Butterworth High Pass"; }
            };

            struct BandPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Butterworth Band Pass"; }
            };

            struct BandStopDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Butterworth Band Stop"; }
            };

            struct LowShelfDescription
            {
                static Kind getKind() { return kindLowShelf; }
                static const char* getName() { return "Butterworth Low Shelf"; }
            };

            struct HighShelfDescription
            {
                static Kind getKind() { return kindHighShelf; }
                static const char* getName() { return "Butterworth High Shelf"; }
            };

            struct BandShelfDescription
            {
                static Kind getKind() { return kindBandShelf; }
                static const char* getName() { return "Butterworth Band Shelf"; }
            };

            // This glues on the Order parameter
            template <int MaxOrder,
                template <class> class TypeClass,
                template <int> class FilterClass>
            struct OrderBase : TypeClass <FilterClass <MaxOrder> >
            {
                const ParamInfo getParamInfo_1() const
                {
                    return ParamInfo(idOrder, "Order", "Order",
                        1, MaxOrder, 2,
                        &ParamInfo::Int_toControlValue,
                        &ParamInfo::Int_toNativeValue,
                        &ParamInfo::Int_toString);

                }
            };

            //------------------------------------------------------------------------------

            //
            // Design filters
            //

            template <int MaxOrder>
            struct LowPass : OrderBase <MaxOrder, TypeI, Butterworth::LowPass>,
                LowPassDescription
            {
            };

            template <int MaxOrder>
            struct HighPass : OrderBase <MaxOrder, TypeI, Butterworth::HighPass>,
                HighPassDescription
            {
            };

            template <int MaxOrder>
            struct BandPass : OrderBase <MaxOrder, TypeII, Butterworth::BandPass>,
                BandPassDescription
            {
            };

            template <int MaxOrder>
            struct BandStop : OrderBase <MaxOrder, TypeII, Butterworth::BandStop>,
                BandStopDescription
            {
            };

            template <int MaxOrder>
            struct LowShelf : OrderBase <MaxOrder, TypeIII, Butterworth::LowShelf>,
                LowShelfDescription
            {
            };

            template <int MaxOrder>
            struct HighShelf : OrderBase <MaxOrder, TypeIII, Butterworth::HighShelf>,
                HighShelfDescription
            {
            };

            template <int MaxOrder>
            struct BandShelf : OrderBase <MaxOrder, TypeIV, Butterworth::BandShelf>,
                BandShelfDescription
            {
            };

        }

    }

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_BESSEL_H
#define DSPFILTERS_BESSEL_H


namespace Dsp {

    /*
     * Filters with Bessel response characteristics
     *
     */

    namespace Bessel {

        // A Workspace is necessary to find roots

        struct WorkspaceBase
        {
            WorkspaceBase(RootFinderBase* rootsBase)
                : roots(*rootsBase)
            {
            }

            RootFinderBase& roots;

        private:
            WorkspaceBase(WorkspaceBase&);
            WorkspaceBase& operator= (WorkspaceBase&);
        };

        template <int MaxOrder>
        struct Workspace : WorkspaceBase
        {
            Workspace()
                : WorkspaceBase(&m_roots)
            {
            }

        private:
            RootFinder <MaxOrder> m_roots;
        };

        //------------------------------------------------------------------------------

        // Half-band analog prototypes (s-plane)

        class AnalogLowPass : public LayoutBase
        {
        public:
            AnalogLowPass();

            void design(const int numPoles,
                WorkspaceBase* w);

        private:
            int m_numPoles;
        };

        //------------------------------------------------------------------------------

        class AnalogLowShelf : public LayoutBase
        {
        public:
            AnalogLowShelf();

            void design(int numPoles,
                double gainDb,
                WorkspaceBase* w);

        private:
            int m_numPoles;
            double m_gainDb;
        };

        //------------------------------------------------------------------------------

        // Factored implementations to reduce template instantiations

        struct LowPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                WorkspaceBase* w);
        };

        struct HighPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                WorkspaceBase* w);
        };

        struct BandPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                WorkspaceBase* w);
        };

        struct BandStopBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                WorkspaceBase* w);
        };

        struct LowShelfBase : PoleFilterBase <AnalogLowShelf>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb,
                WorkspaceBase* w);
        };

        //------------------------------------------------------------------------------

        //
        // Raw filters
        //

        template <int MaxOrder>
        struct LowPass : PoleFilter <LowPassBase, MaxOrder>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency)
            {
                Workspace <MaxOrder> w;
                LowPassBase::setup(order,
                    sampleRate,
                    cutoffFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct HighPass : PoleFilter <HighPassBase, MaxOrder>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency)
            {
                Workspace <MaxOrder> w;
                HighPassBase::setup(order,
                    sampleRate,
                    cutoffFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder * 2>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency)
            {
                Workspace <MaxOrder> w;
                BandPassBase::setup(order,
                    sampleRate,
                    centerFrequency,
                    widthFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder * 2>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency)
            {
                Workspace <MaxOrder> w;
                BandStopBase::setup(order,
                    sampleRate,
                    centerFrequency,
                    widthFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct LowShelf : PoleFilter <LowShelfBase, MaxOrder, MaxOrder * 2>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double gainDb)
            {
                Workspace <MaxOrder> w;
                LowShelfBase::setup(order,
                    sampleRate,
                    cutoffFrequency,
                    gainDb,
                    &w);
            }
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 3
                };

                static int getNumParams()
                {
                    return 3;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
                }
            };

            struct TypeIIIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultGainParam();
                }
            };

            template <class FilterClass>
            struct TypeIII : TypeIIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]),
                        params[0],
                        params[2],
                        params[3]);
                }
            };

            struct TypeIVBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultGainParam();
                }
            };

            template <class FilterClass>
            struct TypeIV : TypeIVBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            // Factored kind and name

            struct LowPassDescription
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "Bessel Low Pass"; }
            };

            struct HighPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Bessel High Pass"; }
            };

            struct BandPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Bessel Band Pass"; }
            };

            struct BandStopDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Bessel Band Stop"; }
            };

            struct LowShelfDescription
            {
                static Kind getKind() { return kindLowShelf; }
                static const char* getName() { return "Bessel Low Shelf"; }
            };

            // This glues on the Order parameter
            template <int MaxOrder,
                template <class> class TypeClass,
                template <int> class FilterClass>
            struct OrderBase : TypeClass <FilterClass <MaxOrder> >
            {
                const ParamInfo getParamInfo_1() const
                {
                    return ParamInfo(idOrder, "Order", "Order",
                        1, MaxOrder, 2,
                        &ParamInfo::Int_toControlValue,
                        &ParamInfo::Int_toNativeValue,
                        &ParamInfo::Int_toString);

                }
            };

            //------------------------------------------------------------------------------

            //
            // Gui-friendly Design layer
            //

            template <int MaxOrder>
            struct LowPass : OrderBase <MaxOrder, TypeI, Bessel::LowPass>,
                LowPassDescription
            {
            };

            template <int MaxOrder>
            struct HighPass : OrderBase <MaxOrder, TypeI, Bessel::HighPass>,
                HighPassDescription
            {
            };

            template <int MaxOrder>
            struct BandPass : OrderBase <MaxOrder, TypeII, Bessel::BandPass>,
                BandPassDescription
            {
            };

            template <int MaxOrder>
            struct BandStop : OrderBase <MaxOrder, TypeII, Bessel::BandStop>,
                BandStopDescription
            {
            };

            /*
             * NOT IMPLEMENTED
             *
             */
            template <int MaxOrder>
            struct LowShelf : OrderBase <MaxOrder, TypeIII, Bessel::LowShelf>,
                LowShelfDescription
            {
            };

        }

    }

}

#endif

/* This is a test of svn:external */

/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_CUSTOM_H
#define DSPFILTERS_CUSTOM_H

namespace Dsp {

    /*
     * Single pole and Biquad with parameters allowing
     * for directly setting the poles and zeros
     *
     */

    namespace Custom {

        //
        // Raw filters
        //

        struct OnePole : Biquad
        {
            void setup(double scale,
                double pole,
                double zero);
        };

        struct TwoPole : Biquad
        {
            void setup(double scale,
                double poleRho,
                double poleTheta,
                double zeroRho,
                double zeroTheta);
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct OnePole : DesignBase, Custom::OnePole
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_1()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultPoleRealParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultZeroRealParam();
                }

                static Kind getKind() { return kindOther; }
                static const char* getName() { return "Custom One-Pole"; }

                void setParams(const Params& params)
                {
                    setup(pow(10., params[1] / 20),
                        params[2],
                        params[3]);
                }
            };

            struct TwoPole : DesignBase, Custom::TwoPole
            {
                enum
                {
                    NumParams = 6
                };

                static int getNumParams()
                {
                    return 6;
                }

                static const ParamInfo getParamInfo_1()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultPoleRhoParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultPoleThetaParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultZeroRhoParam();
                }

                static const ParamInfo getParamInfo_5()
                {
                    return ParamInfo::defaultZeroThetaParam();
                }


                static Kind getKind() { return kindOther; }
                static const char* getName() { return "Custom Two-Pole"; }

                void setParams(const Params& params)
                {
                    setup(pow(10., params[1] / 20),
                        params[2],
                        params[3],
                        params[4],
                        params[5]);
                }
            };

        }

    }

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_ELLIPTIC_H
#define DSPFILTERS_ELLIPTIC_H

namespace Dsp {

    /*
     * Filters with Elliptic response characteristics
     *
     */

    namespace Elliptic {

        // Solves for Jacobi elliptics
        class Solver
        {
        public:
            static double ellipticK(double k);
        };

        // Half-band analog prototype (s-plane)

        class AnalogLowPass : public LayoutBase
        {
        public:
            AnalogLowPass();

            void design(const int numPoles,
                double rippleDb,
                double rolloff);

        private:
            void prodpoly(int sn);
            void calcfz2(int i);
            void calcfz();
            void calcqz();
            double findfact(int t);
            double calcsn(double u);

#if 0
            template<int n>
            struct CalcArray
            {
                double& operator[](size_t index)
                {
                    assert(index < n);
                    return m_a[index];
                }
            private:
                double m_a[n];
            };
#else
#endif

            double m_p0;
            double m_q;
            double m_K;
            double m_Kprime;
            double m_e;
            int m_nin;
            int m_m;
            int m_n2;
            int m_em;
            double m_zeros[100];
            double m_c1[100];
            double m_b1[100];
            double m_a1[100];
            double m_d1[100];
            double m_q1[100];
            double m_z1[100];
            double m_f1[100];
            double m_s1[100];
            double m_p[100];
            double m_zw1[100];
            double m_zf1[100];
            double m_zq1[100];
            double m_rootR[100];
            double m_rootI[100];

            int m_numPoles;
            double m_rippleDb;
            double m_rolloff;
        };

        //------------------------------------------------------------------------------

        // Factored implementations to reduce template instantiations

        struct LowPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double rippleDb,
                double rolloff);
        };

        struct HighPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                double rippleDb,
                double rolloff);
        };

        struct BandPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double rippleDb,
                double rolloff);
        };

        struct BandStopBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                double rippleDb,
                double rolloff);
        };

        //------------------------------------------------------------------------------

        //
        // Raw filters
        //

        template <int MaxOrder>
        struct LowPass : PoleFilter <LowPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct HighPass : PoleFilter <HighPassBase, MaxOrder>
        {
        };

        template <int MaxOrder>
        struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder * 2>
        {
        };

        template <int MaxOrder>
        struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder * 2>
        {
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 5
                };

                static int getNumParams()
                {
                    return 5;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultRippleDbParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultRolloffParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 6
                };

                static int getNumParams()
                {
                    return 6;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }

                static const ParamInfo getParamInfo_4()
                {
                    return ParamInfo::defaultRippleDbParam();
                }

                static const ParamInfo getParamInfo_5()
                {
                    return ParamInfo::defaultRolloffParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3], params[4], params[5]);
                }
            };

            // Factored kind and name

            struct LowPassDescription
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "Elliptic Low Pass"; }
            };

            struct HighPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Elliptic High Pass"; }
            };

            struct BandPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Elliptic Band Pass"; }
            };

            struct BandStopDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Elliptic Band Stop"; }
            };

            // This glues on the Order parameter
            template <int MaxOrder,
                template <class> class TypeClass,
                template <int> class FilterClass>
            struct OrderBase : TypeClass <FilterClass <MaxOrder> >
            {
                const ParamInfo getParamInfo_1() const
                {
                    return ParamInfo(idOrder, "Order", "Order",
                        1, MaxOrder, 2,
                        &ParamInfo::Int_toControlValue,
                        &ParamInfo::Int_toNativeValue,
                        &ParamInfo::Int_toString);

                }
            };
            //------------------------------------------------------------------------------

            //
            // Design filters
            //

            template <int MaxOrder>
            struct LowPass : OrderBase <MaxOrder, TypeI, Elliptic::LowPass>,
                LowPassDescription
            {
            };

            template <int MaxOrder>
            struct HighPass : OrderBase <MaxOrder, TypeI, Elliptic::HighPass>,
                HighPassDescription
            {
            };

            template <int MaxOrder>
            struct BandPass : OrderBase <MaxOrder, TypeII, Elliptic::BandPass>,
                BandPassDescription
            {
            };

            template <int MaxOrder>
            struct BandStop : OrderBase <MaxOrder, TypeII, Elliptic::BandStop>,
                BandStopDescription
            {
            };

        }

    }

}

#endif

/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_LEGENDRE_H
#define DSPFILTERS_LEGENDRE_H

namespace Dsp {

    /*
     * Filters with Legendre / "Optimum-L" response characteristics
     *
     */

    namespace Legendre {

        // Numerical computation of Legendre "Optimum-L" polynomials

        class PolynomialFinderBase
        {
        public:
            void solve(int n);

            double* coef()
            {
                return m_w;
            }

        private:
            void legendre(double* p, int n);

        protected:
            int m_maxN;
            double* m_w;
            double* m_a;
            double* m_p;
            double* m_s;
            double* m_v;
            double* m_aa;
            double* m_bb;
        };

        template <int maxN>
        class PolynomialFinder : public PolynomialFinderBase
        {
        public:
            PolynomialFinder()
            {
                m_maxN = maxN;
                m_w = m_ws;
                m_a = m_as;
                m_p = m_ps;
                m_s = m_ss;
                m_v = m_vs;
                m_aa = m_aas;
                m_bb = m_bbs;
            }

            void solve(int n)
            {
                assert(n <= maxN);
                PolynomialFinderBase::solve(n);
            }

        private:
            double m_ws[2 * maxN + 1];
            double m_as[maxN + 1];
            double m_ps[2 * maxN + 1];
            double m_ss[2 * maxN + 1];
            double m_vs[2 * maxN + 4];
            double m_aas[maxN + 1];
            double m_bbs[maxN + 1];
        };

        //------------------------------------------------------------------------------

        // A Workspace is necessary to construct the polynomial and find its roots

        struct WorkspaceBase
        {
            WorkspaceBase(PolynomialFinderBase* polyBase,
                RootFinderBase* rootsBase)
                : poly(*polyBase)
                , roots(*rootsBase)
            {
            }

            PolynomialFinderBase& poly;
            RootFinderBase& roots;

        private:
            WorkspaceBase(WorkspaceBase&);
            WorkspaceBase& operator= (WorkspaceBase&);
        };

        template <int MaxOrder>
        struct Workspace : WorkspaceBase
        {
            Workspace()
                : WorkspaceBase(&m_poly, &m_roots)
            {
            }

        private:
            PolynomialFinder <MaxOrder> m_poly;
            RootFinder <MaxOrder * 2> m_roots;
        };

        //------------------------------------------------------------------------------

        // Half-band analog prototypes (s-plane)

        class AnalogLowPass : public LayoutBase
        {
        public:
            AnalogLowPass();

            void design(const int numPoles, WorkspaceBase* w);

        private:
            int m_numPoles;
        };

        //------------------------------------------------------------------------------

        // Factored implementations to reduce template instantiations

        struct LowPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                WorkspaceBase* w);
        };

        struct HighPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency,
                WorkspaceBase* w);
        };

        struct BandPassBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                WorkspaceBase* w);
        };

        struct BandStopBase : PoleFilterBase <AnalogLowPass>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency,
                WorkspaceBase* w);
        };

        //------------------------------------------------------------------------------

        //
        // Raw filters
        //

        template <int MaxOrder>
        struct LowPass : PoleFilter <LowPassBase, MaxOrder>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency)
            {
                Workspace <MaxOrder> w;
                LowPassBase::setup(order,
                    sampleRate,
                    cutoffFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct HighPass : PoleFilter <HighPassBase, MaxOrder>
        {
            void setup(int order,
                double sampleRate,
                double cutoffFrequency)
            {
                Workspace <MaxOrder> w;
                HighPassBase::setup(order,
                    sampleRate,
                    cutoffFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct BandPass : PoleFilter <BandPassBase, MaxOrder, MaxOrder * 2>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency)
            {
                Workspace <MaxOrder> w;
                BandPassBase::setup(order,
                    sampleRate,
                    centerFrequency,
                    widthFrequency,
                    &w);
            }
        };

        template <int MaxOrder>
        struct BandStop : PoleFilter <BandStopBase, MaxOrder, MaxOrder * 2>
        {
            void setup(int order,
                double sampleRate,
                double centerFrequency,
                double widthFrequency)
            {
                Workspace <MaxOrder> w;
                BandStopBase::setup(order,
                    sampleRate,
                    centerFrequency,
                    widthFrequency,
                    &w);
            }
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 3
                };

                static int getNumParams()
                {
                    return 3;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthHzParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(int(params[1]), params[0], params[2], params[3]);
                }
            };

            // Factored kind and name

            struct LowPassDescription
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "Legendre Low Pass"; }
            };

            struct HighPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Legendre High Pass"; }
            };

            struct BandPassDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Legendre Band Pass"; }
            };

            struct BandStopDescription
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "Legendre Band Stop"; }
            };

            // This glues on the Order parameter
            template <int MaxOrder,
                template <class> class TypeClass,
                template <int> class FilterClass>
            struct OrderBase : TypeClass <FilterClass <MaxOrder> >
            {
                const ParamInfo getParamInfo_1() const
                {
                    return ParamInfo(idOrder, "Order", "Order",
                        1, MaxOrder, 2,
                        &ParamInfo::Int_toControlValue,
                        &ParamInfo::Int_toNativeValue,
                        &ParamInfo::Int_toString);

                }
            };

            //------------------------------------------------------------------------------

            //
            // Design filters
            //

            template <int MaxOrder>
            struct LowPass : OrderBase <MaxOrder, TypeI, Legendre::LowPass>,
                LowPassDescription
            {
            };

            template <int MaxOrder>
            struct HighPass : OrderBase <MaxOrder, TypeI, Legendre::HighPass>,
                HighPassDescription
            {
            };

            template <int MaxOrder>
            struct BandPass : OrderBase <MaxOrder, TypeII, Legendre::BandPass>,
                BandPassDescription
            {
            };

            template <int MaxOrder>
            struct BandStop : OrderBase <MaxOrder, TypeII, Legendre::BandStop>,
                BandStopDescription
            {
            };

        }

    }

}

#endif
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_RBJ_H
#define DSPFILTERS_RBJ_H


namespace Dsp {

    /*
     * Filter realizations based on Robert Bristol-Johnson formulae:
     *
     * http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt
     *
     */

    namespace RBJ {

        //
        // Raw filters
        //

        struct LowPass : BiquadBase
        {
            void setup(double sampleRate,
                double cutoffFrequency,
                double q);
        };

        struct HighPass : BiquadBase
        {
            void setup(double sampleRate,
                double cutoffFrequency,
                double q);
        };

        struct BandPass1 : BiquadBase
        {
            // (constant skirt gain, peak gain = Q)
            void setup(double sampleRate,
                double centerFrequency,
                double bandWidth);
        };

        struct BandPass2 : BiquadBase
        {
            // (constant 0 dB peak gain)
            void setup(double sampleRate,
                double centerFrequency,
                double bandWidth);
        };

        struct BandStop : BiquadBase
        {
            void setup(double sampleRate,
                double centerFrequency,
                double bandWidth);
        };

        struct LowShelf : BiquadBase
        {
            void setup(double sampleRate,
                double cutoffFrequency,
                double gainDb,
                double shelfSlope);
        };

        struct HighShelf : BiquadBase
        {
            void setup(double sampleRate,
                double cutoffFrequency,
                double gainDb,
                double shelfSlope);
        };

        struct BandShelf : BiquadBase
        {
            void setup(double sampleRate,
                double centerFrequency,
                double gainDb,
                double bandWidth);
        };

        struct AllPass : BiquadBase
        {
            void setup(double sampleRate,
                double phaseFrequency,
                double q);
        };

        //------------------------------------------------------------------------------

        //
        // Gui-friendly Design layer
        //

        namespace Design {

            struct TypeIBase : DesignBase
            {
                enum
                {
                    NumParams = 3
                };

                static int getNumParams()
                {
                    return 3;
                }

                static const ParamInfo getParamInfo_1()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultQParam();
                }
            };

            template <class FilterClass>
            struct TypeI : TypeIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(params[0], params[1], params[2]);
                }
            };

            struct TypeIIBase : DesignBase
            {
                enum
                {
                    NumParams = 3
                };

                static int getNumParams()
                {
                    return 3;
                }

                static const ParamInfo getParamInfo_1()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultBandwidthParam();
                }
            };

            template <class FilterClass>
            struct TypeII : TypeIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(params[0], params[1], params[2]);
                }
            };

            struct TypeIIIBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_1()
                {
                    return ParamInfo::defaultCutoffFrequencyParam();
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultSlopeParam();
                }
            };

            template <class FilterClass>
            struct TypeIII : TypeIIIBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(params[0], params[1], params[2], params[3]);
                }
            };

            struct TypeIVBase : DesignBase
            {
                enum
                {
                    NumParams = 4
                };

                static int getNumParams()
                {
                    return 4;
                }

                static const ParamInfo getParamInfo_1()
                {
                    return ParamInfo::defaultCenterFrequencyParam();
                }

                static const ParamInfo getParamInfo_2()
                {
                    return ParamInfo::defaultGainParam();
                }

                static const ParamInfo getParamInfo_3()
                {
                    return ParamInfo::defaultBandwidthParam();
                }
            };

            template <class FilterClass>
            struct TypeIV : TypeIVBase, FilterClass
            {
                void setParams(const Params& params)
                {
                    FilterClass::setup(params[0], params[1], params[2], params[3]);
                }
            };

            //------------------------------------------------------------------------------

            struct LowPass : TypeI <RBJ::LowPass>
            {
                static Kind getKind() { return kindLowPass; }
                static const char* getName() { return "RBJ Low Pass"; }
            };

            struct HighPass : TypeI <RBJ::HighPass>
            {
                static Kind getKind() { return kindHighPass; }
                static const char* getName() { return "RBJ High Pass"; }
            };

            struct BandPass1 : TypeII <RBJ::BandPass1>
            {
                static Kind getKind() { return kindBandPass; }
                static const char* getName() { return "RBJ Band Pass 1"; }
            };

            struct BandPass2 : TypeII <RBJ::BandPass2>
            {
                static Kind getKind() { return kindBandPass; }
                static const char* getName() { return "RBJ Band Pass 2"; }
            };

            struct BandStop : TypeII <RBJ::BandStop>
            {
                static Kind getKind() { return kindBandStop; }
                static const char* getName() { return "RBJ Band Stop"; }
            };

            struct LowShelf : TypeIII <RBJ::LowShelf>
            {
                static Kind getKind() { return kindLowShelf; }
                static const char* getName() { return "RBJ Low Shelf"; }
            };

            struct HighShelf : TypeIII <RBJ::HighShelf>
            {
                static Kind getKind() { return kindHighShelf; }
                static const char* getName() { return "RBJ High Shelf"; }
            };

            struct BandShelf : TypeIV <RBJ::BandShelf>
            {
                static Kind getKind() { return kindBandShelf; }
                static const char* getName() { return "RBJ Band Shelf"; }
            };

            struct AllPass : TypeI <RBJ::AllPass>
            {
                static Kind getKind() { return kindOther; }
                static const char* getName() { return "RBJ All Pass"; }
            };

        }

    }

}

#endif

/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.cpp for contact information, notes, and bibliography.

--------------------------------------------------------------------------------

License: MIT License (http://www.opensource.org/licenses/mit-license.php)
Copyright (c) 2009 by Vinnie Falco

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

*******************************************************************************/

#ifndef DSPFILTERS_SMOOTHEDFILTER_H
#define DSPFILTERS_SMOOTHEDFILTER_H


namespace Dsp {

    /*
     * Implements smooth modulation of time-varying filter parameters
     *
     */
    template <class DesignClass,
        int Channels,
        class StateType = DirectFormII>
        class SmoothedFilterDesign
        : public FilterDesign <DesignClass,
        Channels,
        StateType>
    {
    public:
        typedef FilterDesign <DesignClass, Channels, StateType> filter_type_t;

        SmoothedFilterDesign(int transitionSamples)
            : m_transitionSamples(transitionSamples)
            , m_remainingSamples(-1) // first time flag
        {
        }

        // Process a block of samples.
        template <typename Sample>
        void processBlock(int numSamples,
            Sample* const* destChannelArray)
        {
            const int numChannels = this->getNumChannels();

            // If this goes off it means setup() was never called
            assert(m_remainingSamples >= 0);

            // first handle any transition samples
            int remainingSamples = std::min(m_remainingSamples, numSamples);

            if (remainingSamples > 0)
            {
                // interpolate parameters for each sample
                const double t = 1. / m_remainingSamples;
                double dp[maxParameters];
                for (int i = 0; i < DesignClass::NumParams; ++i)
                    dp[i] = (this->getParams()[i] - m_transitionParams[i]) * t;

                for (int n = 0; n < remainingSamples; ++n)
                {
                    for (int i = DesignClass::NumParams; --i >= 0;)
                        m_transitionParams[i] += dp[i];

                    m_transitionFilter.setParams(m_transitionParams);

                    for (int i = numChannels; --i >= 0;)
                    {
                        Sample* dest = destChannelArray[i] + n;
                        *dest = this->m_state[i].process(*dest, m_transitionFilter);
                    }
                }

                m_remainingSamples -= remainingSamples;

                if (m_remainingSamples == 0)
                    m_transitionParams = this->getParams();
            }

            // do what's left
            if (numSamples - remainingSamples > 0)
            {
                // no transition
                for (int i = 0; i < numChannels; ++i)
                    this->m_design.process(numSamples - remainingSamples,
                        destChannelArray[i] + remainingSamples,
                        this->m_state[i]);
            }
        }

        void process(int numSamples, float* const* arrayOfChannels)
        {
            processBlock(numSamples, arrayOfChannels);
        }

        void process(int numSamples, double* const* arrayOfChannels)
        {
            processBlock(numSamples, arrayOfChannels);
        }

    protected:
        void doSetParams(const Params& parameters)
        {
            if (m_remainingSamples >= 0)
            {
                m_remainingSamples = m_transitionSamples;
            }
            else
            {
                // first time
                m_remainingSamples = 0;
                m_transitionParams = parameters;
            }

            filter_type_t::doSetParams(parameters);
        }

    protected:
        Params m_transitionParams;
        DesignClass m_transitionFilter;
        int m_transitionSamples;

        int m_remainingSamples;        // remaining transition samples
    };

}

#endif



