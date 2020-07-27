
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

namespace Dsp {

    namespace Bessel {

        // returns fact(n) = n!
        static double fact(int n)
        {
            if (n == 0)
                return 1;

            double y = n;
            for (double m = n; --m;)
                y *= m;

            return y;
        }

        // returns the k-th zero based coefficient of the reverse bessel polynomial of degree n
        static double reversebessel(int k, int n)
        {
            return fact(2 * n - k) /
                ((fact(n - k) * fact(k)) * pow(2., n - k));
        }

        //------------------------------------------------------------------------------

        AnalogLowPass::AnalogLowPass()
            : m_numPoles(-1)
        {
            setNormal(0, 1);
        }

        void AnalogLowPass::design(int numPoles,
            WorkspaceBase* w)
        {
            if (m_numPoles != numPoles)
            {
                m_numPoles = numPoles;

                reset();

                RootFinderBase& solver(w->roots);
                for (int i = 0; i < numPoles + 1; ++i)
                    solver.coef()[i] = reversebessel(i, numPoles);
                solver.solve(numPoles);

                const int pairs = numPoles / 2;
                for (int i = 0; i < pairs; ++i)
                {
                    complex_t c = solver.root()[i];
                    addPoleZeroConjugatePairs(c, infinity());
                }

                if (numPoles & 1)
                    add(solver.root()[pairs].real(), infinity());
            }
        }

        //------------------------------------------------------------------------------

        AnalogLowShelf::AnalogLowShelf()
            : m_numPoles(-1)
        {
            setNormal(doublePi, 1);
        }

        void AnalogLowShelf::design(int numPoles,
            double gainDb,
            WorkspaceBase* w)
        {
            if (m_numPoles != numPoles ||
                m_gainDb != gainDb)
            {
                m_numPoles = numPoles;
                m_gainDb = gainDb;

                reset();

                const double G = pow(10., gainDb / 20) - 1;

                RootFinderBase& poles(w->roots);
                for (int i = 0; i < numPoles + 1; ++i)
                    poles.coef()[i] = reversebessel(i, numPoles);
                poles.solve(numPoles);

                RootFinder<50> zeros;
                for (int i = 0; i < numPoles + 1; ++i)
                    zeros.coef()[i] = reversebessel(i, numPoles);
                double a0 = reversebessel(0, numPoles);
                zeros.coef()[0] += G * a0;
                zeros.solve(numPoles);

                const int pairs = numPoles / 2;
                for (int i = 0; i < pairs; ++i)
                {
                    complex_t p = poles.root()[i];
                    complex_t z = zeros.root()[i];
                    addPoleZeroConjugatePairs(p, z);
                }

                if (numPoles & 1)
                    add(poles.root()[pairs].real(), zeros.root()[pairs].real());
            }
        }

        //------------------------------------------------------------------------------

        void LowPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandPassBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandStopBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            BandStopTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void LowShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, gainDb, w);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

    }

}
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

namespace Dsp {

    BiquadPoleState::BiquadPoleState(const BiquadBase& s)
    {
        const double a0 = s.getA0();
        const double a1 = s.getA1();
        const double a2 = s.getA2();
        const double b0 = s.getB0();
        const double b1 = s.getB1();
        const double b2 = s.getB2();

        if (a2 == 0 && b2 == 0)
        {
            // single pole
            poles.first = -a1;
            zeros.first = -b0 / b1;
            poles.second = 0;
            zeros.second = 0;
        }
        else
        {
            {
                const complex_t c = sqrt(complex_t(a1 * a1 - 4 * a0 * a2, 0));
                double d = 2. * a0;
                poles.first = -(a1 + c) / d;
                poles.second = (c - a1) / d;
                assert(!poles.is_nan());
            }

            {
                const complex_t c = sqrt(complex_t(
                    b1 * b1 - 4 * b0 * b2, 0));
                double d = 2. * b0;
                zeros.first = -(b1 + c) / d;
                zeros.second = (c - b1) / d;
                assert(!zeros.is_nan());
            }
        }

        gain = b0 / a0;
    }

    //------------------------------------------------------------------------------

    complex_t BiquadBase::response(double normalizedFrequency) const
    {
        const double a0 = getA0();
        const double a1 = getA1();
        const double a2 = getA2();
        const double b0 = getB0();
        const double b1 = getB1();
        const double b2 = getB2();

        const double w = 2 * doublePi * normalizedFrequency;
        const complex_t czn1 = std::polar(1., -w);
        const complex_t czn2 = std::polar(1., -2 * w);
        complex_t ch(1);
        complex_t cbot(1);

        complex_t ct(b0 / a0);
        complex_t cb(1);
        ct = addmul(ct, b1 / a0, czn1);
        ct = addmul(ct, b2 / a0, czn2);
        cb = addmul(cb, a1 / a0, czn1);
        cb = addmul(cb, a2 / a0, czn2);
        ch *= ct;
        cbot *= cb;

        return ch / cbot;
    }

    std::vector<PoleZeroPair> BiquadBase::getPoleZeros() const
    {
        std::vector<PoleZeroPair> vpz;
        BiquadPoleState bps(*this);
        vpz.push_back(bps);
        return vpz;
    }

    void BiquadBase::setCoefficients(double a0, double a1, double a2,
        double b0, double b1, double b2)
    {
        assert(!Dsp::is_nan(a0) && !Dsp::is_nan(a1) && !Dsp::is_nan(a2) &&
            !Dsp::is_nan(b0) && !Dsp::is_nan(b1) && !Dsp::is_nan(b2));

        m_a0 = a0;
        m_a1 = a1 / a0;
        m_a2 = a2 / a0;
        m_b0 = b0 / a0;
        m_b1 = b1 / a0;
        m_b2 = b2 / a0;
    }

    void BiquadBase::setOnePole(complex_t pole, complex_t zero)
    {
#if 0
        pole = adjust_imag(pole);
        zero = adjust_imag(zero);
#else
        assert(pole.imag() == 0);
        assert(zero.imag() == 0);
#endif

        const double a0 = 1;
        const double a1 = -pole.real();
        const double a2 = 0;
        const double b0 = -zero.real();
        const double b1 = 1;
        const double b2 = 0;

        setCoefficients(a0, a1, a2, b0, b1, b2);
    }

    void BiquadBase::setTwoPole(complex_t pole1, complex_t zero1,
        complex_t pole2, complex_t zero2)
    {
#if 0
        pole1 = adjust_imag(pole1);
        pole2 = adjust_imag(pole2);
        zero1 = adjust_imag(zero1);
        zero2 = adjust_imag(zero2);
#endif

        const double a0 = 1;
        double a1;
        double a2;

        if (pole1.imag() != 0)
        {
            assert(pole2 == std::conj(pole1));

            a1 = -2 * pole1.real();
            a2 = std::norm(pole1);
        }
        else
        {
            assert(pole2.imag() == 0);

            a1 = -(pole1.real() + pole2.real());
            a2 = pole1.real() * pole2.real();
        }

        const double b0 = 1;
        double b1;
        double b2;

        if (zero1.imag() != 0)
        {
            assert(zero2 == std::conj(zero1));

            b1 = -2 * zero1.real();
            b2 = std::norm(zero1);
        }
        else
        {
            assert(zero2.imag() == 0);

            b1 = -(zero1.real() + zero2.real());
            b2 = zero1.real() * zero2.real();
        }

        setCoefficients(a0, a1, a2, b0, b1, b2);
    }

    void BiquadBase::setPoleZeroForm(const BiquadPoleState& bps)
    {
        setPoleZeroPair(bps);
        applyScale(bps.gain);
    }

    void BiquadBase::setIdentity()
    {
        setCoefficients(1, 0, 0, 1, 0, 0);
    }

    void BiquadBase::applyScale(double scale)
    {
        m_b0 *= scale;
        m_b1 *= scale;
        m_b2 *= scale;
    }

    //------------------------------------------------------------------------------

    Biquad::Biquad()
    {
    }

    // Construct a second order section from a pair of poles and zeroes
    Biquad::Biquad(const BiquadPoleState& bps)
    {
        setPoleZeroForm(bps);
    }

    //------------------------------------------------------------------------------

}
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

namespace Dsp {

    namespace Butterworth {

        AnalogLowPass::AnalogLowPass()
            : m_numPoles(-1)
        {
            setNormal(0, 1);
        }

        void AnalogLowPass::design(int numPoles)
        {
            if (m_numPoles != numPoles)
            {
                m_numPoles = numPoles;

                reset();

                const double n2 = 2 * numPoles;
                const int pairs = numPoles / 2;
                for (int i = 0; i < pairs; ++i)
                {
                    complex_t c = std::polar(1., doublePi_2 + (2 * i + 1) * doublePi / n2);
                    addPoleZeroConjugatePairs(c, infinity());
                }

                if (numPoles & 1)
                    add(-1, infinity());
            }
        }

        //------------------------------------------------------------------------------

        AnalogLowShelf::AnalogLowShelf()
            : m_numPoles(-1)
        {
            setNormal(doublePi, 1);
        }

        void AnalogLowShelf::design(int numPoles, double gainDb)
        {
            if (m_numPoles != numPoles ||
                m_gainDb != gainDb)
            {
                m_numPoles = numPoles;
                m_gainDb = gainDb;

                reset();

                const double n2 = numPoles * 2;
                const double g = pow(pow(10., gainDb / 20), 1. / n2);
                const double gp = -1. / g;
                const double gz = -g;

                const int pairs = numPoles / 2;
                for (int i = 1; i <= pairs; ++i)
                {
                    const double theta = doublePi * (0.5 - (2 * i - 1) / n2);
                    addPoleZeroConjugatePairs(std::polar(gp, theta), std::polar(gz, theta));
                }

                if (numPoles & 1)
                    add(gp, gz);
            }
        }

        //------------------------------------------------------------------------------

        void LowPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency)
        {
            m_analogProto.design(order);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency)
        {
            m_analogProto.design(order);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandPassBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency)
        {
            m_analogProto.design(order);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandStopBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency)
        {
            m_analogProto.design(order);

            BandStopTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void LowShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb)
        {
            m_analogProto.design(order, gainDb);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb)
        {
            m_analogProto.design(order, gainDb);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandShelfBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double gainDb)
        {
            m_analogProto.design(order, gainDb);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            // HACK!
            m_digitalProto.setNormal(((centerFrequency / sampleRate) < 0.25) ? doublePi : 0, 1);

            Cascade::setLayout(m_digitalProto);
        }

    }

}
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

namespace Dsp {

    Cascade::Cascade()
        : m_numStages(0)
        , m_maxStages(0)
        , m_stageArray(0)
    {
    }

    void Cascade::setCascadeStorage(const Storage& storage)
    {
        m_numStages = 0;
        m_maxStages = storage.maxStages;
        m_stageArray = storage.stageArray;
    }

    complex_t Cascade::response(double normalizedFrequency) const
    {
        double w = 2 * doublePi * normalizedFrequency;
        const complex_t czn1 = std::polar(1., -w);
        const complex_t czn2 = std::polar(1., -2 * w);
        complex_t ch(1);
        complex_t cbot(1);

        const Biquad* stage = m_stageArray;
        for (int i = m_numStages; --i >= 0; ++stage)
        {
            complex_t cb(1);
            complex_t ct(stage->getB0() / stage->getA0());
            ct = addmul(ct, stage->getB1() / stage->getA0(), czn1);
            ct = addmul(ct, stage->getB2() / stage->getA0(), czn2);
            cb = addmul(cb, stage->getA1() / stage->getA0(), czn1);
            cb = addmul(cb, stage->getA2() / stage->getA0(), czn2);
            ch *= ct;
            cbot *= cb;
        }

        return ch / cbot;
    }

    std::vector<PoleZeroPair> Cascade::getPoleZeros() const
    {
        std::vector<PoleZeroPair> vpz;
        vpz.reserve(m_numStages);

        const Stage* stage = m_stageArray;
        for (int i = m_numStages; --i >= 0;)
        {
            BiquadPoleState bps(*stage++);
            assert(!bps.isSinglePole() || i == 0);
            vpz.push_back(bps);
        }

        return vpz;
    }

    void Cascade::applyScale(double scale)
    {
        // For higher order filters it might be helpful
        // to spread this factor between all the stages.
        assert(m_numStages > 0);
        m_stageArray->applyScale(scale);
    }

    void Cascade::setLayout(const LayoutBase& proto)
    {
        const int numPoles = proto.getNumPoles();
        m_numStages = (numPoles + 1) / 2;
        assert(m_numStages <= m_maxStages);

        Biquad* stage = m_stageArray;
        for (int i = 0; i < m_numStages; ++i, ++stage)
            stage->setPoleZeroPair(proto[i]);

        applyScale(proto.getNormalGain() /
            std::abs(response(proto.getNormalW() / (2 * doublePi))));
    }

}

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

namespace Dsp {

    namespace ChebyshevI {

        AnalogLowPass::AnalogLowPass()
            : m_numPoles(-1)
        {
        }

        void AnalogLowPass::design(int numPoles,
            double rippleDb)
        {
            if (m_numPoles != numPoles ||
                m_rippleDb != rippleDb)
            {
                m_numPoles = numPoles;
                m_rippleDb = rippleDb;

                reset();

                const double eps = std::sqrt(1. / std::exp(-rippleDb * 0.1 * doubleLn10) - 1);
                const double v0 = asinh(1 / eps) / numPoles;
                const double sinh_v0 = -sinh(v0);
                const double cosh_v0 = cosh(v0);

                const double n2 = 2 * numPoles;
                const int pairs = numPoles / 2;
                for (int i = 0; i < pairs; ++i)
                {
                    const int k = 2 * i + 1 - numPoles;
                    double a = sinh_v0 * cos(k * doublePi / n2);
                    double b = cosh_v0 * sin(k * doublePi / n2);

                    //addPoleZero (complex_t (a, b), infinity());
                    //addPoleZero (complex_t (a, -b), infinity());
                    addPoleZeroConjugatePairs(complex_t(a, b), infinity());
                }

                if (numPoles & 1)
                {
                    add(complex_t(sinh_v0, 0), infinity());
                    setNormal(0, 1);
                }
                else
                {
                    setNormal(0, pow(10, -rippleDb / 20.));
                }
            }
        }

        //------------------------------------------------------------------------------

        //
        // Chebyshev Type I low pass shelf prototype
        // From "High-Order Digital Parametric Equalizer Design"
        // Sophocles J. Orfanidis
        // http://www.ece.rutgers.edu/~orfanidi/ece521/hpeq.pdf
        //

        AnalogLowShelf::AnalogLowShelf()
        {
            setNormal(doublePi, 1);
        }

        void AnalogLowShelf::design(int numPoles,
            double gainDb,
            double rippleDb)
        {
            if (m_numPoles != numPoles ||
                m_rippleDb != rippleDb ||
                m_gainDb != gainDb)
            {
                m_numPoles = numPoles;
                m_rippleDb = rippleDb;
                m_gainDb = gainDb;

                reset();

                gainDb = -gainDb;

                if (rippleDb >= fabs(gainDb))
                    rippleDb = fabs(gainDb);
                if (gainDb < 0)
                    rippleDb = -rippleDb;

                const double G = std::pow(10., gainDb / 20.0);
                const double Gb = std::pow(10., (gainDb - rippleDb) / 20.0);
                const double G0 = 1;
                const double g0 = pow(G0, 1. / numPoles);

                double eps;
                if (Gb != G0)
                    eps = sqrt((G * G - Gb * Gb) / (Gb * Gb - G0 * G0));
                else
                    eps = G - 1; // This is surely wrong

                const double b = pow(G / eps + Gb * sqrt(1 + 1 / (eps * eps)), 1. / numPoles);
                const double u = log(b / g0);
                const double v = log(pow(1. / eps + sqrt(1 + 1 / (eps * eps)), 1. / numPoles));

                const double sinh_u = sinh(u);
                const double sinh_v = sinh(v);
                const double cosh_u = cosh(u);
                const double cosh_v = cosh(v);
                const double n2 = 2 * numPoles;
                const int pairs = numPoles / 2;
                for (int i = 1; i <= pairs; ++i)
                {
                    const double a = doublePi * (2 * i - 1) / n2;
                    const double sn = sin(a);
                    const double cs = cos(a);
                    addPoleZeroConjugatePairs(complex_t(-sn * sinh_u, cs * cosh_u),
                        complex_t(-sn * sinh_v, cs * cosh_v));
                }

                if (numPoles & 1)
                    add(-sinh_u, -sinh_v);
            }
        }

        //------------------------------------------------------------------------------

        void LowPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double rippleDb)
        {
            m_analogProto.design(order, rippleDb);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double rippleDb)
        {
            m_analogProto.design(order, rippleDb);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandPassBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double rippleDb)
        {
            m_analogProto.design(order, rippleDb);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandStopBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double rippleDb)
        {
            m_analogProto.design(order, rippleDb);

            BandStopTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void LowShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb,
            double rippleDb)
        {
            m_analogProto.design(order, gainDb, rippleDb);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb,
            double rippleDb)
        {
            m_analogProto.design(order, gainDb, rippleDb);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandShelfBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double gainDb,
            double rippleDb)
        {
            m_analogProto.design(order, gainDb, rippleDb);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            m_digitalProto.setNormal(((centerFrequency / sampleRate) < 0.25) ? doublePi : 0, 1);

            Cascade::setLayout(m_digitalProto);
        }

    }

}
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

namespace Dsp {

    namespace ChebyshevII {

        // "Chebyshev Filter Properties"
        // http://cnx.org/content/m16906/latest/

        AnalogLowPass::AnalogLowPass()
            : m_numPoles(-1)
        {
            setNormal(0, 1);
        }

        void AnalogLowPass::design(int numPoles,
            double stopBandDb)
        {
            if (m_numPoles != numPoles ||
                m_stopBandDb != stopBandDb)
            {
                m_numPoles = numPoles;
                m_stopBandDb = stopBandDb;

                reset();

                const double eps = std::sqrt(1. / (std::exp(stopBandDb * 0.1 * doubleLn10) - 1));
                const double v0 = asinh(1 / eps) / numPoles;
                const double sinh_v0 = -sinh(v0);
                const double cosh_v0 = cosh(v0);
                const double fn = doublePi / (2 * numPoles);

                int k = 1;
                for (int i = numPoles / 2; --i >= 0; k += 2)
                {
                    const double a = sinh_v0 * cos((k - numPoles) * fn);
                    const double b = cosh_v0 * sin((k - numPoles) * fn);
                    const double d2 = a * a + b * b;
                    const double im = 1 / cos(k * fn);
                    addPoleZeroConjugatePairs(complex_t(a / d2, b / d2),
                        complex_t(0, im));
                }

                if (numPoles & 1)
                {
                    add(1 / sinh_v0, infinity());
                }
            }
        }

        //------------------------------------------------------------------------------

        //
        // Chebyshev Type I low pass shelf prototype
        // From "High-Order Digital Parametric Equalizer Design"
        // Sophocles J. Orfanidis
        // http://www.ece.rutgers.edu/~orfanidi/ece521/hpeq.pdf
        //

        AnalogLowShelf::AnalogLowShelf()
            : m_numPoles(-1)
        {
            setNormal(doublePi, 1);
        }

        void AnalogLowShelf::design(int numPoles,
            double gainDb,
            double stopBandDb)
        {
            if (m_numPoles != numPoles ||
                m_stopBandDb != stopBandDb ||
                m_gainDb != gainDb)
            {
                m_numPoles = numPoles;
                m_stopBandDb = stopBandDb;
                m_gainDb = gainDb;

                reset();

                gainDb = -gainDb;

                if (stopBandDb >= fabs(gainDb))
                    stopBandDb = fabs(gainDb);
                if (gainDb < 0)
                    stopBandDb = -stopBandDb;

                const double G = std::pow(10., gainDb / 20.0);
                const double Gb = std::pow(10., (gainDb - stopBandDb) / 20.0);
                const double G0 = 1;
                const double g0 = pow(G0, 1. / numPoles);

                double eps;
                if (Gb != G0)
                    eps = sqrt((G * G - Gb * Gb) / (Gb * Gb - G0 * G0));
                else
                    eps = G - 1; // This is surely wrong

                const double b = pow(G / eps + Gb * sqrt(1 + 1 / (eps * eps)), 1. / numPoles);
                const double u = log(b / g0);
                const double v = log(pow(1. / eps + sqrt(1 + 1 / (eps * eps)), 1. / numPoles));

                const double sinh_u = sinh(u);
                const double sinh_v = sinh(v);
                const double cosh_u = cosh(u);
                const double cosh_v = cosh(v);
                const double n2 = 2 * numPoles;
                const int pairs = numPoles / 2;
                for (int i = 1; i <= pairs; ++i)
                {
                    const double a = doublePi * (2 * i - 1) / n2;
                    const double sn = sin(a);
                    const double cs = cos(a);
                    addPoleZeroConjugatePairs(complex_t(-sn * sinh_u, cs * cosh_u),
                        complex_t(-sn * sinh_v, cs * cosh_v));
                }

                if (numPoles & 1)
                    add(-sinh_u, -sinh_v);
            }
        }

        //------------------------------------------------------------------------------

        void LowPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double stopBandDb)
        {
            m_analogProto.design(order, stopBandDb);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double stopBandDb)
        {
            m_analogProto.design(order, stopBandDb);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandPassBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double stopBandDb)
        {
            m_analogProto.design(order, stopBandDb);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandStopBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double stopBandDb)
        {
            m_analogProto.design(order, stopBandDb);

            BandStopTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void LowShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb,
            double stopBandDb)
        {
            m_analogProto.design(order, gainDb, stopBandDb);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighShelfBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double gainDb,
            double stopBandDb)
        {
            m_analogProto.design(order, gainDb, stopBandDb);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandShelfBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double gainDb,
            double stopBandDb)
        {
            m_analogProto.design(order, gainDb, stopBandDb);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            m_digitalProto.setNormal(((centerFrequency / sampleRate) < 0.25) ? doublePi : 0, 1);

            Cascade::setLayout(m_digitalProto);
        }

    }

}
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

namespace Dsp {

    namespace Custom {

        void OnePole::setup(double scale,
            double pole,
            double zero)
        {
            setOnePole(pole, zero);
            applyScale(scale);
        }

        void TwoPole::setup(double scale,
            double poleRho,
            double poleTheta,
            double zeroRho,
            double zeroTheta)
        {
            complex_t pole = std::polar(poleRho, poleTheta);
            complex_t zero = std::polar(zeroRho, zeroTheta);

            setTwoPole(pole, zero, std::conj(pole), std::conj(zero));
            applyScale(scale);
        }

    }

}
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

namespace Dsp {

}
/*******************************************************************************

"A Collection of Useful C++ Classes for Digital Signal Processing"
 By Vinnie Falco

Official project location:
https://github.com/vinniefalco/DSPFilters

See Documentation.h for contact information, notes, and bibliography.

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

********************************************************************************

="A Collection of Useful C++ Classes for Digital Signal Processing"=
_By Vincent Falco_

_"Techniques for digital signal processing are well guarded and held
close to the chest, as they have valuable applications for multimedia
content. The black art of Infinite Impulse Response ("IIR") filtering
has remained shrouded in secrecy with little publicly available source
code....until now."_

----
Building on the work of cherished luminaries such as Sophocles Orfanidis, Andreas Antoniou, Martin Holters, and Udo Zolzer, this library harnesses the power of C++ templates to solve a useful problem in digital signal processing: the realization of multichannel IIR filters of arbitrary order and prescribed specifications with various properties such as Butterworth, Chebyshev, Elliptic, and Optimum-L (Legendre) responses. The library is provided under the MIT license and is therefore fully compatible with proprietary usage.

Classes are designed as independent re-usable building blocks. Use some or all of the provided features, or extend the functionality by writing your own objects that plug into the robust framework. Only the code that you need will get linked into your application. Here's a list of features:

    * Exclusive focus on IIR filters instead of boring FIR filters
    * Complete implementation of all "RBJ Biquad" Cookbook filter formulas
    * Butterworth, Chebyshev, Elliptic, Bessel, Legendre designs
    * Low Pass, High Pass, Band Pass, Band Stop transformations
    * Low, High, and Band Shelf filter implementations for most types
    * Smooth interpolation of filter settings, pole/zeros, and biquad coefficients to achieve seamless parameter changes
    * Representation of digital filters using poles and zeros
    * Realization using Direct Form I, Direct Form II, or user provided class
    * Fully factored to minimize template instantiations
    * "Design" layer provides runtime introspection into a filter
    * Utility template functions for manipulating buffers of sample data
    * No calls to malloc or new, great for embedded systems
    * No external dependencies, just the standard C++ library!

Using these filters is easy:

{{{
    // Create a Chebyshev type I Band Stop filter of order 3
    // with state for processing 2 channels of audio.
    Dsp::SimpleFilter <Dsp::ChebyshevI::BandStop <3>, 2> f;
    f.setup (3,    // order
             44100,// sample rate
             4000, // center frequency
             880,  // band width
             1);   // ripple dB
    f.process (numSamples, arrayOfChannels);
}}}

An accompanying demonstration program that works on most popular platforms by using the separately licensed Juce application framework (included), exercises all the functionality of the library, including these features:

    * Dynamic interface creates itself using filter introspection capabilities
    * Audio playback with real time application of a selected filter
    * Live time stretching and amplitude modulation without clicks or popping
    * Charts to show magnitude, phase response and pole/zero placement
    * Thread safety "best practices" for audio applications

Here's a screenshot of the DspFilters Demo

http://dspfilterscpp.googlecode.com/files/dspfiltersdemo.png

If you've been searching in futility on the Internet for some source code for implementing high order filters, then look no further because this is it! Whether you are a student of C++ or digital signal processing, a writer of audio plugins, or even a VST synthesizer coder, "A Collection of Useful C++ Classes for Digital Signal Processing" might have something for you!

********************************************************************************

Notes:

  Please direct all comments this DSP and Plug-in Development forum:

  http://www.kvraudio.com/forum/viewforum.php?f=33

Credits

  All of this code was written by the author Vinnie Falco except where marked.

  Some filter ideas are based on a java applet (http://www.falstad.com/dfilter/)
  developed by Paul Falstad.

Bibliography

  "High-Order Digital Parametric Equalizer Design"
   Sophocles J. Orfanidis
   (Journal of the Audio Engineering Society, vol 53. pp 1026-1046)

  http://crca.ucsd.edu/~msp/techniques/v0.08/book-html/node1.html

  "Spectral Transformations for digital filters"
   A. G. Constantinides, B.Sc.(Eng.) Ph.D.
   (Proceedings of the IEEE, vol. 117, pp. 1585-1590, August 1970)

********************************************************************************

DOCUMENTATION

All symbols are in the Dsp namespace.

class Filter

  This is an abstract polymorphic interface that supports any filter. The
  parameters to the filter are passed in the Params structure, which is
  essentially an array of floating point numbers with a hard coded size
  limit (maxParameters). Each filter makes use of the Params as it sees fit.

  Filter::getKind ()
  Filter::getName ()
  Filter::getNumParams ()
  Filter::getParamInfo ()

  Through the use of these functions, the caller can determine the meaning
  of each indexed filter parameter at run-time. The ParamInfo structure
  contains methods that describe information about an individual parameter,
  including convenience functions to map a filter parameter to a "control
  value" in the range 0...1, suitable for presentation by a GUI element such
  as a knob or scrollbar.

  Filter::getDefaultParams ()
  Filter::getParams ()
  Filter::getParam ()
  Filter::setParam ()
  Filter::findParamId ()
  Filter::setParamById ()
  Filter::setParams ()
  Filter::copyParamsFrom ()

  These methods allow the caller to inspect the values of the parameters,
  and set the filter parameters in various ways. When parameters are changed
  they take effect on the filter immediately.

  Filter::getPoleZeros ()
  Filter::response ()

  For analysis, these routines provide insight into the pole/zero arrangement
  in the z-plane, and the complex valued response at a given normalized
  frequency in the range (0..nyquist = 0.5]. From the complex number the
  magnitude and phase can be calculated.

  Filter::getNumChannels()
  Filter::reset()
  Filter::process()

  These functions are for applying the filter to channels of data. If the
  filter was not created with channel state (i.e. Channels==0 in the derived
  class template) then they will throw an exception.

  To create a Filter object, use operator new on a subclass template with
  appropriate parameters based on the type of filter you want. Here are the
  subclasses.



template <class DesignClass, int Channels = 0, class StateType = DirectFormII>
class FilterDesign : public Filter

  This subclass of Filter takes a DesignClass (explained below) representing
  a filter, an optional parameter for the number of channels of data to
  process, and an optional customizable choice of which state realization
  to use for processing samples. Channels may be zero, in which case the
  object can only be used for analysis.

  Because the DesignClass is a member and not inherited, it is in general
  not possible to call members of the DesignClass directly. You must go
  through the Filter interface.



template <class DesignClass, int Channels, class StateType = DirectFormII>
class SmoothedFilterDesign : public Filter

  This subclass of FilterDesign implements a filter of the given DesignClass,
  and also performs smoothing of parameters over time. Specifically, when
  one or more filter parameters (such as cutoff frequency) are changed, the
  class creates a transition over a given number of samples from the original
  values to the new values. This process is invisible and seamless to the
  caller, except that the constructor takes an additional parameter that
  indicates the duration of transitions when parameters change.



template <class FilterClass, int Channels = 0, class StateType = DirectFormII>
class SimpleFilter : public FilterClass

  This is a simple wrapper around a given raw FilterClass (explained below).
  It uses inheritance so all of the members of the FilterClass are available
  to instances of this object. The simple wrapper provides state information
  for processing channels in the given form.

  The wrapper does not support introspection, parameter smoothing, or the
  Params style of applying filter settings. Instead, it uses the interface
  of the given FilterClass, which is typically a function called setup()
  that takes a list of arguments representing the parameters.

  The use of this class bypasses the virtual function overhead of going
  through a Filter object. It is not practical to change filter parameters
  of a SimpleFilter, unless you are re-using the filter for a brand new
  stream of data in which case reset() should be called immediately before
  or after changing parameters, to clear the state and prevent audible
  artifacts.



Filter family namespaces

  Each family of filters is given its own namespace. Currently these namespaces
  include:

  RBJ:          Filters from the RBJ Cookbook
  Butterworth:  Filters with Butterworth response
  ChebyshevI:   Filters using Chebyshev polynomials (ripple in the passband)
  ChebyshevII:  "Inverse Chebyshev" filters (ripple in the stopband)
  Elliptic:     Filters with ripple in both the passband and stopband
  Bessel:       Uses Bessel polynomials, theoretically with linear phase
  Legendre:     "Optimum-L" filters with steepest transition and monotonic passband.
  Custom:       Simple filters that allow poles and zeros to be specified directly

<class FilterClass>

  Within each namespace we have a set of "raw filters" (each one is an example
  of a FilterClass). For example, the raw filters in the Butterworth namespace are:

  Butterworth::LowPass
  Butterworth::HighPass
  Butterworth::BandPass
  Butterworth::BandStop
  Butterworth::LowShelf
  Butterworth::HighShelf
  Butterworth::BandShelf

  When a class template (such as SimpleFilter) requires a FilterClass, it is
  expecting an identifier of a raw filter. For example, Legendre::LowPass. The
  raw filters do not have any support for introspection or the Params style of
  changing filter settings. All they offer is a setup() function for updating
  the IIR coefficients to a given set of parameters.

<class DesignClass>

  Each filter family namespace also has the nested namespace "Design". Inside
  this namespace we have all of the raw filter names repeated, except that
  these classes additional provide the Design interface, which adds
  introspection, polymorphism, the Params style of changing filter settings,
  and in general all of the features necessary to interoperate with the Filter
  virtual base class and its derived classes. For example, the design filters
  from the Butterworth namespace are:

  Butterworth::Design::LowPass
  Butterworth::Design::HighPass
  Butterworth::Design::BandPass
  Butterworth::Design::BandStop
  Butterworth::Design::LowShelf
  Butterworth::Design::HighShelf
  Butterworth::Design::BandShelf

  For any class template that expects a DesignClass, you must pass a suitable
  object from the Design namespace of the desired filter family. For example,
  ChebyshevI::Design::BandPass.

*******************************************************************************/

//
// Usage Examples
//
// This shows you how to operate the filters
//


// This is the only include you need

#include <sstream>
#include <iostream>
#include <iomanip>

namespace {

    void UsageExamples()
    {
        // create a two channel audio buffer
        int numSamples = 2000;
        float* audioData[2];
        audioData[0] = new float[numSamples];
        audioData[1] = new float[numSamples];

        // create a 2-channel RBJ Low Pass with parameter smoothing
        // and apply it to the audio data
        {
            // "1024" is the number of samples over which to fade parameter changes
            Dsp::Filter* f = new Dsp::SmoothedFilterDesign
                <Dsp::RBJ::Design::LowPass, 2>(1024);
            Dsp::Params params;
            params[0] = 44100; // sample rate
            params[1] = 4000; // cutoff frequency
            params[2] = 1.25; // Q
            f->setParams(params);
            f->process(numSamples, audioData);
        }

        // set up a 2-channel RBJ High Pass with parameter smoothing,
        // but bypass virtual function overhead
        {
            // the difference here is that we don't go through a pointer.
            Dsp::SmoothedFilterDesign <Dsp::RBJ::Design::LowPass, 2> f(1024);
            Dsp::Params params;
            params[0] = 44100; // sample rate
            params[1] = 4000; // cutoff frequency
            params[2] = 1.25; // Q
            f.setParams(params);
            f.process(numSamples, audioData);
        }

        // create a 2-channel Butterworth Band Pass of order 4,
        // with parameter smoothing and apply it to the audio data.
        // Output samples are generated using Direct Form II realization.
        {
            Dsp::Filter* f = new Dsp::SmoothedFilterDesign
                <Dsp::Butterworth::Design::BandPass <4>, 2, Dsp::DirectFormII>(1024);
            Dsp::Params params;
            params[0] = 44100; // sample rate
            params[1] = 4; // order
            params[2] = 4000; // center frequency
            params[3] = 880; // band width
            f->setParams(params);
            f->process(numSamples, audioData);
        }

        // create a 2-channel Inverse Chebyshev Low Shelf of order 5
        // and passband ripple 0.1dB, without parameter smoothing and apply it.
        {
            Dsp::Filter* f = new Dsp::FilterDesign
                <Dsp::ChebyshevII::Design::LowShelf <5>, 2>;
            Dsp::Params params;
            params[0] = 44100; // sample rate
            params[1] = 5; // order
            params[2] = 4000; // corner frequency
            params[3] = 6; // shelf gain
            params[4] = 0.1; // passband ripple
            f->setParams(params);
            f->process(numSamples, audioData);
        }

        // create an abstract Butterworth High Pass of order 4.
        // This one can't process channels, it can only be used for analysis
        // (i.e. extract poles and zeros).
        {
            Dsp::Filter* f = new Dsp::FilterDesign
                <Dsp::Butterworth::Design::HighPass <4> >;
            Dsp::Params params;
            params[0] = 44100; // sample rate
            params[1] = 4; // order
            params[2] = 4000; // cutoff frequency
            f->setParams(params);
            // this will cause a runtime assertion
            f->process(numSamples, audioData);
        }

        // Use the simple filter API to create a Chebyshev Band Stop of order 3
        // and 1dB ripple in the passband. The simle API has a smaller
        // footprint, but no introspection or smoothing.
        {
            // Note we use the raw filter instead of the one
            // from the Design namespace.
            Dsp::SimpleFilter <Dsp::ChebyshevI::BandStop <3>, 2> f;
            f.setup(3,    // order
                44100,// sample rate
                4000, // center frequency
                880,  // band width
                1);   // ripple dB
            f.process(numSamples, audioData);
        }

        // Set up a filter, extract the coefficients and print them to standard
        // output. Note that this filter is not capable of processing samples,
        // as it has no state. It only has coefficients.
        {
            Dsp::SimpleFilter <Dsp::RBJ::LowPass> f;
            f.setup(44100, // sample rate Hz
                440, // cutoff frequency Hz
                1); // "Q" (resonance)

            std::ostringstream os;

            os << "a0 = " << f.getA0() << "\n"
                << "a1 = " << f.getA1() << "\n"
                << "a2 = " << f.getA2() << "\n"
                << "b0 = " << f.getB0() << "\n"
                << "b1 = " << f.getB1() << "\n"
                << "b2 = " << f.getB2() << "\n";
            ;

            std::cout << os.str();
        }

        // Create an instance of a raw filter. This is as low as it gets, any
        // lower and we will just have either a Biquad or a Cascade, and you'll
        // be setting the coefficients manually.
        {
            // This is basically like eating uncooked food
            Dsp::RBJ::LowPass f;
            f.setup(44100, 440, 1);

            // calculate response at frequency 440 Hz
            Dsp::complex_t response = f.response(440. / 44100);
        }

        // Extract coefficients from a Cascade
        {
            Dsp::SimpleFilter <Dsp::Butterworth::HighPass <3> > f;
            f.setup(3, 44100, 2000);

            std::ostringstream os;

            os << "numStages = " << f.getNumStages() << "\n"
                << "a0[0] = " << f[0].getA0() << "\n"
                << "a1[0] = " << f[0].getA1() << "\n"
                << "a2[0] = " << f[0].getA2() << "\n"
                << "b0[0] = " << f[0].getB0() << "\n"
                << "b1[0] = " << f[0].getB1() << "\n"
                << "b2[0] = " << f[0].getB2() << "\n"
                << "a0[1] = " << f[1].getA0() << "\n"
                << "a1[1] = " << f[1].getA1() << "\n"
                << "a2[1] = " << f[1].getA2() << "\n"
                << "b0[1] = " << f[1].getB0() << "\n"
                << "b1[1] = " << f[1].getB1() << "\n"
                << "b2[1] = " << f[1].getB2() << "\n"
                ;

            std::cout << os.str();
        }
    }

}
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

namespace Dsp {

    namespace Elliptic {

        // shit ton of math in here

        // approximation to complete elliptic integral of the first kind.
        // fast convergence, peak error less than 2e-16.
        double Solver::ellipticK(double k)
        {
            double m = k * k;
            double a = 1;
            double b = sqrt(1 - m);
            double c = a - b;
            double co;
            do
            {
                co = c;
                c = (a - b) / 2;
                double ao = (a + b) / 2;
                b = sqrt(a * b);
                a = ao;
            } while (c < co);

            return doublePi / (a + a);
        }

        //------------------------------------------------------------------------------

        AnalogLowPass::AnalogLowPass()
            : m_numPoles(-1)
        {
            setNormal(0, 1);
        }

        void AnalogLowPass::design(int numPoles,
            double rippleDb,
            double rolloff)
        {
            if (m_numPoles != numPoles ||
                m_rippleDb != rippleDb ||
                m_rolloff != rolloff)
            {
                m_numPoles = numPoles;
                m_rippleDb = rippleDb;
                m_rolloff = rolloff;

                reset();

                // calculate
                //const double ep = rippleDb; // passband ripple

                const int n = numPoles;

                double e2 = pow(10., rippleDb / 10) - 1;
                //double xi = rolloff + 1;
                double xi = 5 * exp(rolloff - 1) + 1;

                m_K = Solver::ellipticK(1 / xi);
                m_Kprime = Solver::ellipticK(sqrt(1 - 1 / (xi * xi)));

                int ni = ((n & 1) == 1) ? 0 : 1;
                int i;
                double f[100]; // HACK!!!
                for (i = 1; i <= n / 2; i++)
                {
                    double u = (2 * i - ni) * m_K / n;
                    double sn = calcsn(u);
                    sn *= 2 * doublePi / m_K;
                    f[i] = m_zeros[i - 1] = 1 / sn;
                }
                m_zeros[n / 2] = std::numeric_limits<double>::infinity();
                double fb = 1 / (2 * doublePi);
                m_nin = n % 2;
                m_n2 = n / 2;
                for (i = 1; i <= m_n2; i++)
                {
                    double x = f[m_n2 + 1 - i];
                    m_z1[i] = sqrt(1 - 1 / (x * x));
                }
                double ee = e2;//pow(10., rippleDb/20)-1;
                m_e = sqrt(ee);
                double fbb = fb * fb;
                m_m = m_nin + 2 * m_n2;
                m_em = 2 * (m_m / 2);
                double tp = 2 * doublePi;
                calcfz();
                calcqz();
                if (m_m > m_em)
                    m_c1[2 * m_m] = 0;
                for (i = 0; i <= 2 * m_m; i += 2)
                    m_a1[m_m - i / 2] = m_c1[i] + m_d1[i];
                double a0 = findfact(m_m);
                int r = 0;
                while (r < m_em / 2)
                {
                    r++;
                    m_p[r] /= 10;
                    m_q1[r] /= 100;
                    double d = 1 + m_p[r] + m_q1[r];
                    m_b1[r] = (1 + m_p[r] / 2) * fbb / d;
                    m_zf1[r] = fb / pow(d, .25);
                    m_zq1[r] = 1 / sqrt(fabs(2 * (1 - m_b1[r] / (m_zf1[r] * m_zf1[r]))));
                    m_zw1[r] = tp * m_zf1[r];

                    m_rootR[r] = -.5 * m_zw1[r] / m_zq1[r];
                    m_rootR[r + m_em / 2] = m_rootR[r];
                    m_rootI[r] = .5 * sqrt(fabs(m_zw1[r] * m_zw1[r] / (m_zq1[r] * m_zq1[r]) - 4 * m_zw1[r] * m_zw1[r]));
                    m_rootI[r + m_em / 2] = -m_rootI[r];

                    complex_t pole(
                        -.5 * m_zw1[r] / m_zq1[r],
                        .5 * sqrt(fabs(m_zw1[r] * m_zw1[r] / (m_zq1[r] * m_zq1[r]) - 4 * m_zw1[r] * m_zw1[r])));

                    complex_t zero(0, m_zeros[r - 1]);

                    addPoleZeroConjugatePairs(pole, zero);
                }

                if (a0 != 0)
                {
                    m_rootR[r + 1 + m_em / 2] = -sqrt(fbb / (.1 * a0 - 1)) * tp;
                    m_rootI[r + 1 + m_em / 2] = 0;

                    add(-sqrt(fbb / (.1 * a0 - 1)) * tp, infinity());
                }

                setNormal(0, (numPoles & 1) ? 1. : pow(10., -rippleDb / 20.0));
            }
        }

        // generate the product of (z+s1[i]) for i = 1 .. sn and store it in b1[]
        // (i.e. f[z] = b1[0] + b1[1] z + b1[2] z^2 + ... b1[sn] z^sn)
        void AnalogLowPass::prodpoly(int sn)
        {
            m_b1[0] = m_s1[1];
            m_b1[1] = 1;
            int i, j;
            for (j = 2; j <= sn; j++)
            {
                m_a1[0] = m_s1[j] * m_b1[0];
                for (i = 1; i <= j - 1; i++)
                    m_a1[i] = m_b1[i - 1] + m_s1[j] * m_b1[i];
                for (i = 0; i != j; i++)
                    m_b1[i] = m_a1[i];
                m_b1[j] = 1;
            }
        }

        // determine f(z)^2
        void AnalogLowPass::calcfz2(int i)
        {
            int ji = 0;
            int jf = 0;
            if (i < m_em + 2)
            {
                ji = 0;
                jf = i;
            }
            if (i > m_em)
            {
                ji = i - m_em;
                jf = m_em;
            }
            m_c1[i] = 0;
            int j;
            for (j = ji; j <= jf; j += 2)
                m_c1[i] += m_a1[j] * (m_a1[i - j] * pow(10., m_m - i / 2));
        }

        // calculate f(z)
        void AnalogLowPass::calcfz(void)
        {
            int i = 1;
            if (m_nin == 1)
                m_s1[i++] = 1;
            for (; i <= m_nin + m_n2; i++)
                m_s1[i] = m_s1[i + m_n2] = m_z1[i - m_nin];
            prodpoly(m_nin + 2 * m_n2);
            for (i = 0; i <= m_em; i += 2)
                m_a1[i] = m_e * m_b1[i];
            for (i = 0; i <= 2 * m_em; i += 2)
                calcfz2(i);
        }

        // determine q(z)
        void AnalogLowPass::calcqz(void)
        {
            int i;
            for (i = 1; i <= m_nin; i++)
                m_s1[i] = -10;
            for (; i <= m_nin + m_n2; i++)
                m_s1[i] = -10 * m_z1[i - m_nin] * m_z1[i - m_nin];
            for (; i <= m_nin + 2 * m_n2; i++)
                m_s1[i] = m_s1[i - m_n2];
            prodpoly(m_m);
            int dd = ((m_nin & 1) == 1) ? -1 : 1;
            for (i = 0; i <= 2 * m_m; i += 2)
                m_d1[i] = dd * m_b1[i / 2];
        }

        // compute factors
        double AnalogLowPass::findfact(int t)
        {
            int i;
            double a = 0;
            for (i = 1; i <= t; i++)
                m_a1[i] /= m_a1[0];
            m_a1[0] = m_b1[0] = m_c1[0] = 1;
            int i1 = 0;
            for (;;)
            {
                if (t <= 2)
                    break;
                double p0 = 0, q0 = 0;
                i1++;
                for (;;)
                {
                    m_b1[1] = m_a1[1] - p0;
                    m_c1[1] = m_b1[1] - p0;
                    for (i = 2; i <= t; i++)
                        m_b1[i] = m_a1[i] - p0 * m_b1[i - 1] - q0 * m_b1[i - 2];
                    for (i = 2; i < t; i++)
                        m_c1[i] = m_b1[i] - p0 * m_c1[i - 1] - q0 * m_c1[i - 2];
                    int x1 = t - 1;
                    int x2 = t - 2;
                    int x3 = t - 3;
                    double x4 = m_c1[x2] * m_c1[x2] + m_c1[x3] * (m_b1[x1] - m_c1[x1]);
                    if (x4 == 0)
                        x4 = 1e-3;
                    double ddp = (m_b1[x1] * m_c1[x2] - m_b1[t] * m_c1[x3]) / x4;
                    p0 += ddp;
                    double dq = (m_b1[t] * m_c1[x2] - m_b1[x1] * (m_c1[x1] - m_b1[x1])) / x4;
                    q0 += dq;
                    if (fabs(ddp + dq) < 1e-6)
                        break;
                }
                m_p[i1] = p0;
                m_q1[i1] = q0;
                m_a1[1] = m_a1[1] - p0;
                t -= 2;
                for (i = 2; i <= t; i++)
                    m_a1[i] -= p0 * m_a1[i - 1] + q0 * m_a1[i - 2];
                if (t <= 2)
                    break;
            }

            if (t == 2)
            {
                i1++;
                m_p[i1] = m_a1[1];
                m_q1[i1] = m_a1[2];
            }
            if (t == 1)
                a = -m_a1[1];

            return a;
        }

        double AnalogLowPass::calcsn(double u)
        {
            double sn = 0;
            int j;
            // q = modular constant
            double q = exp(-doublePi * m_Kprime / m_K);
            double v = doublePi * .5 * u / m_K;
            for (j = 0; ; j++)
            {
                double w = pow(q, j + .5);
                sn += w * sin((2 * j + 1) * v) / (1 - w * w);
                if (w < 1e-7)
                    break;
            }
            return sn;
        }

        //------------------------------------------------------------------------------

        void LowPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double rippleDb,
            double rolloff)
        {
            m_analogProto.design(order, rippleDb, rolloff);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            double rippleDb,
            double rolloff)
        {
            m_analogProto.design(order, rippleDb, rolloff);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandPassBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double rippleDb,
            double rolloff)
        {
            m_analogProto.design(order, rippleDb, rolloff);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandStopBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            double rippleDb,
            double rolloff)
        {
            m_analogProto.design(order, rippleDb, rolloff);

            BandStopTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

    }

}
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

namespace Dsp {

    Params Filter::getDefaultParams() const
    {
        Params params;

        params.clear();

        for (int i = 0; i < getNumParams(); ++i)
            params[i] = getParamInfo(i).getDefaultValue();

        return params;
    }

    Filter::~Filter()
    {
    }

    int Filter::findParamId(int paramId)
    {
        int index = -1;

        for (int i = getNumParams(); --i >= 0;)
        {
            if (getParamInfo(i).getId() == paramId)
            {
                index = i;
                break;
            }
        }

        return index;
    }

    void Filter::setParamById(int paramId, double nativeValue)
    {
        for (int i = getNumParams(); --i >= 0;)
        {
            if (getParamInfo(i).getId() == paramId)
            {
                setParam(i, nativeValue);
                return;
            }
        }

        assert(0);
    }

    void Filter::copyParamsFrom(Dsp::Filter const* other)
    {
        // first, set reasonable defaults
        m_params = getDefaultParams();

        if (other)
        {
            // now loop
            for (int i = 0; i < getNumParams(); ++i)
            {
                const ParamInfo& paramInfo = getParamInfo(i);

                // find a match
                for (int j = 0; j < other->getNumParams(); ++j)
                {
                    const ParamInfo& otherParamInfo = other->getParamInfo(j);

                    if (paramInfo.getId() == otherParamInfo.getId())
                    {
                        // match!
                        m_params[i] = paramInfo.clamp(other->getParam(j));
                        break;
                    }
                }
            }
        }

        doSetParams(m_params);
    }

}
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

#include <sstream>
#include <iostream>
#include <iomanip>

namespace Dsp {

    namespace Legendre {

        static inline double m_sqrt2()
        {
            return 1.41421356237309504880;
        }

        //  Optimum 'L' Filter algorithm.
        //  (C) 2004, C. Bond.
        //
        //  Based on discussion in Kuo, "Network Analysis and Synthesis",
        //  pp. 379-383. Original method due to A.Papoulis."On Monotonic
        //  Response Filters", Proc. IRE, 47, Feb. 1959.
        //
        //  Rewritten by Vinnie Falco to change the way temporary
        //  storage is allocated
        //

        //
        //  This routine calculates the coefficients of the Legendre polynomial
        //  of the 1st kind. It uses a recursion relation. The first few polynomials
        //  are hard coded and the rest are found by recursion.
        //
        //  (n+1)Pn+1 = (2n+1)xPn - nPn-1 	Recursion relation.
        //
        void PolynomialFinderBase::legendre(double* p, int n)
        {
            int i, j;

            if (n == 0) {
                p[0] = 1.0;
                return;
            }
            if (n == 1) {
                p[0] = 0.0;
                p[1] = 1.0;
                return;
            }
            p[0] = -0.5;
            p[1] = 0.0;
            p[2] = 1.5;

            if (n == 2) return;

            for (i = 0; i <= n; i++) {
                m_aa[i] = m_bb[i] = 0.0;
            }
            m_bb[1] = 1.0;

            for (i = 3; i <= n; i++) {
                for (j = 0; j <= i; j++) {
                    m_aa[j] = m_bb[j];
                    m_bb[j] = p[j];
                    p[j] = 0.0;
                }
                for (j = i - 2; j >= 0; j -= 2) {
                    p[j] -= (i - 1) * m_aa[j] / i;
                }
                for (j = i - 1; j >= 0; j -= 2) {
                    p[j + 1] += (2 * i - 1) * m_bb[j] / i;
                }
            }
        }

        //
        //
        //  In the following routine n = 2k + 1 for odd 'n' and n = 2k + 2 for
        //  even 'n'.
        //
        //
        //      n   k
        //      -----
        //      1   0
        //      2   0
        //      3   1
        //      4   1
        //      5   2
        //      6   2
        //

        void PolynomialFinderBase::solve(int n)
        {
            assert(n <= m_maxN);

            double c0, c1;
            int i, j, k;

            k = (n - 1) / 2;
            //
            //  form vector of 'a' constants
            //
            if (n & 1) {                // odd
                for (i = 0; i <= k; i++) {
                    m_a[i] = (2.0 * i + 1.0) / (m_sqrt2() * (k + 1.0));
                }
            }                           // even
            else {
                for (i = 0; i < k + 1; i++) {
                    m_a[i] = 0.0;
                }
                if (k & 1) {
                    for (i = 1; i <= k; i += 2) {
                        m_a[i] = (2 * i + 1) / sqrt(double((k + 1) * (k + 2)));
                    }
                }
                else {
                    for (i = 0; i <= k; i += 2) {
                        m_a[i] = (2 * i + 1) / sqrt(double((k + 1) * (k + 2)));
                    }
                }
            }
            for (i = 0; i <= n; i++) {
                m_s[i] = 0.0;
                m_w[i] = 0.0;
            }
            //
            // form s[] = sum of a[i]*P[i]
            //
            m_s[0] = m_a[0];
            m_s[1] = m_a[1];
            for (i = 2; i <= k; i++) {
                legendre(m_p, i);
                for (j = 0; j <= i; j++) {
                    m_s[j] += m_a[i] * m_p[j];
                }
            }
            //
            //  form v[] = square of s[]
            //
            for (i = 0; i <= 2 * k + 2; i++) {
                m_v[i] = 0.0;
            }
            for (i = 0; i <= k; i++) {
                for (j = 0; j <= k; j++) {
                    m_v[i + j] += m_s[i] * m_s[j];
                }
            }
            //
            //  modify integrand for even 'n'
            //
            m_v[2 * k + 1] = 0.0;
            if ((n & 1) == 0) {
                for (i = n; i >= 0; i--) {
                    m_v[i + 1] += m_v[i];
                }
            }
            //
            //  form integral of v[]
            //
            for (i = n + 1; i >= 0; i--) {
                m_v[i + 1] = m_v[i] / (double)(i + 1.0);
            }
            m_v[0] = 0.0;
            //
            // clear s[] for use in computing definite integral
            //
            for (i = 0; i < (n + 2); i++) {
                m_s[i] = 0.0;
            }
            m_s[0] = -1.0;
            m_s[1] = 2.0;
            //
            //  calculate definite integral
            //
            for (i = 1; i <= n; i++) {
                if (i > 1) {
                    c0 = -m_s[0];
                    for (j = 1; j < i + 1; j++) {
                        c1 = -m_s[j] + 2.0 * m_s[j - 1];
                        m_s[j - 1] = c0;
                        c0 = c1;
                    }
                    c1 = 2.0 * m_s[i];
                    m_s[i] = c0;
                    m_s[i + 1] = c1;
                }
                for (j = i; j > 0; j--) {
                    m_w[j] += (m_v[i] * m_s[j]);
                }
            }
            if ((n & 1) == 0) m_w[1] = 0.0;
        }

        //------------------------------------------------------------------------------

        AnalogLowPass::AnalogLowPass()
            : m_numPoles(-1)
        {
            setNormal(0, 1);
        }

        void AnalogLowPass::design(int numPoles,
            WorkspaceBase* w)
        {
            if (m_numPoles != numPoles)
            {
                m_numPoles = numPoles;

                reset();

                PolynomialFinderBase& poly(w->poly);
                RootFinderBase& poles(w->roots);

                poly.solve(numPoles);
                int degree = numPoles * 2;

                poles.coef()[0] = 1 + poly.coef()[0];
                poles.coef()[1] = 0;
                for (int i = 1; i <= degree; ++i)
                {
                    poles.coef()[2 * i] = poly.coef()[i] * ((i & 1) ? -1 : 1);
                    poles.coef()[2 * i + 1] = 0;
                }
                poles.solve(degree);

                int j = 0;
                for (int i = 0; i < degree; ++i)
                    if (poles.root()[i].real() <= 0)
                        poles.root()[j++] = poles.root()[i];
                // sort descending imag() and cut degree in half
                poles.sort(degree / 2);

                const int pairs = numPoles / 2;
                for (int i = 0; i < pairs; ++i)
                {
                    complex_t c = poles.root()[i];
                    addPoleZeroConjugatePairs(c, infinity());
                }

                if (numPoles & 1)
                    add(poles.root()[pairs].real(), infinity());
            }
        }

        //------------------------------------------------------------------------------

        void LowPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            LowPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void HighPassBase::setup(int order,
            double sampleRate,
            double cutoffFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            HighPassTransform(cutoffFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandPassBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            BandPassTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

        void BandStopBase::setup(int order,
            double sampleRate,
            double centerFrequency,
            double widthFrequency,
            WorkspaceBase* w)
        {
            m_analogProto.design(order, w);

            BandStopTransform(centerFrequency / sampleRate,
                widthFrequency / sampleRate,
                m_digitalProto,
                m_analogProto);

            Cascade::setLayout(m_digitalProto);
        }

    }

}
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
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace Dsp {

    ParamInfo::ParamInfo()
    {
        throw std::logic_error("invalid usage of ParamInfo");
    }

    double ParamInfo::clamp(double nativeValue) const
    {
        const double minValue = toNativeValue(0);
        const double maxValue = toNativeValue(1);
        if (nativeValue < minValue)
            nativeValue = minValue;
        else if (nativeValue > maxValue)
            nativeValue = maxValue;
        return nativeValue;
    }

    //------------------------------------------------------------------------------

    double ParamInfo::Int_toControlValue(double nativeValue) const
    {
        return (nativeValue - m_arg1) / (m_arg2 - m_arg1);
    }

    double ParamInfo::Int_toNativeValue(double controlValue) const
    {
        return std::floor(m_arg1 + controlValue * (m_arg2 - m_arg1) + 0.5);
    }

    //------------------------------------------------------------------------------

    double ParamInfo::Real_toControlValue(double nativeValue) const
    {
        return (nativeValue - m_arg1) / (m_arg2 - m_arg1);
    }

    double ParamInfo::Real_toNativeValue(double controlValue) const
    {
        return m_arg1 + controlValue * (m_arg2 - m_arg1);
    }

    //------------------------------------------------------------------------------

    double ParamInfo::Log_toControlValue(double nativeValue) const
    {
        const double base = 1.5;
        double l0 = log(m_arg1) / log(base);
        double l1 = log(m_arg2) / log(base);
        return (log(nativeValue) / log(base) - l0) / (l1 - l0);
    }

    double ParamInfo::Log_toNativeValue(double controlValue) const
    {
        const double base = 1.5;
        double l0 = log(m_arg1) / log(base);
        double l1 = log(m_arg2) / log(base);
        return pow(base, l0 + controlValue * (l1 - l0));
    }

    //------------------------------------------------------------------------------

    double ParamInfo::Pow2_toControlValue(double nativeValue) const
    {
        return ((log(nativeValue) / log(2.)) - m_arg1) / (m_arg2 - m_arg1);
    }

    double ParamInfo::Pow2_toNativeValue(double controlValue) const
    {
        return pow(2., (controlValue * (m_arg2 - m_arg1)) + m_arg1);
    }

    //------------------------------------------------------------------------------

    std::string ParamInfo::Int_toString(double nativeValue) const
    {
        std::ostringstream os;
        os << int(nativeValue);
        return os.str();
    }

    std::string ParamInfo::Hz_toString(double nativeValue) const
    {
        std::ostringstream os;
        os << int(nativeValue) << " Hz";
        return os.str();
    }

    std::string ParamInfo::Real_toString(double nativeValue) const
    {
        std::ostringstream os;
        os << std::fixed << std::setprecision(3) << nativeValue;
        return os.str();
    }

    std::string ParamInfo::Db_toString(double nativeValue) const
    {
        const double af = fabs(nativeValue);
        int prec;
        if (af < 1)
            prec = 3;
        else if (af < 10)
            prec = 2;
        else
            prec = 1;
        std::ostringstream os;
        os << std::fixed << std::setprecision(prec) << nativeValue << " dB";
        return os.str();
    }

    //------------------------------------------------------------------------------

    ParamInfo ParamInfo::defaultSampleRateParam()
    {
        return ParamInfo(idSampleRate, "Fs", "Sample Rate",
            11025, 192000, 44100,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Hz_toString);
    }

    ParamInfo ParamInfo::defaultCutoffFrequencyParam()
    {
        return ParamInfo(idFrequency, "Fc", "Cutoff Frequency",
            10, 22040, 2000,
            &ParamInfo::Log_toControlValue,
            &ParamInfo::Log_toNativeValue,
            &ParamInfo::Hz_toString);
    }

    ParamInfo ParamInfo::defaultCenterFrequencyParam()
    {
        return ParamInfo(idFrequency, "Fc", "Center Frequency",
            10, 22040, 2000,
            &ParamInfo::Log_toControlValue,
            &ParamInfo::Log_toNativeValue,
            &ParamInfo::Hz_toString);
    }

    ParamInfo ParamInfo::defaultQParam()
    {
        return ParamInfo(idQ, "Q", "Resonance",
            -4, 4, 1,
            &ParamInfo::Pow2_toControlValue,
            &ParamInfo::Pow2_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultBandwidthParam()
    {
        return ParamInfo(idBandwidth, "BW", "Bandwidth (Octaves)",
            -4, 4, 1,
            &ParamInfo::Pow2_toControlValue,
            &ParamInfo::Pow2_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultBandwidthHzParam()
    {
        return ParamInfo(idBandwidthHz, "BW", "Bandwidth (Hz)",
            10, 22040, 1720,
            &ParamInfo::Log_toControlValue,
            &ParamInfo::Log_toNativeValue,
            &ParamInfo::Hz_toString);
    }

    ParamInfo ParamInfo::defaultGainParam()
    {
        return ParamInfo(idGain, "Gain", "Gain",
            -24, 24, -6,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Db_toString);
    }

    ParamInfo ParamInfo::defaultSlopeParam()
    {
        return ParamInfo(idSlope, "Slope", "Slope",
            -2, 2, 1,
            &ParamInfo::Pow2_toControlValue,
            &ParamInfo::Pow2_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultRippleDbParam()
    {
        return ParamInfo(idRippleDb, "Ripple", "Ripple dB",
            0.001, 12, 0.01,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Db_toString);
    }

    ParamInfo ParamInfo::defaultStopDbParam()
    {
        return ParamInfo(idStopDb, "Stop", "Stopband dB",
            3, 60, 48,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Db_toString);
    }

    ParamInfo ParamInfo::defaultRolloffParam()
    {
        return ParamInfo(idRolloff, "W", "Transition Width",
            -16, 4, 0,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultPoleRhoParam()
    {
        return ParamInfo(idPoleRho, "Pd", "Pole Distance",
            0, 1, 0.5,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultPoleThetaParam()
    {
        return ParamInfo(idPoleTheta, "Pa", "Pole Angle",
            0, doublePi, doublePi / 2,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultZeroRhoParam()
    {
        return ParamInfo(idZeroRho, "Pd", "Zero Distance",
            0, 1, 0.5,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultZeroThetaParam()
    {
        return ParamInfo(idZeroTheta, "Pa", "Zero Angle",
            0, doublePi, doublePi / 2,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultPoleRealParam()
    {
        return ParamInfo(idPoleReal, "A1", "Pole Real",
            -1, 1, 0.25,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }

    ParamInfo ParamInfo::defaultZeroRealParam()
    {
        return ParamInfo(idZeroReal, "B1", "Zero Real",
            -1, 1, -0.25,
            &ParamInfo::Real_toControlValue,
            &ParamInfo::Real_toNativeValue,
            &ParamInfo::Real_toString);
    }
}

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

namespace Dsp {

    //------------------------------------------------------------------------------

    complex_t LowPassTransform::transform(complex_t c)
    {
        if (c == infinity())
            return complex_t(-1, 0);

        // frequency transform
        c = f * c;

        // bilinear low pass transform
        return (1. + c) / (1. - c);
    }

    LowPassTransform::LowPassTransform(double fc,
        LayoutBase& digital,
        LayoutBase const& analog)
    {
        digital.reset();

        // prewarp
        f = tan(doublePi * fc);

        const int numPoles = analog.getNumPoles();
        const int pairs = numPoles / 2;
        for (int i = 0; i < pairs; ++i)
        {
            const PoleZeroPair& pair = analog[i];
            digital.addPoleZeroConjugatePairs(transform(pair.poles.first),
                transform(pair.zeros.first));
        }

        if (numPoles & 1)
        {
            const PoleZeroPair& pair = analog[pairs];
            digital.add(transform(pair.poles.first),
                transform(pair.zeros.first));
        }

        digital.setNormal(analog.getNormalW(),
            analog.getNormalGain());
    }

    //------------------------------------------------------------------------------

    complex_t HighPassTransform::transform(complex_t c)
    {
        if (c == infinity())
            return complex_t(1, 0);

        // frequency transform
        c = f * c;

        // bilinear high pass transform
        return -(1. + c) / (1. - c);
    }

    HighPassTransform::HighPassTransform(double fc,
        LayoutBase& digital,
        LayoutBase const& analog)
    {
        digital.reset();

        // prewarp
        f = 1. / tan(doublePi * fc);

        const int numPoles = analog.getNumPoles();
        const int pairs = numPoles / 2;
        for (int i = 0; i < pairs; ++i)
        {
            const PoleZeroPair& pair = analog[i];
            digital.addPoleZeroConjugatePairs(transform(pair.poles.first),
                transform(pair.zeros.first));
        }

        if (numPoles & 1)
        {
            const PoleZeroPair& pair = analog[pairs];
            digital.add(transform(pair.poles.first),
                transform(pair.zeros.first));
        }

        digital.setNormal(doublePi - analog.getNormalW(),
            analog.getNormalGain());
    }

    //------------------------------------------------------------------------------

    BandPassTransform::BandPassTransform(double fc,
        double fw,
        LayoutBase& digital,
        LayoutBase const& analog)
    {
        // handle degenerate cases efficiently
        // THIS DOESNT WORK because the cascade states won't match
#if 0
        const double fw_2 = fw / 2;
        if (fc - fw_2 < 0)
        {
            LowPassTransform::transform(fc + fw_2, digital, analog);
        }
        else if (fc + fw_2 >= 0.5)
        {
            HighPassTransform::transform(fc - fw_2, digital, analog);
        }
        else
#endif

            digital.reset();

        const double ww = 2 * doublePi * fw;

        // pre-calcs
        wc2 = 2 * doublePi * fc - (ww / 2);
        wc = wc2 + ww;

        // what is this crap?
        if (wc2 < 1e-8)
            wc2 = 1e-8;
        if (wc > doublePi - 1e-8)
            wc = doublePi - 1e-8;

        a = cos((wc + wc2) * 0.5) /
            cos((wc - wc2) * 0.5);
        b = 1 / tan((wc - wc2) * 0.5);
        a2 = a * a;
        b2 = b * b;
        ab = a * b;
        ab_2 = 2 * ab;

        const int numPoles = analog.getNumPoles();
        const int pairs = numPoles / 2;
        for (int i = 0; i < pairs; ++i)
        {
            const PoleZeroPair& pair = analog[i];
            ComplexPair p1 = transform(pair.poles.first);
            ComplexPair z1 = transform(pair.zeros.first);

            //
            // Optimize out the calculations for conjugates for Release builds
            //
#ifndef NDEBUG
            ComplexPair p2 = transform(pair.poles.second);
            ComplexPair z2 = transform(pair.zeros.second);
            assert(p2.first == std::conj(p1.first));
            assert(p2.second == std::conj(p1.second));
#endif

            digital.addPoleZeroConjugatePairs(p1.first, z1.first);
            digital.addPoleZeroConjugatePairs(p1.second, z1.second);
        }

        if (numPoles & 1)
        {
            ComplexPair poles = transform(analog[pairs].poles.first);
            ComplexPair zeros = transform(analog[pairs].zeros.first);

            digital.add(poles, zeros);
        }

        double wn = analog.getNormalW();
        digital.setNormal(2 * atan(sqrt(tan((wc + wn) * 0.5) * tan((wc2 + wn) * 0.5))),
            analog.getNormalGain());
    }

    ComplexPair BandPassTransform::transform(complex_t c)
    {
        if (c == infinity())
            return ComplexPair(-1, 1);

        c = (1. + c) / (1. - c); // bilinear

        complex_t v = 0;
        v = addmul(v, 4 * (b2 * (a2 - 1) + 1), c);
        v += 8 * (b2 * (a2 - 1) - 1);
        v *= c;
        v += 4 * (b2 * (a2 - 1) + 1);
        v = std::sqrt(v);

        complex_t u = -v;
        u = addmul(u, ab_2, c);
        u += ab_2;

        v = addmul(v, ab_2, c);
        v += ab_2;

        complex_t d = 0;
        d = addmul(d, 2 * (b - 1), c) + 2 * (1 + b);

        return ComplexPair(u / d, v / d);
    }

    //------------------------------------------------------------------------------

    BandStopTransform::BandStopTransform(double fc,
        double fw,
        LayoutBase& digital,
        LayoutBase const& analog)
    {
        digital.reset();

        const double ww = 2 * doublePi * fw;

        wc2 = 2 * doublePi * fc - (ww / 2);
        wc = wc2 + ww;

        // this is crap
        if (wc2 < 1e-8)
            wc2 = 1e-8;
        if (wc > doublePi - 1e-8)
            wc = doublePi - 1e-8;

        a = cos((wc + wc2) * .5) /
            cos((wc - wc2) * .5);
        b = tan((wc - wc2) * .5);
        a2 = a * a;
        b2 = b * b;

        const int numPoles = analog.getNumPoles();
        const int pairs = numPoles / 2;
        for (int i = 0; i < pairs; ++i)
        {
            const PoleZeroPair& pair = analog[i];
            ComplexPair p = transform(pair.poles.first);
            ComplexPair z = transform(pair.zeros.first);

            //
            // Optimize out the calculations for conjugates for Release builds
            //
#ifdef NDEBUG
    // trick to get the conjugate
            if (z.second == z.first)
                z.second = std::conj(z.first);

#else
    // Do the full calculation to verify correctness
            ComplexPair pc = transform(analog[i].poles.second);
            ComplexPair zc = transform(analog[i].zeros.second);

            // get the conjugates into pc and zc
            if (zc.first == z.first)
                std::swap(zc.first, zc.second);

            assert(pc.first == std::conj(p.first));
            assert(pc.second == std::conj(p.second));
            assert(zc.first == std::conj(z.first));
            assert(zc.second == std::conj(z.second));

#endif

            digital.addPoleZeroConjugatePairs(p.first, z.first);
            digital.addPoleZeroConjugatePairs(p.second, z.second);
        }

        if (numPoles & 1)
        {
            ComplexPair poles = transform(analog[pairs].poles.first);
            ComplexPair zeros = transform(analog[pairs].zeros.first);

            digital.add(poles, zeros);
        }

        if (fc < 0.25)
            digital.setNormal(doublePi, analog.getNormalGain());
        else
            digital.setNormal(0, analog.getNormalGain());
    }

    ComplexPair BandStopTransform::transform(complex_t c)
    {
        if (c == infinity())
            c = -1;
        else
            c = (1. + c) / (1. - c); // bilinear

        complex_t u(0);
        u = addmul(u, 4 * (b2 + a2 - 1), c);
        u += 8 * (b2 - a2 + 1);
        u *= c;
        u += 4 * (a2 + b2 - 1);
        u = std::sqrt(u);

        complex_t v = u * -.5;
        v += a;
        v = addmul(v, -a, c);

        u *= .5;
        u += a;
        u = addmul(u, -a, c);

        complex_t d(b + 1);
        d = addmul(d, b - 1, c);

        return ComplexPair(u / d, v / d);
    }

}
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


namespace Dsp {

    namespace RBJ {

        void LowPass::setup(double sampleRate,
            double cutoffFrequency,
            double q)
        {
            double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / (2 * q);
            double b0 = (1 - cs) / 2;
            double b1 = 1 - cs;
            double b2 = (1 - cs) / 2;
            double a0 = 1 + AL;
            double a1 = -2 * cs;
            double a2 = 1 - AL;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void HighPass::setup(double sampleRate,
            double cutoffFrequency,
            double q)
        {
            double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / (2 * q);
            double b0 = (1 + cs) / 2;
            double b1 = -(1 + cs);
            double b2 = (1 + cs) / 2;
            double a0 = 1 + AL;
            double a1 = -2 * cs;
            double a2 = 1 - AL;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void BandPass1::setup(double sampleRate,
            double centerFrequency,
            double bandWidth)
        {
            double w0 = 2 * doublePi * centerFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / (2 * bandWidth);
            double b0 = bandWidth * AL;// sn / 2;
            double b1 = 0;
            double b2 = -bandWidth * AL;//-sn / 2;
            double a0 = 1 + AL;
            double a1 = -2 * cs;
            double a2 = 1 - AL;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void BandPass2::setup(double sampleRate,
            double centerFrequency,
            double bandWidth)
        {
            double w0 = 2 * doublePi * centerFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / (2 * bandWidth);
            double b0 = AL;
            double b1 = 0;
            double b2 = -AL;
            double a0 = 1 + AL;
            double a1 = -2 * cs;
            double a2 = 1 - AL;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void BandStop::setup(double sampleRate,
            double centerFrequency,
            double bandWidth)
        {
            double w0 = 2 * doublePi * centerFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / (2 * bandWidth);
            double b0 = 1;
            double b1 = -2 * cs;
            double b2 = 1;
            double a0 = 1 + AL;
            double a1 = -2 * cs;
            double a2 = 1 - AL;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void LowShelf::setup(double sampleRate,
            double cutoffFrequency,
            double gainDb,
            double shelfSlope)
        {
            double A = pow(10, gainDb / 40);
            double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / 2 * ::std::sqrt((A + 1 / A) * (1 / shelfSlope - 1) + 2);
            double sq = 2 * sqrt(A) * AL;
            double b0 = A * ((A + 1) - (A - 1) * cs + sq);
            double b1 = 2 * A * ((A - 1) - (A + 1) * cs);
            double b2 = A * ((A + 1) - (A - 1) * cs - sq);
            double a0 = (A + 1) + (A - 1) * cs + sq;
            double a1 = -2 * ((A - 1) + (A + 1) * cs);
            double a2 = (A + 1) + (A - 1) * cs - sq;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void HighShelf::setup(double sampleRate,
            double cutoffFrequency,
            double gainDb,
            double shelfSlope)
        {
            double A = pow(10, gainDb / 40);
            double w0 = 2 * doublePi * cutoffFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / 2 * ::std::sqrt((A + 1 / A) * (1 / shelfSlope - 1) + 2);
            double sq = 2 * sqrt(A) * AL;
            double b0 = A * ((A + 1) + (A - 1) * cs + sq);
            double b1 = -2 * A * ((A - 1) + (A + 1) * cs);
            double b2 = A * ((A + 1) + (A - 1) * cs - sq);
            double a0 = (A + 1) - (A - 1) * cs + sq;
            double a1 = 2 * ((A - 1) - (A + 1) * cs);
            double a2 = (A + 1) - (A - 1) * cs - sq;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void BandShelf::setup(double sampleRate,
            double centerFrequency,
            double gainDb,
            double bandWidth)
        {
            double A = pow(10, gainDb / 40);
            double w0 = 2 * doublePi * centerFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn * sinh(doubleLn2 / 2 * bandWidth * w0 / sn);
            assert(!Dsp::is_nan(AL));
            double b0 = 1 + AL * A;
            double b1 = -2 * cs;
            double b2 = 1 - AL * A;
            double a0 = 1 + AL / A;
            double a1 = -2 * cs;
            double a2 = 1 - AL / A;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

        void AllPass::setup(double sampleRate,
            double phaseFrequency,
            double q)
        {
            double w0 = 2 * doublePi * phaseFrequency / sampleRate;
            double cs = cos(w0);
            double sn = sin(w0);
            double AL = sn / (2 * q);
            double b0 = 1 - AL;
            double b1 = -2 * cs;
            double b2 = 1 + AL;
            double a0 = 1 + AL;
            double a1 = -2 * cs;
            double a2 = 1 - AL;
            setCoefficients(a0, a1, a2, b0, b1, b2);
        }

    }

}
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
#include <stdexcept>

namespace Dsp {

    void RootFinderBase::solve(int degree,
        bool polish,
        bool doSort)
    {
        assert(degree <= m_maxdegree);

        const double EPS = 1.0e-30;

        int its;
        complex_t x, b, c;

        int m = degree;

        // copy coefficients
        for (int j = 0; j <= m; ++j)
            m_ad[j] = m_a[j];

        // for each root
        for (int j = m - 1; j >= 0; --j)
        {
            // initial guess at 0
            x = 0.0;
            laguerre(j + 1, m_ad, x, its);

            if (fabs(std::imag(x)) <= 2.0 * EPS * fabs(std::real(x)))
                x = complex_t(std::real(x), 0.0);

            m_root[j] = x;

            // deflate
            b = m_ad[j + 1];
            for (int jj = j; jj >= 0; --jj)
            {
                c = m_ad[jj];
                m_ad[jj] = b;
                b = x * b + c;
            }
        }

        if (polish)
            for (int j = 0; j < m; ++j)
                laguerre(degree, m_a, m_root[j], its);

        if (doSort)
            sort(degree);
    }

    void RootFinderBase::sort(int degree)
    {
        for (int j = 1; j < degree; ++j)
        {
            complex_t x = m_root[j];

            int i;
            for (i = j - 1; i >= 0; --i)
            {
                if (m_root[i].imag() >= x.imag())
                    break;

                m_root[i + 1] = m_root[i];
            }

            m_root[i + 1] = x;
        }
    }

    //------------------------------------------------------------------------------

    void RootFinderBase::laguerre(int degree,
        complex_t a[],
        complex_t& x,
        int& its)
    {
        const int MR = 8, MT = 10, MAXIT = MT * MR;
        const double EPS = std::numeric_limits<double>::epsilon();

        static const double frac[MR + 1] =
        { 0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0 };

        complex_t dx, x1, b, d, f, g, h, sq, gp, gm, g2;

        int m = degree;
        for (int iter = 1; iter <= MAXIT; ++iter)
        {
            its = iter;
            b = a[m];
            double err = std::abs(b);
            d = f = 0.0;
            double abx = std::abs(x);
            for (int j = m - 1; j >= 0; --j)
            {
                f = x * f + d;
                d = x * d + b;
                b = x * b + a[j];
                err = std::abs(b) + abx * err;
            }
            err *= EPS;
            if (std::abs(b) <= err)
                return;
            g = d / b;
            g2 = g * g;
            h = g2 - 2.0 * f / b;

            sq = sqrt(double(m - 1) * (double(m) * h - g2));
            gp = g + sq;
            gm = g - sq;

            double abp = std::abs(gp);
            double abm = std::abs(gm);
            if (abp < abm)
                gp = gm;
            dx = std::max(abp, abm) > 0.0 ? double(m) / gp : std::polar(1 + abx, double(iter));
            x1 = x - dx;
            if (x == x1)
                return;
            if (iter % MT != 0)
                x = x1;
            else
                x -= frac[iter / MT] * dx;
        }

        throw std::logic_error("laguerre failed");
    }

    //------------------------------------------------------------------------------

    complex_t RootFinderBase::eval(int degree,
        const complex_t& x)
    {
        complex_t y;

        if (x != 0.)
        {
            for (int i = 0; i <= degree; ++i)
                y += m_a[i] * pow(x, double(i));
        }
        else
        {
            y = m_a[0];
        }

        return y;
    }


}
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



namespace Dsp {

    //------------------------------------------------------------------------------

}









