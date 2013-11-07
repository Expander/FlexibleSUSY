// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include <cmath>
#include "pv.hpp"

#ifdef ENABLE_LOOPTOOLS
#  include <clooptools.h>
#elif defined(ENABLE_FFLITE)
#  include "fflite.hpp"
#else
#  include "numerics.h"
#endif

namespace flexiblesusy {

namespace passarino_veltman {

using namespace std;

#ifdef ENABLE_LOOPTOOLS

namespace {

struct Initialize_looptools {
    Initialize_looptools() {
	ltini();
    }
    ~Initialize_looptools() {
	ltexi();
    }
} initialize_looptools;

}

complex<double> A0(double m2, double scl2)
{
    setmudim(scl2);
    return ::A0(m2);
}

complex<double> B0(double p2, double m2a, double m2b, double scl2)
{
    setmudim(scl2);
    return ::B0(p2, m2a, m2b);
}

complex<double> B1(double p2, double m2a, double m2b, double scl2)
{
    setmudim(scl2);
    return ::B1(p2, m2a, m2b);
}

complex<double> B00(double p2, double m2a, double m2b, double scl2)
{
    setmudim(scl2);
    return ::B00(p2, m2a, m2b);
}

complex<double> A0(complex<double> m2, double scl2)
{
    setmudim(scl2);
    return ::A0C(m2);
}

complex<double> B0
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2)
{
    setmudim(scl2);
    return ::B0C(p2, m2a, m2b);
}

complex<double> B1
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2)
{
    setmudim(scl2);
    return ::B1C(p2, m2a, m2b);
}

complex<double> B00
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2)
{
    setmudim(scl2);
    return ::B00C(p2, m2a, m2b);
}

#elif defined(ENABLE_FFLITE)

namespace {

const double eps = 1e-22;	// see src/include/ff.h in LoopTools

struct Initialize_looptools {
    Initialize_looptools() {
	ltini_();
    }
    ~Initialize_looptools() {
	ltexi_();
    }
} initialize_looptools;

}

complex<double> A0(double m2, double scl2)
{
    complex<double> ca0;
    int ier;
    ljffxa0_(ca0, 0, scl2, m2, ier);
    return ca0;
}

complex<double> B0(double p2, double m2a, double m2b, double scl2)
{
    // see src/B/Bcoeff.F in LoopTools
    if (fabs(p2) + fabs(m2a) + fabs(m2b) < eps) return 0.0;

    complex<double> cb0;
    int ier;
    ljffxb0_(cb0, 0, scl2, p2, m2a, m2b, ier);
    return cb0;
}

complex<double> B1(double p2, double m2a, double m2b, double scl2)
{
    // see src/B/Bcoeff.F in LoopTools
    if (fabs(p2) + fabs(m2a) + fabs(m2b) < eps) return 0.0;

    complex<double> cb1;
    int ier;
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    double piDpj[9];
    ljffdot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffxb1_(cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb1;
}

complex<double> B00(double p2, double m2a, double m2b, double scl2)
{
    complex<double> cb2i[2];
    int ier;
    complex<double> cb1 = B1(p2, m2a, m2b, scl2);
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    double piDpj[9];
    ljffdot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffxb2p_(cb2i, cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb2i[1];
}

complex<double> A0
(complex<double> m2, double scl2)
{
    complex<double> ca0;
    int ier;
    ljffca0_(ca0, 0, scl2, m2, ier);
    return ca0;
}

complex<double> B0
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2)
{
    complex<double> cb0;
    int ier;
    ljffcb0_(cb0, 0, scl2, p2, m2a, m2b, ier);
    return cb0;
}

complex<double> B1
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2)
{
    complex<double> cb1;
    int ier;
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    complex<double> piDpj[9];
    ljffcot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffcb1_(cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb1;
}

complex<double> B00
(complex<double> p2, complex<double> m2a, complex<double> m2b, double scl2)
{
    complex<double> cb2i[2];
    int ier;
    complex<double> cb1 = B1(p2, m2a, m2b, scl2);
    complex<double> cb0 = B0(p2, m2a, m2b, scl2);
    complex<double> ca0i[2] = { A0(m2a, scl2), A0(m2b, scl2) };
    complex<double> piDpj[9];
    ljffcot2_(piDpj, p2, m2a, m2b, m2a-p2, m2b-p2, m2a-m2b, ier);
    ljffcb2p_(cb2i, cb1, cb0, ca0i, p2, m2a, m2b, piDpj, ier);
    return cb2i[1];
}

#endif // defined(ENABLE_FFLITE)

#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)

// CHECK: are the following correct complexifications of B22, H0, F0, G0?

complex<double> B22(double p2, double m2a, double m2b, double scl2)
{
    return B00(p2, m2a, m2b, scl2) - A0(m2a, scl2)/4.0 - A0(m2b, scl2)/4.0;
}

complex<double> H0(double p2, double m2a, double m2b, double scl2)
{
    return 4.0*B00(p2, m2a, m2b, scl2) + G0(p2, m2a, m2b, scl2);
}

complex<double> F0(double p2, double m2a, double m2b, double scl2)
{
    return A0(m2a, scl2) - 2.0*A0(m2b, scl2)
	   - (2*p2 + 2*m2a - m2b) * B0(p2, m2a, m2b, scl2);
}

complex<double> G0(double p2, double m2a, double m2b, double scl2)
{
    return (p2 - m2a - m2b) * B0(p2, m2a, m2b, scl2)
	   - A0(m2a, scl2) - A0(m2b, scl2);
}

#endif // defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)

double ReA0(double m2, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return A0(m2, scl2).real();
#else
    return a0(sqrt(m2), sqrt(scl2));
#endif
}

double ReB0(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B0(p2, m2a, m2b, scl2).real();
#else
    return b0(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
#endif
}

double ReB1(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B1(p2, m2a, m2b, scl2).real();
#else
    return -b1(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
#endif
}

double ReB00(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B00(p2, m2a, m2b, scl2).real();
#else
    return b22(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
#endif
}

double ReB22(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return B22(p2, m2a, m2b, scl2).real();
#else
    return ReB00(p2, m2a, m2b, scl2) - ReA0(m2a, scl2)/4 - ReA0(m2b, scl2)/4;
#endif
}

double ReH0(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return H0(p2, m2a, m2b, scl2).real();
#else
    return 4*ReB00(p2, m2a, m2b, scl2) + ReG0(p2, m2a, m2b, scl2);
#endif
}

double ReF0(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return F0(p2, m2a, m2b, scl2).real();
#else
    return ReA0(m2a, scl2) - 2*ReA0(m2b, scl2)
	   - (2*p2 + 2*m2a - m2b) * ReB0(p2, m2a, m2b, scl2);
#endif
}

double ReG0(double p2, double m2a, double m2b, double scl2)
{
#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)
    return G0(p2, m2a, m2b, scl2).real();
#else
    return (p2 - m2a - m2b) * ReB0(p2, m2a, m2b, scl2)
	   - ReA0(m2a, scl2) - ReA0(m2b, scl2);
#endif
}

} // passarino_veltman

} // flexiblesusy
