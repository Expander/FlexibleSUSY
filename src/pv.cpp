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

#ifdef ENABLE_LOOPTOOLS
#  include <clooptools.h>
#else
#  include "numerics.h"
#endif

#include <cmath>
#include "pv.hpp"

namespace flexiblesusy {

namespace passarino_veltman {

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

#endif // ENABLE_LOOPTOOLS

using namespace std;

double ReA0(double m2, double scl2)
{
#ifdef ENABLE_LOOPTOOLS
    setmudim(scl2);
    return ::A0(m2).real();
#else
    return a0(sqrt(m2), sqrt(scl2));
#endif
}

double ReB0(double p2, double m21, double m22, double scl2)
{
#ifdef ENABLE_LOOPTOOLS
    setmudim(scl2);
    return ::B0(p2, m21, m22).real();
#else
    return b0(sqrt(p2), sqrt(m21), sqrt(m22), sqrt(scl2));
#endif
}

double ReB1(double p2, double m21, double m22, double scl2)
{
#ifdef ENABLE_LOOPTOOLS
    setmudim(scl2);
    return ::B1(p2, m21, m22).real();
#else
    return -b1(sqrt(p2), sqrt(m21), sqrt(m22), sqrt(scl2));
#endif
}

double ReB00(double p2, double m21, double m22, double scl2)
{
#ifdef ENABLE_LOOPTOOLS
    setmudim(scl2);
    return ::B00(p2, m21, m22).real();
#else
    return b22(p2, m21, m22, sqrt(scl2));
#endif
}

double ReB22(double p2, double m21, double m22, double scl2)
{
    return ReB00(p2, m21, m22, scl2) - ReA0(m21, scl2)/4 - ReA0(m22, scl2)/4;
}

double ReH0(double p2, double m21, double m22, double scl2)
{
    return 4*ReB00(p2, m21, m22, scl2) + ReG0(p2, m21, m22, scl2);
}

double ReF0(double p2, double m21, double m22, double scl2)
{
    return ReA0(m21, scl2) - 2*ReA0(m22, scl2)
	   - (2*p2 + 2*m21 - m22) * ReB0(p2, m21, m22, scl2);
}

double ReG0(double p2, double m21, double m22, double scl2)
{
    return (p2 - m21 - m22) * ReB0(p2, m21, m22, scl2)
	   - ReA0(m21, scl2) - ReA0(m22, scl2);
}

} // passarino_veltman

} // flexiblesusy
