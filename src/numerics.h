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

/** \file numerics.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/

   \brief loop functions
*/

#ifndef NUMERICS_H
#define NUMERICS_H

#include <cmath>
#include "def.h"
#include "utils.h"

namespace softsusy {

/// Passarino-Veltman function definition
double b0(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b1(double p, double m1, double m2, double q);
/// Passarino-Veltman function definition
double b22(double p,  double m1, double m2, double q);
/// Passarino-Veltman function definition
double c0(double m1, double m2, double m3);
/// Passarino-Veltman function definition
double d27(double m1, double m2, double m3, double m4);
/// Passarino-Veltman function definition
double d0(double m1, double m2, double m3, double m4);

// inlined PV functions
inline double a0(double m, double q) {
   using std::fabs;
   using std::log;
  if (fabs(m) < softsusy::EPSTOL) return 0.;
  return sqr(m) * (1.0 - 2. * log(abs(m / q)));
}

inline double ffn(double p, double m1, double m2, double q) {
  return a0(m1, q) - 2.0 * a0(m2, q) -
    (2.0 * sqr(p) + 2.0 * sqr(m1) - sqr(m2)) *
    b0(p, m1, m2, q);
}

inline double gfn(double p, double m1, double m2, double q) {
  return (sqr(p) - sqr(m1) - sqr(m2)) * b0(p, m1, m2, q) - a0(m1, q)
    - a0(m2, q);
}

inline double hfn(double p, double m1, double m2, double q) {
  return 4.0 * b22(p, m1, m2, q) + gfn(p, m1, m2, q);
}

inline double b22bar(double p, double m1, double m2, double q) {
  return b22(p, m1, m2, q) - 0.25 * a0(m1, q) - 0.25 * a0(m2, q);
}

} // namespace softsusy

#endif
