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

#include "ffunctions.hpp"
#include "logger.hpp"
#include "numerics.h"
#include "numerics2.hpp"

#include <cmath>

namespace flexiblesusy {

namespace {
inline double sqr(double x)  { return x*x; }
inline double cube(double x) { return x*x*x; }
inline double quad(double x) { return sqr(x)*sqr(x); }
}

namespace gm2 {

double F1C(double x) {
   if(is_equal(x, 1.))
      ERROR("F1C: x must not be 1 !");

   return 2. / quad(1. - x) * (2. + 3. * x - 6. * sqr(x)
                           + cube(x) + 6. * x * log(x));
}

double F2C(double x) {
   if(is_equal(x, 1.))
      ERROR("F2C: x must not be 1 !");

   return 3. / (2. * cube(1. - x)) * (- 3. + 4. * x - sqr(x) - 2. * log(x));
}

double F3C(double x) {
   if(is_equal(x, 1.))
      ERROR("F3C: x must not be 1 !");

   return ( 4. / (141. * quad(1. - x)) * ((1. - x) * (151. * sqr(x) - 335. * x + 592.)
             + 6. * (21. * cube(x) - 108. * sqr(x) - 93. * x + 50.) * log(x)
             - 54. * x * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 108. * x * (sqr(x) - 2. * x + 12.) * dilog(1.- x)) );
}

double F4C(double x) {
   if(is_equal(x, 1.))
      ERROR("F4C: x must not be 1 !");

   return ( - 9. / (122. * cube(1. - x)) * (8. * (sqr(x) - 3. * x + 2.)
             +  (11. * sqr(x) - 40. * x + 5.) * log(x)
             - 2. * (sqr(x) - 2. * x - 2.) * sqr(log(x))
             - 4. * (sqr(x) - 2. * x + 9.) * dilog(1.- x)) );
}

double F1N(double x) {
   if(is_equal(x, 1.))
      ERROR("F1N: x must not be 1 !");

   return 2. / quad(1. - x) * (1. - 6. * x + 3. * sqr(x)
                          + 2. * cube(x) - 6. * sqr(x) * log(x));
}

double F2N(double x) {
   if(is_equal(x, 1.))
      ERROR("F2N: x must not be 1 !");

   return 3. / cube(1. - x) * (1. - sqr(x) + 2. * x * log(x));
}

double F3N(double x) {
   if(is_equal(x, 1.))
      ERROR("F3N: x must not be 1 !");

   return 4. / (105. * quad(1. - x)) * ((1. - x) * (- 97. * sqr(x) - 529. * x + 2.)
            + 6. * sqr(x) * (13. * x + 81.) * log(x)
            + 108. * x * (7. * x + 4.) * dilog(1. - x));
}

double F4N(double x) {
   if(is_equal(x, 1.))
      ERROR("F4N: x must not be 1 !");

   return - 2.25 / cube(1. - x) * ((x + 3.) * (x * log(x) + x - 1.)
                                  + (6. * x + 2.) * dilog(1. - x));
}

double Fa(double x, double y) {
   if(is_equal(x, y))
      ERROR("Fa: x must not be equal y!");

   return - (G3(x) - G3(y)) / (x - y);
}

double Fb(double x, double y) {
   if(is_equal(x, y))
      ERROR("Fb: x must not be equal y!");

   return - (G4(x) - G4(y)) / (x - y);
}

double G3(double x) {
   if(is_equal(x, 1.))
      ERROR("G3: x must not be 1 !");

   return 1. / (2. * cube(x - 1.)) * ((x - 1.) * (x - 3.) + 2. * log(x));
}

double G4(double x) {
   if(is_equal(x, 1.))
      ERROR("G4: x must not be 1 !");

   return 1. / (2. * cube(x - 1.)) * ((x - 1.) * (x + 1.) - 2. * x * log(x));
}

double Iabc(double a, double b, double c) {

   return ( (sqr(a * b) * log(sqr(a / b))
           + sqr(b * c) * log(sqr(b / c))
           + sqr(c * a) * log(sqr(c / a)))
           / ((sqr(a) - sqr(b)) * (sqr(b) - sqr(c)) * (sqr(a) - sqr(c))) );
}

double f_PS(double z) {
   double result = 0.;
   if(z < 0.25) {
      double y = sqrt(1. - 4. * z);
      result = 2. * z / y * (dilog(1. - 0.5 * (1. - y) / z) - dilog(1. - 0.5 * (1. + y) / z));
   } else {
      Complex y = sqrt(Complex(1. - 4. * z, 0.));
      Complex zc(z, 0.);
      result = real(2. * zc / y * (dilog(1. - 0.5 * (1. - y) / zc) - dilog(1. - 0.5 * (1. + y) / zc)));
   }

   return result;
}

double f_S(double z) {
   if(z < 0.)
      ERROR("f_S: z must not be negativ!");

   return (2. * z - 1.) * f_PS(z) - 2. * z * (2. + log(z));
}

double f_sferm(double z) {
   if(z < 0.)
      ERROR("f_sferm: z must not be negativ!");

   return 0.5 * z * (2. + log(z) - f_PS(z));
}

} // namespace gm2

} // namespace flexiblesusy
