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

#include "effective_couplings.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

namespace effective_couplings {

std::complex<double> scaling_function(double tau)
{
   std::complex<double> result;

   if (tau > 1) {
      result = -0.25 * Sqr(Log((1.0 + Sqrt(1.0 - 1.0 / tau)) /
         (1.0 - Sqrt(1.0 - 1.0 / tau))) - std::complex<double>(0,1) * Pi);
   } else if (tau == 1.0) {
      result = 0.25 * Sqr(Pi);
   } else {
      result = Sqr(ArcSin(Sqrt(tau)));
   }

   return result;
}

std::complex<double> AS0(double tau)
{
   return -(tau - scaling_function(tau)) / Sqr(tau);
}

std::complex<double> AS12(double tau)
{
   return 2.0 * (tau + (tau - 1.0) * scaling_function(tau)) / Sqr(tau);
}

std::complex<double> AS1(double tau)
{
   return -(2.0 * Sqr(tau) + 3.0 * tau + 3.0 * (2.0 * tau - 1.0)
      * scaling_function(tau)) / Sqr(tau);
}

std::complex<double> AP12(double tau)
{
   return scaling_function(tau) / tau;
}

} // namespace effective_couplings

} // namespace flexiblesusy
