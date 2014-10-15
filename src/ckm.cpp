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

#include "ckm.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

CKM_parameters::CKM_parameters()
{
   set_angles(Electroweak_constants::CKM_THETA12,
              Electroweak_constants::CKM_THETA13,
              Electroweak_constants::CKM_THETA23,
              Electroweak_constants::CKM_DELTA);
}

/**
 * Initializes Wolfenstein parameters from V_CKM angles (see
 * hep-ph/0406184)
 */
void CKM_parameters::set_angles(double theta_12, double theta_13,
                                double theta_23, double delta)
{
   const double sin_12 = Sin(theta_12);
   const double sin_13 = Sin(theta_13);
   const double sin_23 = Sin(theta_23);
   const std::complex<double> V13(sin_13 * std::polar(1.0, - delta));

   // Eq. (7)
   lambda  = sin_12;
   A       = sin_23 / Sqr(lambda);

   const double rho =   Re(V13 / (A * Power(lambda, 3)));
   const double eta = - Im(V13 / (A * Power(lambda, 3)));

   // Eq. (15)
   rho_bar = rho - 0.5 * rho * Sqr(lambda)
      + (0.5 * Sqr(A) * rho - 0.125 * rho - Sqr(A) * (Sqr(rho) - Sqr(eta)))
      * Power(lambda,4);

   // Eq. (16)
   eta_bar = eta - 0.5 * eta * Sqr(lambda)
      + (0.5 * Sqr(A) * eta - 0.125 * eta - 2.0 * Sqr(A) * rho * eta)
      * Power(lambda,4);
}

} // namespace flexiblesusy
