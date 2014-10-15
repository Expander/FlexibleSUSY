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
   : theta_12(Electroweak_constants::CKM_THETA12)
   , theta_13(Electroweak_constants::CKM_THETA13)
   , theta_23(Electroweak_constants::CKM_THETA23)
   , delta(Electroweak_constants::CKM_DELTA)
{
}

/**
 * Calculates V_CKM angles from Wolfenstein parameters (see
 * hep-ph/0406184)
 */
void CKM_parameters::set_from_wolfenstein(double lambdaW, double aCkm,
                                          double rhobar, double etabar)
{
   theta_12 = ArcSin(lambdaW);
   theta_23 = ArcSin(aCkm * Sqr(lambdaW));

   const double lambdaW3 = lambdaW * Sqr(lambdaW);
   const double lambdaW4 = lambdaW * lambdaW3;

   const std::complex<double> rpe(rhobar, etabar);
   const std::complex<double> V13conj = aCkm * lambdaW3 * rpe
      * Sqrt(1.0 - Sqr(aCkm) * lambdaW4) /
      Sqrt(1.0 - Sqr(lambdaW)) / (1.0 - Sqr(aCkm) * lambdaW4 * rpe);

   theta_13 = ArcSin(Abs(V13conj));
   delta = arg(V13conj);
}

/**
 * Calculates Wolfenstein parameters from V_CKM angles (see
 * hep-ph/0406184)
 */
void CKM_parameters::get_wolfenstein(double& lambdaW, double& aCkm,
                                     double& rhobar, double& etabar)
{
   const double sin_12 = Sin(theta_12);
   const double sin_13 = Sin(theta_13);
   const double sin_23 = Sin(theta_23);

   // Eq. (11.4) from PDG
   lambdaW  = sin_12;
   aCkm     = sin_23 / Sqr(lambdaW);

   const double c = Sqrt((1.0 - Sqr(sin_23)) / (1.0 - Sqr(lambdaW)));
   const std::complex<double> eid(std::polar(1.0, delta));
   const std::complex<double> r(sin_13 * eid /
      (c * aCkm * lambdaW * Sqr(lambdaW) + sin_13 * eid * Sqr(sin_23)));

   rhobar = Re(r);
   etabar = Im(r);
}

} // namespace flexiblesusy
