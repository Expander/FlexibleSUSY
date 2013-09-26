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

#include "snmssm_two_scale.hpp"

namespace flexiblesusy {

SNmssm<Two_scale>::SNmssm()
   : Two_scale_model()
   , precision(1.0e-5)
{
}

SNmssm<Two_scale>::SNmssm(const SoftParsNmssm& softPars)
{
   setSoftPars(softPars);
}

SNmssm<Two_scale>::~SNmssm()
{
}

void SNmssm<Two_scale>::calculate_spectrum()
{
   run_to(maximum(displayMsusy(), MZ));
   physical(3);
   runto(displayMz());
}

void SNmssm<Two_scale>::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   runto(scale, eps);
}

}
