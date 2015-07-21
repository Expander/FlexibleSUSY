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

#include "logger.hpp"
#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "MSSMNoFV_onshell.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

int main()
{
   MSSMNoFV_mass_eigenstates model;
   gm2os::MSSMNoFV_onshell osmodel(model);

   const double gm2_1l = gm2os::calculate_gm2_1loop(osmodel);
   const double gm2_2l = amu2LFSfapprox(osmodel)
                        + amuChipmPhotonic(osmodel) + amuChi0Photonic(osmodel);
   const double gm2_2l_tanb_approx =  + (tan_beta_cor(osmodel) - 1.) * gm2_1l;

   return 0;
}
