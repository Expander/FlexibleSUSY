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

#include "observables.hpp"

namespace flexiblesusy {

Observables::Observables()
   : a_muon(0.)
   , a_muon_gm2calc(0.)
   , a_muon_gm2calc_uncertainty(0.)
{
}

Eigen::ArrayXd Observables::get() const
{
   Eigen::ArrayXd vec(NUMBER_OF_OBSERVABLES);
   vec(0) = a_muon;
   vec(1) = a_muon_gm2calc;
   vec(2) = a_muon_gm2calc_uncertainty;
   return vec;
}

std::vector<std::string> Observables::get_names()
{
   std::vector<std::string> names(Observables::NUMBER_OF_OBSERVABLES);
   names[0] = "a_muon";
   names[1] = "a_muon_gm2calc";
   names[2] = "a_muon_gm2calc_uncertainty";
   return names;
}

void Observables::clear()
{
   a_muon = 0.;
   a_muon_gm2calc = 0.;
   a_muon_gm2calc_uncertainty = 0.;
}

void Observables::set(const Eigen::ArrayXd& vec)
{
   assert(vec.rows() == NUMBER_OF_OBSERVABLES);
   a_muon = vec(0);
   a_muon_gm2calc = vec(1);
   a_muon_gm2calc_uncertainty = vec(2);
}

} // namespace flexiblesusy
