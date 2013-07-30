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

#include "mssm_two_scale_low_scale_constraint.hpp"
#include "mssm_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

Mssm_low_scale_constraint::Mssm_low_scale_constraint(const Mssm_parameter_point& pp_)
   : Constraint<Two_scale>()
   , mssm(NULL)
   , scale(MZ)
   , pp(pp_)
{
}

Mssm_low_scale_constraint::~Mssm_low_scale_constraint()
{
}

void Mssm_low_scale_constraint::apply()
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");

   update_scale();
   mssm->sparticleThresholdCorrections(pp.tanBeta);
}

double Mssm_low_scale_constraint::get_scale() const
{
   return scale;
}

void Mssm_low_scale_constraint::set_model(Two_scale_model* model)
{
   mssm = cast_model<Mssm<Two_scale> >(model);
}

void Mssm_low_scale_constraint::update_scale()
{
   scale = mssm->displayMz();
}

}
