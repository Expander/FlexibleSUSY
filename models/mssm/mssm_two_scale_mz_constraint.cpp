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

#include "mssm_two_scale_mz_constraint.hpp"
#include "mssm_two_scale.hpp"

#include <cassert>

Mssm_mz_constraint::Mssm_mz_constraint(double tanBeta_)
   : Constraint<Two_scale>()
   , mssm(NULL)
   , tanBeta(tanBeta_)
   , scale(MZ)
{
}

Mssm_mz_constraint::~Mssm_mz_constraint()
{
}

void Mssm_mz_constraint::apply()
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");

   update_scale();
   mssm->sparticleThresholdCorrections(tanBeta);
}

double Mssm_mz_constraint::get_scale() const
{
   return scale;
}

void Mssm_mz_constraint::set_model(Two_scale_model* model)
{
#ifdef DEBUG
   Mssm<Two_scale>* tmp = dynamic_cast<Mssm<Two_scale>*>(model);
   if (tmp) {
      mssm = tmp;
   } else {
      FATAL("<Mssm_mz_constraint::set_model>: model pointer "
            << model << " is not of type Mssm<Two_scale>*");
   }
#else
   mssm = static_cast<Mssm<Two_scale>*>(model);
#endif
}

void Mssm_mz_constraint::update_scale()
{
   scale = mssm->displayMz();
}
