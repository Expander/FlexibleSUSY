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

#include "mssm_two_scale_susy_scale_constraint.hpp"
#include "mssm_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

/**
 * Constructor
 *
 * @param pp_ Mssm parameter point
 */
Mssm_susy_scale_constraint::Mssm_susy_scale_constraint(const Mssm_parameter_point& pp_)
   : Constraint<Two_scale>()
   , mssm(NULL)
   , scale(pp_.msGuess)
   , pp(pp_)
{
}

Mssm_susy_scale_constraint::~Mssm_susy_scale_constraint()
{
}

void Mssm_susy_scale_constraint::apply()
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");

   mssm->calcDrBarPars();
   update_scale();
   const double mtrun = mssm->displayDrBarPars().mt;
   mssm->rewsb(pp.signMu, mtrun, pp.get_soft_pars());
}

double Mssm_susy_scale_constraint::get_scale() const
{
   return scale;
}

void Mssm_susy_scale_constraint::set_model(Two_scale_model* model)
{
   mssm = cast_model<Mssm<Two_scale> >(model);
}

void Mssm_susy_scale_constraint::update_scale()
{
   mssm->setMsusy(mssm->calcMs());
   scale = mssm->displayMsusy();
}

}
