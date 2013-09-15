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

#include "snmssm_two_scale_susy_scale_constraint.hpp"
#include "snmssm_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

/**
 * Constructor
 *
 * @param pp_ SNmssm parameter point
 */
SNmssm_susy_scale_constraint::SNmssm_susy_scale_constraint(const SNmssm_parameter_point& pp_)
   : Constraint<Two_scale>()
   , snmssm(NULL)
   , scale(pp_.msGuess)
   , pp(pp_)
{
}

SNmssm_susy_scale_constraint::~SNmssm_susy_scale_constraint()
{
}

void SNmssm_susy_scale_constraint::apply()
{
   assert(snmssm && "Error: pointer to SNmssm<Two_scale> cannot be zero");

   snmssm->calcDrBarPars();
   update_scale();
   const double mtrun = snmssm->displayDrBarPars().mt;
   snmssm->rewsb(pp.signMu, mtrun);
}

double SNmssm_susy_scale_constraint::get_scale() const
{
   return scale;
}

void SNmssm_susy_scale_constraint::set_model(Two_scale_model* model)
{
   snmssm = cast_model<SNmssm<Two_scale> >(model);
}

void SNmssm_susy_scale_constraint::update_scale()
{
   snmssm->setMsusy(snmssm->calcMs());
   scale = snmssm->displayMsusy();
}

}
