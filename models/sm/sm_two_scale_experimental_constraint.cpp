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

#include "sm_two_scale_experimental_constraint.hpp"
#include "sm_two_scale.hpp"
#include "logger.hpp"

#include <cassert>

namespace flexiblesusy {

StandardModel_exp_constraint::StandardModel_exp_constraint()
   : Constraint<Two_scale>()
   , sm(NULL)
{
}

StandardModel_exp_constraint::~StandardModel_exp_constraint()
{
}

/**
 * Apply all experimental constraints to the Standard Model class.  It
 * is assumed that the Standard Model class is at the Z mass scale.
 * If this is not the case, a warning is printed.
 */
void StandardModel_exp_constraint::apply()
{
   assert(sm && "pointer to StandardModel<Two_scale> must not be zero");

   VERBOSE_MSG("Applying SM experimental constraints at scale "
               << sm->get_scale());
   if (std::fabs(Electroweak_constants::MZ - sm->get_scale()) > 1.0)
      WARNING("Applying the experimental constraints "
              "of StandardModel<Two_scale> at scale " << sm->get_scale()
              << " != MZ is not save!");

   sm->setYukawaElement(StandardModel<Two_scale>::YU, 3, 3, Electroweak_constants::yt);
   sm->setYukawaElement(StandardModel<Two_scale>::YD, 3, 3, Electroweak_constants::yb);
   sm->setYukawaElement(StandardModel<Two_scale>::YE, 3, 3, Electroweak_constants::ytau);
   sm->setGaugeCoupling(1, Electroweak_constants::g1);
   sm->setGaugeCoupling(2, Electroweak_constants::g2);
   sm->setGaugeCoupling(3, Electroweak_constants::g3);
}

/**
 * Returns the on-shell Z mass MZ
 * @return MZ
 */
double StandardModel_exp_constraint::get_scale() const
{
   return Electroweak_constants::MZ;
}

void StandardModel_exp_constraint::set_model(Two_scale_model* model)
{
   sm = cast_model<StandardModel<Two_scale>*>(model);
}

}
