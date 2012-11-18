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

StandardModelExpConstraint::StandardModelExpConstraint(StandardModel<Two_scale>* sm_)
   : Constraint<Two_scale>()
   , sm(sm_)
{
   assert(sm && "pointer to StandardModel<Two_scale> must not be zero");
}

StandardModelExpConstraint::~StandardModelExpConstraint()
{
}

void StandardModelExpConstraint::apply()
{
   VERBOSE_MSG("Applying SM experimental constraints at scale "
               << sm->displayMu())
   if (std::fabs(ewConstants::MZ - sm->displayMu()) > 1.0)
      WARNING("Applying the experimental constraints "
              "of StandardModel<Two_scale> at scale " << sm->displayMu()
              << " != MZ is not save!")

   sm->setYukawaElement(StandardModel<Two_scale>::YU, 3, 3, ewConstants::yt);
   sm->setYukawaElement(StandardModel<Two_scale>::YD, 3, 3, ewConstants::yb);
   sm->setYukawaElement(StandardModel<Two_scale>::YE, 3, 3, ewConstants::ytau);
   sm->setGaugeCoupling(1, sqrt(4 * M_PI * ewConstants::alpha1));
   sm->setGaugeCoupling(2, sqrt(4 * M_PI * ewConstants::alpha2));
   sm->setGaugeCoupling(3, sqrt(4 * M_PI * ewConstants::alpha3));
}

double StandardModelExpConstraint::get_scale() const
{
   return ewConstants::MZ;
}
