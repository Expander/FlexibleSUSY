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

#ifndef SMCW_GUT_CONSTRAINT_H
#define SMCW_GUT_CONSTRAINT_H

#include "two_scale_constraint.hpp"
#include "smcw_two_scale.hpp"
#include "gut_scale_calculator.hpp"

class StandardModelCWGUTConstraint : public Constraint<Two_scale> {
public:
   StandardModelCWGUTConstraint(StandardModelCW<Two_scale>*, double);
   virtual ~StandardModelCWGUTConstraint();
   virtual void apply();
   virtual double estimate_scale() const;
   virtual void update_scale();

private:
   double estimated_scale;
   StandardModelCW<Two_scale>* smcw;
   GUT_scale_calculator<StandardModelCW<Two_scale> > gut_scale_calculator;
};

#endif
