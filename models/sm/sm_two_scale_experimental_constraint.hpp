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

#ifndef SM_TWO_SCALE_EXP_CONSTRAINT_H
#define SM_TWO_SCALE_EXP_CONSTRAINT_H

#include "two_scale_constraint.hpp"
#include "ew_input.hpp"
#include <cmath>

namespace flexiblesusy {

class Two_scale;

template <class T>
class StandardModel;

/**
 * @class StandardModel_exp_constraint
 * @brief Experimental constraints of the Standard Model
 *
 * This class applies all experimental constraints to the Standard
 * Model at the Z mass scale.
 */
class StandardModel_exp_constraint : public Constraint<Two_scale> {
public:
   StandardModel_exp_constraint();
   virtual ~StandardModel_exp_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

private:
   StandardModel<Two_scale>* sm; ///< model to apply the constraints to
};

}

#endif
