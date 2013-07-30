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

#ifndef SM_TWO_SCALE_CONVERGENCE_TESTER_H
#define SM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "two_scale_convergence_tester.hpp"

namespace flexiblesusy {

template <class T>
class StandardModel;

/**
 * @class StandardModel_convergence_tester
 * @brief Tests the StandardModel<Two_scale> class for convergence
 *
 * This class tests the StandardModel<Two_scale> class for convergence
 * during an iteration between a high and a low scale.  The Standard
 * Model is assumed to be a fixed theory (all parameters are known),
 * the function accuracy_goal_reached() returns allways true.
 */
class StandardModel_convergence_tester : public Convergence_tester<Two_scale> {
public:
   StandardModel_convergence_tester(StandardModel<Two_scale>*);
   virtual ~StandardModel_convergence_tester();
   virtual bool accuracy_goal_reached();

private:
   StandardModel<Two_scale>* sm; ///< the model to test for convergence
};

}

#endif
