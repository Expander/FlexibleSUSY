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

#ifndef SMCW_TWO_SCALE_CONVERGENCE_TESTER_H
#define SMCW_TWO_SCALE_CONVERGENCE_TESTER_H

#include "two_scale_convergence_tester.hpp"
#include "smcw_two_scale.hpp"

class StandardModelCW_convergence_tester : public Convergence_tester<Two_scale> {
public:
   StandardModelCW_convergence_tester(StandardModelCW<Two_scale>*, double);
   virtual ~StandardModelCW_convergence_tester();
   virtual bool accuracy_goal_reached();
   virtual unsigned int max_iterations() const;

private:
   StandardModelCW<Two_scale>* smcw;
   StandardModelCW<Two_scale> last_iteration;
   unsigned long it_count;
   double accuracy_goal;

   double scale_difference() const;
   double rel_scale_difference() const;
   bool scale_has_changed() const;
};

#endif
