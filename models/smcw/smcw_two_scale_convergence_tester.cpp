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

#include "smcw_two_scale_convergence_tester.hpp"

namespace flexiblesusy {

StandardModelCW_convergence_tester::StandardModelCW_convergence_tester(StandardModelCW<Two_scale>* smcw_, double accuracy_goal_)
   : Convergence_tester_DRbar<StandardModelCW<Two_scale> >(smcw_, accuracy_goal_)
{
}

StandardModelCW_convergence_tester::~StandardModelCW_convergence_tester()
{
}

double StandardModelCW_convergence_tester::max_rel_diff() const
{
   const StandardModelCW<Two_scale>* model = get_model();
   const StandardModelCW<Two_scale>* last_iteration_model = get_last_iteration_model();

   const double dg4 = std::fabs(model->displayGaugeCoupling(4)
                                - last_iteration_model->displayGaugeCoupling(4));
   const double dlamda = std::fabs(model->displayLambda()
                                   - last_iteration_model->displayLambda());
   const double max_diff = std::max(dg4, dlamda);

   return max_diff;
}

unsigned int StandardModelCW_convergence_tester::max_iterations() const
{
   return 10;
}

}
