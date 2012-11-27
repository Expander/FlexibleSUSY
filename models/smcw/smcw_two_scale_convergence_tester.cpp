
#include "smcw_two_scale_convergence_tester.hpp"
#include "logger.hpp"

#include <limits>

StandardModelCW_convergence_tester::StandardModelCW_convergence_tester(StandardModelCW<Two_scale>* smcw_, double accuracy_goal_)
   : Convergence_tester<Two_scale>()
   , smcw(smcw_)
   , last_iteration()
   , it_count(0)
   , accuracy_goal(accuracy_goal_)
{
}

StandardModelCW_convergence_tester::~StandardModelCW_convergence_tester()
{
}

bool StandardModelCW_convergence_tester::accuracy_goal_reached()
{
   bool precision_reached;
   if (it_count == 0) {
      precision_reached = false;
   } else {
      if (scale_has_changed())
         WARNING("scale has changed by " << scale_difference()
                 << " GeV (" << rel_scale_difference()
                 << "%), parameter comparison might fail");

      const double dg4 = std::fabs(smcw->displayGaugeCoupling(4) - last_iteration.displayGaugeCoupling(4));
      const double dlamda = std::fabs(smcw->displayLambda() - last_iteration.displayLambda());
      const double max_diff = std::max(dg4, dlamda);
      precision_reached = max_diff < accuracy_goal;
   }

   // save old model parameters
   last_iteration = *smcw;
   ++it_count;

   return precision_reached;
}

bool StandardModelCW_convergence_tester::scale_has_changed() const
{
   return !is_zero(scale_difference());
}

double StandardModelCW_convergence_tester::scale_difference() const
{
   return smcw->getScale() - last_iteration.getScale();
}

double StandardModelCW_convergence_tester::rel_scale_difference() const
{
   const double diff = scale_difference();
   if (!is_zero(last_iteration.getScale()))
      return diff / last_iteration.getScale();
   return std::numeric_limits<double>::infinity();
}
