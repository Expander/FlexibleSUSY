
#include "sm_two_scale_convergence_tester.hpp"
#include <cassert>

StandardModel_convergence_tester::StandardModel_convergence_tester(StandardModel<Two_scale>* sm_)
   : Convergence_tester<Two_scale>()
   , sm(sm_)
{
   assert(sm && "pointer to StandardModel<Two_scale> must not be zero");
}

StandardModel_convergence_tester::~StandardModel_convergence_tester()
{
}

/**
 * Returns allways true, because the Standard Model is assumed to be a
 * fixed theory.
 *
 * @return true
 */
bool StandardModel_convergence_tester::accuracy_goal_reached()
{
   return true;
}
