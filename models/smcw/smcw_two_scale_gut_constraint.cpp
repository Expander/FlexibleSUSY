
#include "smcw_two_scale_gut_constraint.hpp"
#include <cassert>

namespace flexiblesusy {

StandardModelCWGUTConstraint::StandardModelCWGUTConstraint(double estimated_scale_, double lambda_at_mgut_)
   : Constraint<Two_scale>()
   , estimated_scale(estimated_scale_)
   , smcw(NULL)
   , gut_scale_calculator()
   , lambda_at_mgut(lambda_at_mgut_)
{
}

StandardModelCWGUTConstraint::~StandardModelCWGUTConstraint()
{
}

void StandardModelCWGUTConstraint::apply()
{
   assert(smcw && "Error: pointer to StandardModelCW<Two_scale> cannot be zero");

   update_scale();

   const double g1 = smcw->displayGaugeCoupling(1);
   const double g2 = smcw->displayGaugeCoupling(2);
   const double g_mean = 0.5 * (g1 + g2);
   smcw->setGaugeCoupling(1, g_mean);
   smcw->setGaugeCoupling(2, g_mean);
   smcw->setGaugeCoupling(4, g_mean);
   smcw->setLambda(lambda_at_mgut);
}

double StandardModelCWGUTConstraint::get_scale() const
{
   return estimated_scale;
}

void StandardModelCWGUTConstraint::set_model(Two_scale_model* model)
{
   smcw = cast_model<StandardModelCW<Two_scale> >(model);
}

void StandardModelCWGUTConstraint::update_scale()
{
   estimated_scale = gut_scale_calculator.calculateGUTScale(*smcw);
}

}
