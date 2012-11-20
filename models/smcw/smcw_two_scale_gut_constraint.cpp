
#include "smcw_two_scale_gut_constraint.hpp"
#include <cassert>

StandardModelCWGUTConstraint::StandardModelCWGUTConstraint(StandardModelCW<Two_scale>* smcw_, double estimated_scale_)
   : Constraint<Two_scale>()
   , estimated_scale(estimated_scale_)
   , smcw(smcw_)
   , gut_scale_calculator()
{
   assert(smcw && "Error: pointer to StandardModelCW<Two_scale> cannot be zero");
}

StandardModelCWGUTConstraint::~StandardModelCWGUTConstraint()
{
}

void StandardModelCWGUTConstraint::apply()
{
   const double g1 = smcw->displayGaugeCoupling(1);
   const double g2 = smcw->displayGaugeCoupling(2);
   const double g_mean = 0.5 * (g1 + g2);
   smcw->setGaugeCoupling(1, g_mean);
   smcw->setGaugeCoupling(2, g_mean);
}

double StandardModelCWGUTConstraint::get_scale() const
{
   return estimated_scale;
}

void StandardModelCWGUTConstraint::update_scale()
{
   estimated_scale = gut_scale_calculator.calculateGUTScale(*smcw);
}
