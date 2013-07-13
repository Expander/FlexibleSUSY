
#include "mssm_two_scale_sugra_constraint.hpp"
#include <cassert>

Mssm_sugra_constraint::Mssm_sugra_constraint(const Mssm_parameter_point& pp_)
   : Constraint<Two_scale>()
   , mx_guess(pp_.mxGuess)
   , mssm(NULL)
   , pp(pp_)
   , gut_scale_calculator()
{
}

Mssm_sugra_constraint::~Mssm_sugra_constraint()
{
}

void Mssm_sugra_constraint::apply()
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");

   update_scale();
   mssm->setSugraBcs(pp.m0, pp.m12, pp.a0);
}

double Mssm_sugra_constraint::get_scale() const
{
   return mx_guess;
}

void Mssm_sugra_constraint::set_model(Two_scale_model* model)
{
   mssm = cast_model<Mssm<Two_scale> >(model);
}

void Mssm_sugra_constraint::update_scale()
{
   mx_guess = gut_scale_calculator.calculateGUTScale(*mssm);
}
