
#include "mssm_two_scale_sugra_constraint.hpp"
#include <cassert>

Mssm_sugra_constraint::Mssm_sugra_constraint(double mx_guess_, double m0_, double m12_, double a0_, int signMu_)
   : Constraint<Two_scale>()
   , mx_guess(mx_guess_)
   , mssm(NULL)
   , gut_scale_calculator()
   , m0(m0_)
   , m12(m12_)
   , a0(a0_)
   , signMu(signMu_)
{
}

Mssm_sugra_constraint::~Mssm_sugra_constraint()
{
}

void Mssm_sugra_constraint::apply()
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");

   update_scale();
   mssm->setSugraBcs(m0, m12, a0);
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
