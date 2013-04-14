
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

   const double g1 = mssm->displayGaugeCoupling(1);
   const double g2 = mssm->displayGaugeCoupling(2);
   const double g_mean = 0.5 * (g1 + g2);
   mssm->setGaugeCoupling(1, g_mean);
   mssm->setGaugeCoupling(2, g_mean);
   mssm->setSugraBcs(m0, m12, a0);
}

double Mssm_sugra_constraint::get_scale() const
{
   return mx_guess;
}

void Mssm_sugra_constraint::set_model(Two_scale_model* model)
{
#ifdef DEBUG
   Mssm<Two_scale>* tmp = dynamic_cast<Mssm<Two_scale>*>(model);
   if (tmp) {
      mssm = tmp;
   } else {
      FATAL("<Mssm_sugra_constraint::set_model>: model pointer "
            << model << " is not of type Mssm<Two_scale>*");
   }
#else
   mssm = static_cast<Mssm<Two_scale>*>(model);
#endif
}

void Mssm_sugra_constraint::update_scale()
{
   mx_guess = gut_scale_calculator.calculateGUTScale(*mssm);
}
