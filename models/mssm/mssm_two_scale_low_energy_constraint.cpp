
#include "mssm_two_scale_low_energy_constraint.hpp"
#include <cassert>

Mssm_low_energy_constraint::Mssm_low_energy_constraint(Mssm<Two_scale>* mssm_,
                                                       const QedQcd& oneset_,
                                                       double tanBeta_, double scale_)
   : Constraint<Two_scale>()
   , mssm(mssm_)
   , oneset(oneset_)
   , tanBeta(tanBeta_)
   , scale(scale_)
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");
}

Mssm_low_energy_constraint::~Mssm_low_energy_constraint()
{
}

void Mssm_low_energy_constraint::apply()
{
   mssm->calcDrBarPars();
}

double Mssm_low_energy_constraint::get_scale() const
{
   return scale;
}

void Mssm_low_energy_constraint::update_scale()
{
   drBarPars tree(mssm->displayDrBarPars());
   double tmp_scale = sqrt(tree.mu(2, 3) * tree.mu(1, 3));
   if (tmp_scale > 0.0)
      scale = tmp_scale;
}
