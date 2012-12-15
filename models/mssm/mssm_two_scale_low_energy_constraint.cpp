
#include "mssm_two_scale_low_energy_constraint.hpp"
#include <cassert>

Mssm_mz_constraint::Mssm_mz_constraint(Mssm<Two_scale>* mssm_,
                                       double tanBeta_)
   : Constraint<Two_scale>()
   , mssm(mssm_)
   , tanBeta(tanBeta_)
   , scale(MZ)
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");
}

Mssm_mz_constraint::~Mssm_mz_constraint()
{
}

void Mssm_mz_constraint::apply()
{
   mssm->sparticleThresholdCorrections(tanBeta);
}

double Mssm_mz_constraint::get_scale() const
{
   return scale;
}

void Mssm_mz_constraint::update_scale()
{
   scale = mssm->displayMz();
}

Mssm_msusy_constraint::Mssm_msusy_constraint(Mssm<Two_scale>* mssm_,
                                             const DoubleVector& pars_,
                                             double scale_, int sgnMu_)
   : Constraint<Two_scale>()
   , mssm(mssm_)
   , pars(pars_)
   , scale(scale_)
   , sgnMu(sgnMu_)
{
   assert(mssm && "Error: pointer to Mssm<Two_scale> cannot be zero");
}

Mssm_msusy_constraint::~Mssm_msusy_constraint()
{
}

void Mssm_msusy_constraint::apply()
{
   mssm->calcDrBarPars();
   double mtrun = mssm->displayDrBarPars().mt;
   mssm->rewsb(sgnMu, mtrun, pars);
}

double Mssm_msusy_constraint::get_scale() const
{
   return scale;
}

void Mssm_msusy_constraint::update_scale()
{
   mssm->setMsusy(mssm->calcMs());
   scale = mssm->displayMsusy();
   // drBarPars tree(mssm->displayDrBarPars());
   // double tmp_scale = sqrt(tree.mu(2, 3) * tree.mu(1, 3));
   // if (tmp_scale > 0.0)
   //    scale = tmp_scale;
}
