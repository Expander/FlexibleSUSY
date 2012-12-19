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

#include "mssm_two_scale_msusy_constraint.hpp"
#include "mssm_two_scale.hpp"

#include <cassert>

/**
 * Constructor
 *
 * @param mssm_ Mssm class to apply the constraint to
 * @param pars_ vector with soft parameters m0, m12, a0
 * @param scale_ first guess for the susy scale
 * @param sgnMu_ sign of superpotential parameter \f$\mu\f$
 */
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
   update_scale();
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
}
