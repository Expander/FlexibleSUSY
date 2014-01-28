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

#ifndef MSSM_SOLVER_H
#define MSSM_SOLVER_H

#include "rg_flow.hpp"
#include "softsusy.h"

namespace flexiblesusy {

class SoftSusy_t;

template<>
class RGFlow<SoftSusy_t> {
public:
   typedef void (*THighScaleBoundaryCondition)(MssmSoftsusy&, const DoubleVector&);

   RGFlow();
   ~RGFlow();

   void solve();

   const sPhysical& displayPhys() const;
   const sProblem& displayProblem() const;

   void setHighScaleBoundaryCondition(THighScaleBoundaryCondition);
   void setMxGuess(double);
   void setSoftHighScalePars(const DoubleVector&);
   void setSignMu(int);
   void setTanBeta(double);
   void setLowScaleBoundaryConditions(const QedQcd&);
   void setGaugeUnification(bool);

private:
   MssmSoftsusy                 fMssmSoftSusy;
   THighScaleBoundaryCondition  fHighScaleBoundaryCondition;
   double                       fMxGuess;
   DoubleVector                 fSoftHighScalePars;
   int                          fSignMu;
   double                       fTanBeta;
   QedQcd                       fLowScaleBoundaryContitions;
   bool                         fGaugeUnification;
};

RGFlow<SoftSusy_t>::RGFlow()
   : fMssmSoftSusy()
   , fHighScaleBoundaryCondition(sugraBcs)
   , fMxGuess(1.9e16)
   , fSoftHighScalePars(3)
   , fSignMu(1)
   , fTanBeta(10)
   , fLowScaleBoundaryContitions()
   , fGaugeUnification(true)
{
}

RGFlow<SoftSusy_t>::~RGFlow()
{
}

void RGFlow<SoftSusy_t>::solve()
{
   fMssmSoftSusy.lowOrg(fHighScaleBoundaryCondition, fMxGuess, fSoftHighScalePars,
                        fSignMu, fTanBeta, fLowScaleBoundaryContitions,
                        fGaugeUnification);
}

const sPhysical& RGFlow<SoftSusy_t>::displayPhys() const
{
   return fMssmSoftSusy.displayPhys();
}

const sProblem& RGFlow<SoftSusy_t>::displayProblem() const
{
   return fMssmSoftSusy.displayProblem();
}

void RGFlow<SoftSusy_t>::setHighScaleBoundaryCondition(THighScaleBoundaryCondition bc)
{
   fHighScaleBoundaryCondition = bc;
}

void RGFlow<SoftSusy_t>::setMxGuess(double mxGuess)
{
   fMxGuess = mxGuess;
}

void RGFlow<SoftSusy_t>::setSoftHighScalePars(const DoubleVector& softHighScalePars)
{
   fSoftHighScalePars = softHighScalePars;
}

void RGFlow<SoftSusy_t>::setSignMu(int signMu)
{
   fSignMu = signMu;
}

void RGFlow<SoftSusy_t>::setTanBeta(double tanBeta)
{
   fTanBeta = tanBeta;
}

void RGFlow<SoftSusy_t>::setLowScaleBoundaryConditions(const QedQcd& lowScaleBC)
{
   fLowScaleBoundaryContitions = lowScaleBC;
}

void RGFlow<SoftSusy_t>::setGaugeUnification(bool gaugeUnification)
{
   fGaugeUnification = gaugeUnification;
}

}

#endif
