
#ifndef MSSM_SOLVER_H
#define MSSM_SOLVER_H

#include "rg_flow.hpp"
#include "softsusy.h"

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

typedef RGFlow<SoftSusy_t> MssmSolver;

RGFlow<SoftSusy_t>::RGFlow()
   : fMssmSoftSusy()
   , fHighScaleBoundaryCondition(sugraBcs)
   , fMxGuess(1.0e16)
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

#endif
