
#ifndef MSSM_SOLVER_H
#define MSSM_SOLVER_H

#include "rge_flow.hpp"
#include "softsusy.h"

class SoftSusy_t;

template<>
class RGEFlow<SoftSusy_t> {
public:
   typedef void (*THighScaleBoundaryCondition)(MssmSoftsusy&, const DoubleVector&);

   RGEFlow();
   ~RGEFlow();

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

typedef RGEFlow<SoftSusy_t> MssmSolver;

RGEFlow<SoftSusy_t>::RGEFlow()
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

RGEFlow<SoftSusy_t>::~RGEFlow()
{
}

void RGEFlow<SoftSusy_t>::solve()
{
   fMssmSoftSusy.lowOrg(fHighScaleBoundaryCondition, fMxGuess, fSoftHighScalePars,
                        fSignMu, fTanBeta, fLowScaleBoundaryContitions,
                        fGaugeUnification);
}

const sPhysical& RGEFlow<SoftSusy_t>::displayPhys() const
{
   return fMssmSoftSusy.displayPhys();
}

const sProblem& RGEFlow<SoftSusy_t>::displayProblem() const
{
   return fMssmSoftSusy.displayProblem();
}

void RGEFlow<SoftSusy_t>::setHighScaleBoundaryCondition(THighScaleBoundaryCondition bc)
{
   fHighScaleBoundaryCondition = bc;
}

void RGEFlow<SoftSusy_t>::setMxGuess(double mxGuess)
{
   fMxGuess = mxGuess;
}

void RGEFlow<SoftSusy_t>::setSoftHighScalePars(const DoubleVector& softHighScalePars)
{
   fSoftHighScalePars = softHighScalePars;
}

void RGEFlow<SoftSusy_t>::setSignMu(int signMu)
{
   fSignMu = signMu;
}

void RGEFlow<SoftSusy_t>::setTanBeta(double tanBeta)
{
   fTanBeta = tanBeta;
}

void RGEFlow<SoftSusy_t>::setLowScaleBoundaryConditions(const QedQcd& lowScaleBC)
{
   fLowScaleBoundaryContitions = lowScaleBC;
}

void RGEFlow<SoftSusy_t>::setGaugeUnification(bool gaugeUnification)
{
   fGaugeUnification = gaugeUnification;
}

#endif
