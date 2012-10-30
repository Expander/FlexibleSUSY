
#ifndef MSSM_SOLVER_H
#define MSSM_SOLVER_H

#include "rge_flow.hpp"
#include "softsusy.h"

class MssmSolver : public RGEFlow {
public:
   typedef void (*THighScaleBoundaryCondition)(MssmSoftsusy&, const DoubleVector&);

   MssmSolver();
   virtual ~MssmSolver();

   virtual void lowOrg();

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

#endif
