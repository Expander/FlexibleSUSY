
#include "mssm_solver.h"

MssmSolver::MssmSolver()
   : RGEFlow()
   , fMssmSoftSusy()
   , fHighScaleBoundaryCondition(sugraBcs)
   , fMxGuess(1.0e16)
   , fSoftHighScalePars(3)
   , fSignMu(1)
   , fTanBeta(10)
   , fGaugeUnification(true)
{
   QedQcd oneset;
   const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;

   oneset.setAlpha(ALPHAS, alphasMZ);
   oneset.setPoleMt(mtop);
   oneset.setMass(mBottom, mbmb);
   oneset.toMz();      ///< Runs SM fermion masses to MZ

   fLowScaleBoundaryContitions = oneset;
}

MssmSolver::~MssmSolver()
{
}

void MssmSolver::lowOrg()
{
   fMssmSoftSusy.lowOrg(fHighScaleBoundaryCondition, fMxGuess, fSoftHighScalePars,
                        fSignMu, fTanBeta, fLowScaleBoundaryContitions,
                        fGaugeUnification);
}


const sPhysical& MssmSolver::displayPhys() const
{
   return fMssmSoftSusy.displayPhys();
}

const sProblem& MssmSolver::displayProblem() const
{
   return fMssmSoftSusy.displayProblem();
}

void MssmSolver::setHighScaleBoundaryCondition(THighScaleBoundaryCondition bc)
{
   fHighScaleBoundaryCondition = bc;
}

void MssmSolver::setMxGuess(double mxGuess)
{
   fMxGuess = mxGuess;
}

void MssmSolver::setSoftHighScalePars(const DoubleVector& softHighScalePars)
{
   fSoftHighScalePars = softHighScalePars;
}

void MssmSolver::setSignMu(int signMu)
{
   fSignMu = signMu;
}

void MssmSolver::setTanBeta(double tanBeta)
{
   fTanBeta = tanBeta;
}

void MssmSolver::setLowScaleBoundaryConditions(const QedQcd& lowScaleBC)
{
   fLowScaleBoundaryContitions = lowScaleBC;
}

void MssmSolver::setGaugeUnification(bool gaugeUnification)
{
   fGaugeUnification = gaugeUnification;
}
