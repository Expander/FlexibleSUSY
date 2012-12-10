
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale.hpp"
#include "softsusy.h"

#include <cassert>

Mssm_initial_guesser::Mssm_initial_guesser(Mssm<Two_scale>* mssm_, const QedQcd& oneset_, double mxGuess_, double tanb_, int sgnMu_, const DoubleVector& pars_)
   : Initial_guesser<Two_scale>()
   , mssm(mssm_)
   , oneset(oneset_)
   , mxGuess(mxGuess_)
   , tanb(tanb_)
   , sgnMu(sgnMu_)
   , pars(pars_)
{
   assert(mssm && "Mssm_initial_guesser: Error: pointer to Mssm"
          " cannot be zero");
}

Mssm_initial_guesser::~Mssm_initial_guesser()
{
}

void Mssm_initial_guesser::guess()
{
   double mx = 0.0;
   const static MssmSoftsusy empty;
   double m32 = mssm->displayGravitino();
   // double muCondFirst = mssm->displayMuCond();
   // double maCondFirst = mssm->displayMaCond();

   mssm->setSoftsusy(empty);
   /// These are things that are re-written by the new initialisation
   mssm->setData(oneset);
   mssm->setMw(MW);
   mssm->setM32(m32);
   // mssm->setMuCond(muCondFirst);
   // mssm->setMaCond(maCondFirst);

   double mz = mssm->displayMz();

   mx = mxGuess;

   if (oneset.displayMu() != mz) {
      cout << "WARNING: lowOrg in softsusy.cpp called with oneset at scale\n"
	   << oneset.displayMu() << "\ninstead of " << mz << endl;
   }

   MssmSusy t(mssm->guessAtSusyMt(tanb, oneset));

   t.setLoops(2); /// 2 loops should protect against ht Landau pole
   t.runto(mx);

   mssm->setSusy(t);

   /// Initial guess: B=0, mu=1st parameter, need better guesses
   double m0 = pars.display(1);
   double m12 = pars.display(2);
   double a0 = pars.display(3);

   /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
   /// trilinear scalar coupling to be a0
   ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
   mssm->setSugraBcs(m0, m12, a0);

   mssm->setSusyMu(sgnMu * 1.0);
   mssm->setM3Squared(0.);

   mssm->setScale(mx);
   mssm->run_to(mz);

   if (sgnMu == 1 || sgnMu == -1) mssm->rewsbTreeLevel(sgnMu);

   mssm->physical(0);

   mssm->setThresholds(3);
   mssm->setLoops(2);
}
