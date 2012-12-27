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

#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale.hpp"
#include "softsusy.h"
#include "logger.hpp"

#include <cassert>

Mssm_initial_guesser::Mssm_initial_guesser(Mssm<Two_scale>* mssm_, const QedQcd& oneset_, double mxGuess_, double tanb_, int sgnMu_, const DoubleVector& pars_, bool ewsbBCscale_)
   : Initial_guesser<Two_scale>()
   , mssm(mssm_)
   , oneset(oneset_)
   , mxGuess(mxGuess_)
   , tanb(tanb_)
   , sgnMu(sgnMu_)
   , pars(pars_)
   , ewsbBCscale(ewsbBCscale_)
{
   assert(mssm && "Mssm_initial_guesser: Error: pointer to Mssm"
          " cannot be zero");
}

Mssm_initial_guesser::~Mssm_initial_guesser()
{
}

void Mssm_initial_guesser::guess()
{
   const static MssmSoftsusy empty;

   double mx = 0.0;
   const double muFirst = mssm->displaySusyMu(); /// Remember initial values
   const bool setTbAtMXflag = mssm->displaySetTbAtMX();
   const bool altFlag = mssm->displayAltEwsb();
   const double m32 = mssm->displayGravitino();
   const double muCondFirst = mssm->displayMuCond();
   const double maCondFirst = mssm->displayMaCond();

   mssm->setSoftsusy(empty); /// Always starts from an empty object
   mssm->setSetTbAtMX(setTbAtMXflag);
   if (altFlag)
      mssm->useAlternativeEwsb();
   mssm->setData(oneset);
   mssm->setMw(MW);
   mssm->setM32(m32);
   mssm->setMuCond(muCondFirst);
   mssm->setMaCond(maCondFirst);

   double mz = mssm->displayMz();

   if (mxGuess > 0.0) {
      mx = mxGuess;
   } else {
      string ii("Trying to use negative mx in MssmSoftsusy::lowOrg.\n");
      ii = ii + "Now illegal! Use positive mx for first guess of mx.\n";
      throw ii;
   }

   if (oneset.displayMu() != mz) {
      WARNING("lowOrg in softsusy.cpp called with oneset at scale\n"
              << oneset.displayMu() << "instead of " << mz);
   }

   MssmSusy t(mssm->guessAtSusyMt(tanb, oneset));
   t.setLoops(2); /// 2 loops should protect against ht Landau pole
   t.runto(mx);

   mssm->setSusy(t);

   /// Initial guess: B=0, mu=1st parameter, need better guesses
   const double m0 = pars.display(1);
   const double m12 = pars.display(2);
   const double a0 = pars.display(3);
   mssm->standardSugra(m0, m12, a0);

   if ((sgnMu == 1 || sgnMu == -1) && !ewsbBCscale) {
      mssm->setSusyMu(sgnMu * 1.0);
      mssm->setM3Squared(0.);
   } else {
      if (mssm->displayAltEwsb()) {
         mssm->setSusyMu(mssm->displayMuCond());
         mssm->setM3Squared(mssm->displayMaCond());
      } else {
         mssm->setSusyMu(muFirst);
         mssm->setM3Squared(muFirst);
      }
   }

   mssm->run(mx, mz);

   if (sgnMu == 1 || sgnMu == -1)
      mssm->rewsbTreeLevel(sgnMu);

   mssm->physical(0);
   mssm->setThresholds(3);
   mssm->setLoops(2);

   // test for problems
   sProblem problem = mssm->displayProblem();
   if (problem.test()) {
      if (problem.testSeriousProblem())
         WARNING("<Mssm_initial_guesser::guess()>: serious problem: " << problem);
      else
         WARNING("<Mssm_initial_guesser::guess()>: problem: " << problem);
   }
}
