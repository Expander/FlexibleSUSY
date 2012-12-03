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

#include "mssm_two_scale.hpp"

Mssm<Two_scale>::Mssm()
   : Two_scale_model()
   , mssm()
{
}

Mssm<Two_scale>::Mssm(const SoftParsMssm& softPars)
{
   mssm.setSoftPars(softPars);
}

Mssm<Two_scale>::~Mssm()
{
}

void Mssm<Two_scale>::init(const QedQcd& oneset, double mxGuess, double tanb, int sgnMu, DoubleVector pars)
{
   double mx = 0.0;
   const static MssmSoftsusy empty;
   double m32 = mssm.displayGravitino();
   double muCondFirst = mssm.displayMuCond();
   double maCondFirst = mssm.displayMaCond();

   mssm.setSoftsusy(empty);
   /// These are things that are re-written by the new initialisation
   mssm.setData(oneset);
   mssm.setMw(MW);
   mssm.setM32(m32);
   mssm.setMuCond(muCondFirst);
   mssm.setMaCond(maCondFirst);

   double mz = mssm.displayMz();

   mx = mxGuess;

   if (oneset.displayMu() != mz) {
      cout << "WARNING: lowOrg in softsusy.cpp called with oneset at scale\n"
	   << oneset.displayMu() << "\ninstead of " << mz << endl;
   }

   MssmSusy t(mssm.guessAtSusyMt(tanb, oneset));

   t.setLoops(2); /// 2 loops should protect against ht Landau pole
   t.runto(mx);

   mssm.setSusy(t);

   /// Initial guess: B=0, mu=1st parameter, need better guesses
   double m0 = pars.display(1);
   double m12 = pars.display(2);
   double a0 = pars.display(3);

   /// Sets scalar soft masses equal to m0, fermion ones to m12 and sets the
   /// trilinear scalar coupling to be a0
   ///  if (m0 < 0.0) m.flagTachyon(true); Deleted on request from A Pukhov
   mssm.standardSugra(m0, m12, a0);

   mssm.setSusyMu(sgnMu * 1.0);
   mssm.setM3Squared(0.);

   mssm.run(mx, mz);

   if (sgnMu == 1 || sgnMu == -1) mssm.rewsbTreeLevel(sgnMu);

   mssm.physical(0);

   mssm.setThresholds(3);
   mssm.setLoops(2);
}

int Mssm<Two_scale>::run_to(double scale)
{
   return mssm.runto(scale);
}
