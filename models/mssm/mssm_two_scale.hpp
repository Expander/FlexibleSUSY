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

#ifndef MSSM_TWO_SCALE_H
#define MSSM_TWO_SCALE_H

#include "mssm.hpp"
#include "two_scale_model.hpp"
#include "softsusy.h"

class Two_scale;

template<>
class Mssm<Two_scale>: public Two_scale_model {
public:
   Mssm();
   virtual ~Mssm();

   void init(const QedQcd&, double, double, int, DoubleVector);
   virtual std::string name() const { return "Mssm"; }
   virtual int run_to(double);
   virtual void print(std::ostream& s) const { s << mssm; }

   double displayGaugeCoupling(int i) const { return mssm.displayGaugeCoupling(i); }
   const drBarPars& displayDrBarPars() const { return mssm.displayDrBarPars(); }
   double getScale() const { return mssm.displayMu(); }
   const sPhysical& displayPhys() const { return mssm.displayPhys(); }
   double displayGravitino() const { return mssm.displayGravitino(); }
   double displayMz() const { return mssm.displayMz(); }

   void setScale(double scale) { mssm.setMu(scale); }
   void setGaugeCoupling(int i, double g) { mssm.setGaugeCoupling(i, g); }
   void setData(const QedQcd& d) { mssm.setData(d); }
   void setTanb(double tanb) { mssm.setTanb(tanb); }
   void setSusyMu(double mu) { mssm.setSusyMu(mu); }
   void setSugraBcs(double m0, double m12, double a0) { mssm.standardSugra(m0, m12, a0); }
   void setMw(double mw) { mssm.setMw(mw); }
   void setM32(double m32) { mssm.setM32(m32); }
   void setSusy(const MssmSusy& s) { mssm.setSusy(s); }

   Mssm calcBeta() const { return mssm.beta2(); }
   void calcDrBarPars() { mssm.calcDrBarPars(); }
   MssmSusy guessAtSusyMt(double tanb, const QedQcd& oneset) { return mssm.guessAtSusyMt(tanb, oneset); }

private:
   Mssm(const SoftParsMssm&);

   MssmSoftsusy mssm;
};

#endif
