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
class Mssm<Two_scale>: public Two_scale_model, public MssmSoftsusy {
public:
   Mssm();
   virtual ~Mssm();

   virtual void calculate_spectrum();
   virtual std::string name() const { return "Mssm"; }
   virtual int run_to(double, double eps = -1.0);
   virtual void print(std::ostream& s) const { s << static_cast<MssmSoftsusy>(*this); }

   void setScale(double scale) { setMu(scale); }
   double getScale() const { return displayMu(); }
   Mssm calcBeta() const { return beta2(); }
   void setSugraBcs(double m0, double m12, double a0) { standardSugra(m0, m12, a0); }

private:
   Mssm(const SoftParsMssm&);
};

#endif
