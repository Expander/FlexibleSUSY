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

#ifndef MSSM_PARAMETER_POINT_H
#define MSSM_PARAMETER_POINT_H

#include "linalg.h"
#include "lowe.h"

struct Mssm_parameter_point {
   Mssm_parameter_point()
      : m0(125.0)
      , m12(500.0)
      , a0(0.0)
      , mxGuess(1.0e16)
      , signMu(1)
      , tanBeta(10.0)
      , oneset()
   {
      const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
      oneset.setAlpha(ALPHAS, alphasMZ);
      oneset.setPoleMt(mtop);
      oneset.setMass(mBottom, mbmb);
      oneset.toMz();
   }
   DoubleVector get_soft_pars() const {
      DoubleVector highScaleSoftPars(3);
      highScaleSoftPars(1) = m0;
      highScaleSoftPars(2) = m12;
      highScaleSoftPars(3) = a0;
      return highScaleSoftPars;
   }
   friend std::ostream& operator<<(std::ostream& os, const Mssm_parameter_point& pp) {
      os << "CMSSM parameter point:"
         << " m0=" << pp.m0
         << ", m12=" << pp.m12
         << ", a0=" << pp.a0
         << ", mxGuess=" << pp.mxGuess
         << ", signMu=" << pp.signMu
         << ", tanBeta=" << pp.tanBeta
         << ", QedQcd=" << pp.oneset
         << '\n';
      return os;
   }

   double m0, m12, a0, mxGuess, signMu, tanBeta;
   QedQcd oneset;
};

#endif
