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

#ifndef SM_TWO_SCALE_EXP_CONSTRAINT_H
#define SM_TWO_SCALE_EXP_CONSTRAINT_H

#include "two_scale_constraint.hpp"
#include <cmath>

class Two_scale;

namespace ewConstants {
   namespace {
      const double vev = 246;
      const double root2 = sqrt(2.0);
      const double mtoprun = 165;
      const double mbrun = 2.9;
      const double mtau = 1.77699;
      const double yt = mtoprun * root2 / vev;
      const double yb = mbrun * root2 / vev;
      const double ytau = mtau * root2 / vev;
      const double MZ = 91.1876;
      const double aem = 1.0 / 127.918; // at MZ
      const double sinthWsq = 0.23122;
      const double alpha1 = 5.0 * aem / (3.0 * (1.0 - sinthWsq));
      const double alpha2 = aem / sinthWsq;
      const double alpha3 = 0.1187; // at MZ
      const double g1 = sqrt(4.0 * M_PI * alpha1);
      const double g2 = sqrt(4.0 * M_PI * alpha2);
      const double g3 = sqrt(4.0 * M_PI * alpha3);
   }
}

template <class T>
class StandardModel;

class StandardModelExpConstraint : public Constraint<Two_scale> {
public:
   StandardModelExpConstraint(StandardModel<Two_scale>*);
   virtual ~StandardModelExpConstraint();
   virtual void apply();
   virtual double estimate_scale() const;
   virtual void update_scale();

private:
   StandardModel<Two_scale>* sm;
};

#endif
