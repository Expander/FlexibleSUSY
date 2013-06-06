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

#ifndef MSSM_MSUSY_CONSTRAINT_H
#define MSSM_MSUSY_CONSTRAINT_H

#include "two_scale_constraint.hpp"
#include "linalg.h"

class Two_scale;
template<class T> class Mssm;

/**
 * @class Mssm_msusy_constraint
 * @brief MSSM EWSB constraint at the Susy mass scale
 *
 * This class represents the electroweak symmetry breaking (EWSB)
 * constraint of the MSSM at the Susy mass scale MSusy.  The apply()
 * function calculates the MSSM \f$\overline{DR}\f$ parameters and
 * does the EWSB.
 */

class Mssm_msusy_constraint : public Constraint<Two_scale> {
public:
   Mssm_msusy_constraint(const DoubleVector& pars_, double scale_, int sgnMu_);
   virtual ~Mssm_msusy_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

private:
   Mssm<Two_scale>* mssm;
   DoubleVector pars;
   double scale;
   int sgnMu;

   void update_scale();
};

#endif
