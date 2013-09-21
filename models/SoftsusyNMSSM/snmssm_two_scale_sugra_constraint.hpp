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

#ifndef SoftsusyNMSSM_SUGRA_CONSTRAINT_H
#define SoftsusyNMSSM_SUGRA_CONSTRAINT_H

#include "two_scale_constraint.hpp"
#include "snmssm_two_scale.hpp"
#include "snmssm_parameter_point.hpp"
#include "gut_scale_calculator.hpp"

namespace flexiblesusy {

class SNmssm_sugra_constraint : public Constraint<Two_scale> {
public:
   SNmssm_sugra_constraint(const SNmssm_parameter_point&);
   virtual ~SNmssm_sugra_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

private:
   double mx_guess;
   SNmssm<Two_scale>* mssm;
   SNmssm_parameter_point pp;   ///< SNmssm parameter point
   GUT_scale_calculator<SNmssm<Two_scale> > gut_scale_calculator;

   void update_scale();
};

} // namespace flexiblesusy

#endif
