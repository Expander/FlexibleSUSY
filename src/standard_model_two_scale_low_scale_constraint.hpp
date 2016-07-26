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


#ifndef STANDARD_MODEL_TWO_SCALE_LOW_SCALE_CONSTRAINT_H
#define STANDARD_MODEL_TWO_SCALE_LOW_SCALE_CONSTRAINT_H

#include "standard_model_low_scale_constraint.hpp"
#include "two_scale_constraint.hpp"
#include "lowe.h"
#include <Eigen/Core>

namespace flexiblesusy {

class Two_scale;

namespace standard_model {

template <class T>
class StandardModel;

template<>
class Standard_model_low_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   Standard_model_low_scale_constraint();
   Standard_model_low_scale_constraint(StandardModel<Two_scale>*, const softsusy::QedQcd&);
   virtual ~Standard_model_low_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   void initialize();
   const softsusy::QedQcd& get_sm_parameters() const;
   void set_sm_parameters(const softsusy::QedQcd&);
   void set_threshold_corrections_loop_order(unsigned); ///< threshold corrections loop order

private:
   double scale;
   StandardModel<Two_scale>* model;
   softsusy::QedQcd qedqcd;
};

} // namespace standard_model
} // namespace flexiblesusy

#endif
