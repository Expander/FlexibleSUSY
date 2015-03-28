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

#ifndef MSSMCBS_TWO_SCALE_LOW_SCALE_CONSTRAINT_H
#define MSSMCBS_TWO_SCALE_LOW_SCALE_CONSTRAINT_H

#include "MSSMcbs_low_scale_constraint.hpp"
#include "two_scale_constraint.hpp"
#include "lowe.h"
#include <Eigen/Core>

namespace flexiblesusy {

template <class T>
class MSSMcbs;

class Two_scale;

template<>
class MSSMcbs_low_scale_constraint<Two_scale> : public Constraint<Two_scale> {
public:
   MSSMcbs_low_scale_constraint();
   MSSMcbs_low_scale_constraint(MSSMcbs<Two_scale>*, const QedQcd&);
   virtual ~MSSMcbs_low_scale_constraint();
   virtual void apply();
   virtual double get_scale() const;
   virtual void set_model(Two_scale_model*);

   void clear();
   double get_initial_scale_guess() const;
   void initialize();
   void set_sm_parameters(const QedQcd&);
   void set_threshold_corrections_loop_order(unsigned); ///< threshold corrections loop order

private:
   double scale;
   double initial_scale_guess;
   MSSMcbs<Two_scale>* model;
   QedQcd oneset;
   double MZDRbar;
   double new_g1, new_g2, new_g3;
   unsigned threshold_corrections_loop_order; ///< threshold corrections loop order

   void calculate_DRbar_gauge_couplings();
   void calculate_DRbar_yukawa_couplings();
   void calculate_Yu_DRbar();
   void calculate_Yd_DRbar();
   void calculate_Ye_DRbar();
   double calculate_delta_alpha_em(double) const;
   double calculate_alS5DRbar_over_alS5MSbar(double) const;
   double calculate_zeta_g_QCD_2(double) const;
   double calculate_zeta_g_SUSY_2(double) const;
   void update_scale();
};

} // namespace flexiblesusy

#endif
