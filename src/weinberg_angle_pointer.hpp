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

#ifndef WEINBERG_ANGLE_POINTER_H
#define WEINBERG_ANGLE_POINTER_H

#include "CMSSM_two_scale_model.hpp"

namespace flexiblesusy {

namespace weinberg_angle {

/**
 * @class Weinberg_angle_pointer
 * @brief Class to calculate the DR-bar weak mixing angle from pointer to a model
 */
class Weinberg_angle_pointer {
public:
    Weinberg_angle_pointer(const CMSSM<Two_scale>*);
   ~Weinberg_angle_pointer();

   void set_number_of_iterations(unsigned);         ///< maximum number of iterations
   void set_number_of_loops(unsigned);              ///< set number of loops
   void set_precision_goal(double);                 ///< set precision goal
   void set_model_pointer(const CMSSM<Two_scale>*); ///< set pointer to investigated model

   /// calculates and returns the sinus of the Weinberg angle
   double calculate(double rho_start = 1.0, double sin_start = 0.48);

private:
   /**
    * @class Derived_data
    * @brief Derived model parameters necessary for calculating weak mixing angle
    */
   struct Derived_data {
      Derived_data();

      double alpha_em_drbar;      ///< alpha_em(MZ, DR-bar, SUSY)
      double msel_drbar;          ///< left-handed selectron DR-bar mass
      double msmul_drbar;         ///< left-handed smuon DR-bar mass
      double msve_drbar;          ///< electron-sneutrino DR-bar mass
      double msvm_drbar;          ///< muon-sneutrino DR-bar mass
      double self_energy_z_at_mz; ///< self-energy Z at p = MZ, mt = mt_pole
      double self_energy_w_at_0;  ///< self-energy W at p = 0, mt = mt_pole
      double self_energy_w_at_mw; ///< self-energy W at p = MW, mt = mt_pole
   };

   unsigned number_of_iterations; ///< maximum number of iterations
   unsigned number_of_loops;      ///< number of loops
   double precision_goal;         ///< precision goal
   const CMSSM<Two_scale>* model; ///< pointer to investigated model
   Derived_data derived_data;

   void calculate_derived_data();
   double calculate_self_energy_z_top(double, double);
   double calculate_self_energy_w_top(double, double);

   double calculate_delta_rho(double, double);
   double calculate_delta_r(double, double);
   double calculate_delta_vb(double, double);
   double calculate_delta_vb_sm(double, double);
   double calculate_delta_vb_susy(double);
   static double rho_2(double);
};

} // namespace weinberg_angle

} // namespace flexiblesusy

#endif
