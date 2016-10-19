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

#ifndef SPECTRUM_GENERATOR_SETTINGS_H
#define SPECTRUM_GENERATOR_SETTINGS_H

#include "two_loop_corrections.hpp"
#include <array>
#include <iosfwd>

namespace flexiblesusy {

/**
 * @class Spectrum_generator_settings
 * @brief stores the spectrum generator settings
 *
 * This class stores all spectrum generator settings which can be
 * changed via the SLHA input file.
 */
class Spectrum_generator_settings {
public:
   /// Spectrum generator settings
   enum Settings : unsigned {
      precision,             ///< [0] overall precision goal
      max_iterations,        ///< [1] maximum number of iterations (0 = automatic)
      algorithm,             ///< [2] RG solver algorithm (0 = two-scale)
      calculate_sm_masses,   ///< [3] calculate Standard Model pole masses
      pole_mass_loop_order,  ///< [4] loop-order for calculation of pole masses
      ewsb_loop_order,       ///< [5] loop-order for solving the EWSB eqs.
      beta_loop_order,       ///< [6] loop-order of beta-functions
      threshold_corrections_loop_order, ///< [7]  threshold corrections loop order
      higgs_2loop_correction_at_as,     ///< [8]  Higgs 2-loop correction O(alpha_t alpha_s)
      higgs_2loop_correction_ab_as,     ///< [9]  Higgs 2-loop correction O(alpha_b alpha_s)
      higgs_2loop_correction_at_at,     ///< [10] Higgs 2-loop correction O(alpha_t alpha_t + alpha_t alpha_b + alpha_b alpha_b)
      higgs_2loop_correction_atau_atau, ///< [11] Higgs 2-loop correction O(alpha_tau alpha_tau)
      force_output,          ///< [12] force output
      top_pole_qcd_corrections,         ///< [13] Top-quark pole mass QCD corrections
      beta_zero_threshold,   ///< [14] beta function zero threshold
      calculate_observables, ///< [15] calculate observables (a_muon, ...)
      force_positive_masses, ///< [16] force positive masses of majoran fermions
      pole_mass_scale,       ///< [17] renormalization scale at which the pole masses are calculated
      eft_pole_mass_scale,   ///< [18] renormalization scale at which the pole masses are calculated in the EFT
      eft_matching_scale,    ///< [19] renormalization scale at which the EFT is matched to the full model
      eft_matching_loop_order_up, ///< [20] loop order at which the gauge and Yukawa couplings of the full model are calculated from the EFT ones (upwards matching)
      eft_matching_loop_order_down, ///< [21] loop order at which lambda of the SM is calculated from the full model parameters at the matching scale (downwards matching)
      eft_higgs_index,       ///< [22] index of SM-Higgs in Higgs multiplet
      calculate_bsm_masses,  ///< [23] calculate BSM pole masses
      NUMBER_OF_OPTIONS      ///< number of possible options
   };

   Spectrum_generator_settings();

   double get(Settings) const; ///< get value of spectrum generator setting
   void set(Settings, double); ///< set value of spectrum generator setting
   void reset();               ///< resets all settings to their defaults

   Two_loop_corrections get_two_loop_corrections() const;
   void set_two_loop_corrections(const Two_loop_corrections&);

private:
   std::array<double, NUMBER_OF_OPTIONS> values; ///< spectrum generator settings
};

std::ostream& operator<<(std::ostream&, const Spectrum_generator_settings&);

} // namespace flexiblesusy

#endif
