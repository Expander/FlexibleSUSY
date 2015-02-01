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

#ifndef WEINBERG_ANGLE_H
#define WEINBERG_ANGLE_H

#include <Eigen/Core>

namespace flexiblesusy {

namespace weinberg_angle {

class Weinberg_angle {
public:
   struct Data {
      Data();

      double scale;                  ///< renormalization scale
      double alpha_em_drbar;         ///< alpha_em(MZ, DR-bar, SUSY)
      double fermi_contant;          ///< Fermi constant
      double self_energy_z_at_mz;    ///< self-energy Z at p = MZ
      double self_energy_w_at_0;     ///< self-energy W at p = 0
      double self_energy_w_at_mw;    ///< self-energy W at p = MW
      double mw_pole;
      double mz_pole;
      double mt_pole;
      double mh_drbar;
      double hmix_12;
      double mse_L;
      double msmu_L;
      double msnu_e;
      double msnu_mu;
      Eigen::ArrayXd mneut;
      Eigen::ArrayXd mch;
      Eigen::MatrixXcd zn;
      Eigen::MatrixXcd um;
      Eigen::MatrixXcd up;
      double gY;
      double g2;
      double g3;
      double tan_beta;
      double ymu;
   };

   Weinberg_angle();
   ~Weinberg_angle();

   void set_data(const Data&);
   void set_number_of_iterations(unsigned);
   void set_precision_goal(double);
   double get_rho_hat() const;

   double calculate() const;

private:
   unsigned number_of_iterations; ///< maximum number of iterations
   double precision_goal;         ///< precision goal
   double rho_hat;                ///< output rho-hat parameter
   Data data;

   double rho_2(double) const;

   double calculate_delta_r(double, double, const Data&) const;
   double calculate_delta_rho(double, double, const Data&) const;
   double calculate_delta_vb(double, double, const Data&) const;
};

} // namespace weinberg_angle

} // namespace flexiblesusy

#endif
