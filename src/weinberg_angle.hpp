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
      double mw_pole;                ///< W pole mass
      double mz_pole;                ///< Z pole mass
      double mt_pole;                ///< top quark pole mass
      double mh_drbar;               ///< lightest CP-even Higgs DR-bar mass
      double hmix_12;                ///< CP-even Higgs mixing Cos(alpha)
      double msel_drbar;             ///< left-handed selectron DR-bar mass
      double msmul_drbar;            ///< left-handed smuon DR-bar mass
      double msve_drbar;             ///< electron-sneutrino DR-bar mass
      double msvm_drbar;             ///< muon-sneutrino DR-bar mass
      Eigen::ArrayXd mn_drbar;       ///< Neutralino DR-bar mass
      Eigen::ArrayXd mc_drbar;       ///< Chargino DR-bar mass
      Eigen::MatrixXcd zn;           ///< Neutralino mixing matrix
      Eigen::MatrixXcd um;           ///< Chargino mixing matrix
      Eigen::MatrixXcd up;           ///< Chargino mixing matrix
      double gY;                     ///< U(1)_Y gauge coupling
      double g2;                     ///< SU(2)_L gauge coupling
      double g3;                     ///< SU(3)_c gauge coupling
      double tan_beta;               ///< tan(beta) = vu / vd
      double ymu;                    ///< Myon Yukawa coupling
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
   mutable double rho_hat;        ///< output rho-hat parameter
   Data data;

   void rhohat(double&, double&, const Data&) const;
   static double calculate_delta_r(double, double, const Data&);
   static double calculate_delta_rho(double, double, const Data&);
   static double calculate_delta_vb(double, double, const Data&);
   static double rho_2(double);

   double calculate_sin(double rho_start = 1.0, double sin_start = 0.48) const;
};

} // namespace weinberg_angle

} // namespace flexiblesusy

#endif
