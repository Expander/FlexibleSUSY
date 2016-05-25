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

#ifndef SplitMSSM_FullMSSM_THRESHOLDS_H
#define SplitMSSM_FullMSSM_THRESHOLDS_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace splitmssm_thresholds {

struct Parameters {
   Parameters();

   double g1, g2, g3;  ///< gauge couplings (GUT normalized)
   double gt;          ///< MS-bar top Yukawa coupling of the split-SUSY model
   double At;          ///< trilinear coupling for the stops
   double mu;          ///< bilinear Higgsino coupling
   double mA;          ///< mass of the heavy Higgs doublett
   double m1;          ///< bino mass parameter
   double m2;          ///< wino mass parameter
   double tan_beta;    ///< mixing angle of the heavy Higgs doublett
   double scale;       ///< renormalization scale
   Eigen::Matrix<double,3,3> mq2, mu2, md2, ml2, me2;
};

std::ostream& operator<<(std::ostream&, const Parameters&);

double lambda_tree_level(const Parameters&);
double gYu_tree_level(const Parameters&);
double gYd_tree_level(const Parameters&);
double g2u_tree_level(const Parameters&);
double g2d_tree_level(const Parameters&);

double delta_lambda_1loop_reg(const Parameters&);
double delta_lambda_1loop_phi(const Parameters&);
double delta_lambda_1loop_chi_1(const Parameters&);
double delta_lambda_1loop_chi_1(
   double scale, double mu, double lambda, double gYu, double gYd,
   double g2u, double g2d, double m1, double m2);
double delta_lambda_1loop_chi_2(const Parameters&);
double delta_lambda_1loop_chi_2(
   double scale, double mu, double m2, double g1, double g2, double tan_beta);
double delta_lambda_2loop_phi(const Parameters&);
double delta_lambda_2loop_phi_HSS(const Parameters&);
double delta_gYu_1loop(const Parameters&);
double delta_gYd_1loop(const Parameters&);
double delta_g2u_1loop(const Parameters&);
double delta_g2d_1loop(const Parameters&);
double delta_gt_1loop_chi(
   double scale, double mu, double gYu, double gYd,
   double g2u, double g2d, double m1, double m2);
double delta_m2_1loop_chi(const Parameters&);

} // namespace splitmssm_thresholds
} // namespace flexiblesusy

#endif
