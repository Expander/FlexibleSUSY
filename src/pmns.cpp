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

#include "pmns.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

void PMNS_parameters::reset_to_diagonal()
{
   theta_12 = 0.;
   theta_13 = 0.;
   theta_23 = 0.;
   delta    = 0.;
   alpha_1  = 0.;
   alpha_2  = 0.;
}

void PMNS_parameters::reset_to_observation()
{
   theta_12 = Electroweak_constants::PMNS_THETA12;
   theta_13 = Electroweak_constants::PMNS_THETA13;
   theta_23 = Electroweak_constants::PMNS_THETA23;
   delta    = Electroweak_constants::PMNS_DELTA;
   alpha_1  = Electroweak_constants::PMNS_ALPHA1;
   alpha_2  = Electroweak_constants::PMNS_ALPHA2;
}

Eigen::Matrix<double,3,3> PMNS_parameters::get_real_pmns() const
{
   const std::complex<double> eID(std::polar(1.0, delta));
   const double s12 = Sin(theta_12);
   const double s13 = Sin(theta_13);
   const double s23 = Sin(theta_23);
   const double c12 = Cos(theta_12);
   const double c13 = Cos(theta_13);
   const double c23 = Cos(theta_23);

   // set phase factor e^(i delta) to +1 or -1 depending on the sign
   // of its real part
   const int pf = Sign(Re(eID));

   Eigen::Matrix<double,3,3> pmns_matrix;
   pmns_matrix(0, 0) = c12 * c13;
   pmns_matrix(0, 1) = s12 * c13;
   pmns_matrix(0, 2) = pf * s13;
   pmns_matrix(1, 0) = -s12 * c23 - pf * c12 * s23 * s13;
   pmns_matrix(1, 1) = c12 * c23 - pf * s12 * s23 * s13;
   pmns_matrix(1, 2) = s23 * c13;
   pmns_matrix(2, 0) = s12 * s23 - pf * c12 * c23 * s13;
   pmns_matrix(2, 1) = -c12 * s23 - pf * s12 * c23 * s13;
   pmns_matrix(2, 2) = c23 * c13;

   return pmns_matrix;
}

Eigen::Matrix<std::complex<double>,3,3> PMNS_parameters::get_complex_pmns() const
{
   const std::complex<double> eID(std::polar(1.0, delta));
   const std::complex<double> eIAlpha1(std::polar(1.0, 0.5*alpha_1));
   const std::complex<double> eIAlpha2(std::polar(1.0, 0.5*alpha_2));
   const double s12 = Sin(theta_12);
   const double s13 = Sin(theta_13);
   const double s23 = Sin(theta_23);
   const double c12 = Cos(theta_12);
   const double c13 = Cos(theta_13);
   const double c23 = Cos(theta_23);

   Eigen::Matrix<std::complex<double>,3,3> pmns_matrix;
   pmns_matrix(0, 0) = c12 * c13 * eIAlpha1;
   pmns_matrix(0, 1) = s12 * c13 * eIAlpha2;
   pmns_matrix(0, 2) = s13 / eID;
   pmns_matrix(1, 0) = (-s12 * c23 - c12 * s23 * s13 * eID) * eIAlpha1;
   pmns_matrix(1, 1) = (c12 * c23 - s12 * s23 * s13 * eID) * eIAlpha2;
   pmns_matrix(1, 2) = s23 * c13;
   pmns_matrix(2, 0) = (s12 * s23 - c12 * c23 * s13 * eID) * eIAlpha1;
   pmns_matrix(2, 1) = (-c12 * s23 - s12 * c23 * s13 * eID) * eIAlpha2;
   pmns_matrix(2, 2) = c23 * c13;

   return pmns_matrix;
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<double,3,3>& Vv,
                                        Eigen::Matrix<double,3,3>& Ve,
                                        Eigen::Matrix<double,3,3>& Ue)
{
   Eigen::Matrix<double,3,3> pmns(Ve*Vv.adjoint());
   to_pdg_convention(pmns, Vv, Ve, Ue);
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<double,3,3>& pmns,
                                        Eigen::Matrix<double,3,3>& Vv,
                                        Eigen::Matrix<double,3,3>& Ve,
                                        Eigen::Matrix<double,3,3>& Ue)
{
   Eigen::Matrix<double,3,3> signs_E(Eigen::Matrix<double,3,3>::Identity());
   Eigen::Matrix<double,3,3> signs_V(Eigen::Matrix<double,3,3>::Identity());

   // make 33 element positive
   if (pmns(2, 2) < 0.) {
      signs_E(2, 2) = -1.;
      for (int j = 0; j < 3; ++j) {
         pmns(2, j) *= -1.;
      }
   }

   // make 23 element positive
   if (pmns(1, 2) < 0.) {
      signs_V(2, 2) = -1;
      signs_E(2, 2) *= -1;
      for (int j = 0; j < 3; ++j) {
         pmns(2, j) *= -1;
         pmns(j, 2) *= -1;
      }
   }

   Ve = signs_E * Ve;
   Ue = signs_E * Ue;
   Vv = signs_V * Vv;
}

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<std::complex<double>,3,3>& Vv,
                                        Eigen::Matrix<std::complex<double>,3,3>& Ve,
                                        Eigen::Matrix<std::complex<double>,3,3>& Ue)
{
   Eigen::Matrix<std::complex<double>,3,3> pmns(Ve*Vv.adjoint());
   to_pdg_convention(pmns, Vv, Ve, Ue);
}

namespace {

/// restrict sin or cos to interval [-1,1]
double sanitize_hypot(double sc)
{
   if (sc < -1.) sc = -1.;
   if (sc > 1.) sc = 1.;
   return sc;
}

} // anonymous namespace

void PMNS_parameters::to_pdg_convention(Eigen::Matrix<std::complex<double>,3,3>& pmns,
                                        Eigen::Matrix<std::complex<double>,3,3>& Vv,
                                        Eigen::Matrix<std::complex<double>,3,3>& Ve,
                                        Eigen::Matrix<std::complex<double>,3,3>& Ue)
{
   const double s13 = sanitize_hypot(std::abs(pmns(0,2)));
   const double c13 = std::sqrt(1 - Sqr(s13));

}

} // namespace flexiblesusy
