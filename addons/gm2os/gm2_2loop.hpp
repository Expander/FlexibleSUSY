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

#ifndef GM2_2LOOP_H
#define GM2_2LOOP_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace gm2 {

class Gm2_calculator;

double LogNorm(const Gm2_calculator&);

double tan_beta_cor(const Gm2_calculator&);

double Deltag1(const Gm2_calculator&);
double DeltaYukHiggsino(const Gm2_calculator&);
double DeltaYukBinoHiggsino(const Gm2_calculator&);
double Deltag2(const Gm2_calculator&);
double DeltaYukWinoHiggsino(const Gm2_calculator&);
double DeltaTanBeta(const Gm2_calculator&);

double amuWHnu2L(const Gm2_calculator&);
double amuWHmuL2L(const Gm2_calculator&);
double amuBHmuL2L(const Gm2_calculator&);
double amuBHmuR2L(const Gm2_calculator&);
double amuBmuLmuR2L(const Gm2_calculator&);
double amu2LFSfapprox(const Gm2_calculator&);

double amuChipmPhotonic(const Gm2_calculator&);
double amuChi0Photonic(const Gm2_calculator&);

double tan_alpha(const Gm2_calculator&);
Eigen::Matrix<std::complex<double>,3,3> lambda_mu_cha(const Gm2_calculator&);
Eigen::Matrix<std::complex<double>,2,2> lambda_stop(const Gm2_calculator&);
Eigen::Matrix<std::complex<double>,2,2> lambda_sbot(const Gm2_calculator&);
Eigen::Matrix<std::complex<double>,2,2> lambda_stau(const Gm2_calculator&);
double amua2LSferm(const Gm2_calculator&);
double amua2LCha(const Gm2_calculator&);

} // namespace gm2
} // namespace flexiblesusy

#endif
