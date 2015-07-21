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

#ifndef GM2_1LOOP_H
#define GM2_1LOOP_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace gm2 {

class Gm2_calculator;

double amuWHnu(const Gm2_calculator&);
double amuWHmuL(const Gm2_calculator&);
double amuBHmuL(const Gm2_calculator&);
double amuBHmuR(const Gm2_calculator&);
double amuBmuLmuR(const Gm2_calculator&);
double amu1Lapprox(const Gm2_calculator&);

Eigen::Matrix<std::complex<double>,4,2> n_L(const Gm2_calculator&);
Eigen::Matrix<std::complex<double>,4,2> n_R(const Gm2_calculator&);

Eigen::Array<std::complex<double>,2,1> c_L(const Gm2_calculator&);
Eigen::Array<std::complex<double>,2,1> c_R(const Gm2_calculator&);

Eigen::Array<double,2,1> AAC(const Gm2_calculator&);
Eigen::Matrix<double,4,2> AAN(const Gm2_calculator&);
Eigen::Array<double,2,1> BBC(const Gm2_calculator&);
Eigen::Matrix<double,4,2> BBN(const Gm2_calculator&);

Eigen::Matrix<double,4,2> x_im(const Gm2_calculator&);
Eigen::Array<double,2,1> x_k(const Gm2_calculator&);

double amuChi0(const Gm2_calculator&);
double amuChipm(const Gm2_calculator&);
double calculate_gm2_1loop(const Gm2_calculator&);

} // namespace gm2
} // namespace flexiblesusy

#endif
