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
namespace gm2os {

class MSSMNoFV_onshell;

double amuWHnu(const MSSMNoFV_onshell&);
double amuWHmuL(const MSSMNoFV_onshell&);
double amuBHmuL(const MSSMNoFV_onshell&);
double amuBHmuR(const MSSMNoFV_onshell&);
double amuBmuLmuR(const MSSMNoFV_onshell&);
double amu1Lapprox(const MSSMNoFV_onshell&);

Eigen::Matrix<std::complex<double>,4,2> n_L(const MSSMNoFV_onshell&);
Eigen::Matrix<std::complex<double>,4,2> n_R(const MSSMNoFV_onshell&);

Eigen::Array<std::complex<double>,2,1> c_L(const MSSMNoFV_onshell&);
Eigen::Array<std::complex<double>,2,1> c_R(const MSSMNoFV_onshell&);

Eigen::Array<double,2,1> AAC(const MSSMNoFV_onshell&);
Eigen::Matrix<double,4,2> AAN(const MSSMNoFV_onshell&);
Eigen::Array<double,2,1> BBC(const MSSMNoFV_onshell&);
Eigen::Matrix<double,4,2> BBN(const MSSMNoFV_onshell&);

Eigen::Matrix<double,4,2> x_im(const MSSMNoFV_onshell&);
Eigen::Array<double,2,1> x_k(const MSSMNoFV_onshell&);

double amuChi0(const MSSMNoFV_onshell&);
double amuChipm(const MSSMNoFV_onshell&);
double calculate_gm2_1loop(const MSSMNoFV_onshell&);

} // namespace gm2os
} // namespace flexiblesusy

#endif
