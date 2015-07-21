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

#include "gm2_calculator.hpp"
#include "wrappers.hpp"
#include "logger.hpp"

namespace flexiblesusy {
namespace gm2 {

Gm2_calculator::Gm2_calculator(const MSSMNoFV_mass_eigenstates& model_)
   : MSSMNoFV_mass_eigenstates(model_),
   MW(80.385), MZ(91.1876), TB(10.), EL(0.303), ME(0.00051),
   MM(0.105), ML(1.777), MU(0.04), MC(1.5), MT(173.5), MD(0.04), MS(0.15),
   MB(3.), gY(0.35), EL0(0.30282212), MA0(500.)
{}

void Gm2_calculator::convert_parameters_reverse() {
   set_TB(get_vu() / get_vd());
   double gY = get_gY();
   double g2 = get_g2();
   set_EL(gY * g2 / sqrt(sqr(gY) + sqr(g2)));
   double v = sqrt(sqr(get_vu()) + sqr(get_vd()));
   set_MZ(0.5 * sqrt(sqr(gY) + sqr(g2)) * v);
   set_MW(0.5 * g2 * v);
}


void Gm2_calculator::convert_parameters() {
   set_g1(sqrt(5. / 3.) * EL * MZ / MW);
   set_g2(EL / sqrt(1. - sqr(MW / MZ)));
   double v = 2. * MW / get_g2();
   set_vu(v / sqrt(1. + 1. / sqr(TB)));
   set_vd(get_vu() / TB);
   Eigen::Matrix<double,3,3> Ye_neu(get_Ye());
   Ye_neu(0, 0) = sqrt(2.) * ME / get_vd();
   Ye_neu(1, 1) = sqrt(2.) * MM / get_vd();
   Ye_neu(2, 2) = sqrt(2.) * ML / get_vd();
   set_Ye(Ye_neu);
   Eigen::Matrix<double,3,3> Yu_neu(get_Yu());
   Yu_neu(0, 0) = sqrt(2.) * MU / get_vu();
   Yu_neu(1, 1) = sqrt(2.) * MC / get_vu();
   Yu_neu(2, 2) = sqrt(2.) * MT / get_vu();
   set_Yu(Yu_neu);
   Eigen::Matrix<double,3,3> Yd_neu(get_Yd());
   Yd_neu(0, 0) = sqrt(2.) * MD / get_vd();
   Yd_neu(1, 1) = sqrt(2.) * MS / get_vd();
   Yd_neu(2, 2) = sqrt(2.) * MB / get_vd();
   set_Yd(Yd_neu);
   set_TYe(Ye_neu * Ae);
   set_TYu(Yu_neu * Au);
   set_TYd(Yd_neu * Ad);
   double tan2b = 2. * TB / (1. - sqr(TB));
   set_BMu(0.5 * sqr(MA0) * (tan2b / sqrt(1. + sqr(tan2b))));
}

void Gm2_calculator::calculate_DRbar_parameters() {
   convert_parameters();
   MSSMNoFV_mass_eigenstates::calculate_DRbar_parameters();

   MSmu = get_MSm();
   MStau = get_MStau();
   MSbot = get_MSb();
   MStop = get_MSt();
}

} // gm2
} // namespace flexiblesusy
