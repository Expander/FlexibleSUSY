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

#ifndef NMSSM_TWOLOOPHIGGS_H
#define NMSSM_TWOLOOPHIGGS_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace nmssm_twoloophiggs {

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_at_as_nmssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq,
   double amu, double tanb, double vev2, double gs, double svevS);

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_ab_as_nmssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq,
   double amu, double cotbeta, double vev2, double gs, double svevS);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as, double amu);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as, double amu);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as, double amu);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as, double amu);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as);

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as);

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as);

} // namespace nmssm_twoloophiggs
} // namespace flexiblesusy

#endif
