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

#ifndef MSSM_TWOLOOPHIGGS_H
#define MSSM_TWOLOOPHIGGS_H

#include <Eigen/Core>

namespace flexiblesusy {
namespace mssm_twoloophiggs {

// tadpoles

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq,
   double amu, double tanb, double vev2, double gs);

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   double rmtsq, double rmbsq, double mAsq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2);

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq,
   double amu, double cotbeta, double vev2, double gs);

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   double rmtausq, double mAsq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2);

// self-energies

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double amu,
   double tanb, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   double rmtsq, double rmbsq, double fmasq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double amu,
   double cotbeta, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   double rmtausq, double fmasq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2, int scheme = 0);


Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_as_mssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double amu,
   double tanb, double vev2, double gs);

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_at_mssm(
   double rmtsq, double rmbsq, double fmasq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2);

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double amu,
   double cotbeta, double vev2, double gs);

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_atau_atau_mssm(
   double rmtausq, double fmasq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2);

// self-energies with tadpoles added

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double amu,
   double tanb, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
   double rmtsq, double rmbsq, double fmasq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double amu,
   double cotbeta, double vev2, double gs, int scheme = 0);

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
   double rmtausq, double fmasq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2, int scheme = 0);


double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double amu,
   double tanb, double vev2, double gs);

double self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   double rmtsq, double rmbsq, double fmasq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2);

double self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double amu,
   double cotbeta, double vev2, double gs);

double self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
   double rmtausq, double fmasq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2);

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy

#endif
