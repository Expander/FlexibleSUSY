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

#ifndef ONE_LOOP_DECAY_DIAGRAMS_H
#define ONE_LOOP_DECAY_DIAGRAMS_H

#include "decay_amplitudes.hpp"

#include <complex>

namespace flexiblesusy {

Decay_amplitude_SSS calculate_diagram_SSS_t1g1n1_FFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g2n2_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g3n3_UUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g4n4_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g5n5_SVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g6n6_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g7n7_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g8n8_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g9n9_VVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t1g10n10_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t2g1n11_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t2g2n12_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t3g1n13_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t3g2n14_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t4g1n15_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t4g2n16_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t5g1n17_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t5g2n18_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t6g1n19_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t6g2n20_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t7g1n21_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSS calculate_diagram_SSS_t7g2n22_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g1n23_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g2n24_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g3n25_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g4n26_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g5n27_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g6n28_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g7n29_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g8n30_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g9n31_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t8g10n32_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g1n33_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g2n34_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g3n35_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g4n36_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g5n37_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g6n38_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g7n39_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g8n40_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g9n41_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t9g10n42_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g1n43_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g2n44_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g3n45_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g4n46_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g5n47_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g6n48_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g7n49_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g8n50_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g9n51_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSS calculate_diagram_SSS_t10g10n52_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g1n1_FFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g2n2_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g3n3_UUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g4n4_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g5n5_SVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g6n6_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g7n7_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g8n8_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g9n9_VVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t1g10n10_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t2g1n11_VS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t2g2n12_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t3g1n13_VS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t3g2n14_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t4g1n15_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t4g2n16_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t5g1n17_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t5g2n18_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t6g1n19_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t7g1n20_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SVV calculate_diagram_SVV_t8g6n26_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t8g7n27_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t8g8n28_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t8g9n29_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t8g10n30_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t9g6n36_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t9g7n37_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t9g8n38_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t9g9n39_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t9g10n40_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g1n41_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g2n42_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g3n43_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g4n44_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g5n45_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g6n46_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g7n47_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g8n48_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g9n49_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SVV calculate_diagram_SVV_t10g10n50_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g1n1_FFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g2n2_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g3n3_UUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g4n4_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g5n5_SVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g6n6_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g7n7_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g8n8_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g9n9_VVS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t1g10n10_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t2g1n11_VS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t2g2n12_VV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSV calculate_diagram_SSV_t3g1n13_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t3g2n14_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSV calculate_diagram_SSV_t4g1n15_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t5g1n16_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t5g2n17_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double, double);

Decay_amplitude_SSV calculate_diagram_SSV_t6g1n18_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   double);

Decay_amplitude_SSV calculate_diagram_SSV_t8g6n26_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t8g7n27_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t8g9n29_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t8g10n30_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g1n31_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g2n32_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g3n33_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g4n34_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g5n35_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g6n36_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g7n37_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g8n38_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g9n39_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t9g10n40_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g1n41_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g2n42_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g3n43_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g4n44_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g5n45_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g6n46_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g7n47_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g8n48_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g9n49_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SSV calculate_diagram_SSV_t10g10n50_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g1n1_FFS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g2n2_SSF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g3n3_FFV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g4n4_SVF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g5n5_VSF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t1g6n6_VVF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t2g1n7_SS(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t2g2n8_SV(double, double, double,
   double, double, const std::complex<double>&, const std::complex<double>&,
   const std::complex<double>&, double, double);

Decay_amplitude_SFF calculate_diagram_SFF_t3g1n9_FFS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t3g2n10_FFV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SFF calculate_diagram_SFF_t4g1n11_FFS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t4g2n12_FFV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g1n13_SFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g2n14_SSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g3n15_SUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g4n16_SSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g5n17_SVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double,
   double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g6n18_VFF(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, const
   std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g7n19_VSS(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g8n20_VUU(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g9n21_VSV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

Decay_amplitude_SFF calculate_diagram_SFF_t5g10n22_VVV(double, double, double,
   double, double, double, const std::complex<double>&, const std::complex<
   double>&, const std::complex<double>&, const std::complex<double>&, double);

} // namespace flexiblesusy

#endif
