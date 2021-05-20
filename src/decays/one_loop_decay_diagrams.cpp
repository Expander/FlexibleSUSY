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

#include "one_loop_decay_diagrams.hpp"
#include "wrappers.hpp"
#include "loop_libraries/loop_library.hpp"

namespace flexiblesusy {

namespace {

inline double Den(double x2, double y2) {
   return 1.0/(x2-y2);
}

} // anonymous

/**
 * @brief Evaluates T1G1N1 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F6S3PL coupling Cp[F[5], F[6], S[3]][PL]
 * @param[in] CpF5F6S3PR coupling Cp[F[5], F[6], S[3]][PR]
 * @param[in] CpcF6F4S2PL coupling Cp[-F[6], F[4], S[2]][PL]
 * @param[in] CpcF6F4S2PR coupling Cp[-F[6], F[4], S[2]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g1n1_FFF(
   double mext1, double mext2, double mext3,
   double mLF4, double mLF5, double mLF6,
   const std::complex<double>& CpcF4cF5S1PL, const std::complex<double>&
      CpcF4cF5S1PR, const std::complex<double>& CpF5F6S3PL, const std::complex<
      double>& CpF5F6S3PR, const std::complex<double>& CpcF6F4S2PL, const std::
      complex<double>& CpcF6F4S2PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(2*b0tmp1*(
      CpcF4cF5S1PR*(CpcF6F4S2PR*CpF5F6S3PL*mLF4 + CpcF6F4S2PL*CpF5F6S3PR*mLF5 +
      CpcF6F4S2PL*CpF5F6S3PL*mLF6) + CpcF4cF5S1PL*(CpcF6F4S2PL*CpF5F6S3PR*mLF4
      + CpcF6F4S2PR*CpF5F6S3PL*mLF5 + CpcF6F4S2PR*CpF5F6S3PR*mLF6)) + c1tmp3*
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mLF4*Sqr(mext1) + 3*c2tmp4*
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mLF4*Sqr(mext1) + c1tmp3*CpcF4cF5S1PL
      *CpcF6F4S2PL*CpF5F6S3PR*mLF4*Sqr(mext1) + 3*c2tmp4*CpcF4cF5S1PL*
      CpcF6F4S2PL*CpF5F6S3PR*mLF4*Sqr(mext1) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*
      CpF5F6S3PL*mLF5*Sqr(mext1) + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*
      mLF5*Sqr(mext1) + c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mLF6*Sqr(
      mext1) + 2*c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mLF6*Sqr(mext1) +
      c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PR*mLF6*Sqr(mext1) + 2*c2tmp4*
      CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PR*mLF6*Sqr(mext1) + 3*c1tmp3*
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*mLF4*Sqr(mext2) + c2tmp4*CpcF4cF5S1PR
      *CpcF6F4S2PR*CpF5F6S3PL*mLF4*Sqr(mext2) + 3*c1tmp3*CpcF4cF5S1PL*
      CpcF6F4S2PL*CpF5F6S3PR*mLF4*Sqr(mext2) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PL*
      CpF5F6S3PR*mLF4*Sqr(mext2) + 2*c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PL
      *mLF5*Sqr(mext2) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PL*mLF5*Sqr(
      mext2) + 2*c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*mLF5*Sqr(mext2) +
      c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*mLF5*Sqr(mext2) + c1tmp3*
      CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PL*mLF6*Sqr(mext2) + c1tmp3*CpcF4cF5S1PL
      *CpcF6F4S2PR*CpF5F6S3PR*mLF6*Sqr(mext2) - c1tmp3*CpcF4cF5S1PR*CpcF6F4S2PR
      *CpF5F6S3PL*mLF4*Sqr(mext3) - c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL*
      mLF4*Sqr(mext3) - c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR*mLF4*Sqr(
      mext3) - c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR*mLF4*Sqr(mext3) -
      c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PL*mLF5*Sqr(mext3) - c2tmp4*
      CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6S3PR*mLF5*Sqr(mext3) - c1tmp3*CpcF4cF5S1PR
      *CpcF6F4S2PL*CpF5F6S3PL*mLF6*Sqr(mext3) - c1tmp3*CpcF4cF5S1PL*CpcF6F4S2PR
      *CpF5F6S3PR*mLF6*Sqr(mext3) + c0tmp2*mLF4*(2*(CpcF4cF5S1PR*CpcF6F4S2PL*
      CpF5F6S3PL + CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6S3PR)*mLF4*mLF6 + 2*mLF5*(
      CpcF4cF5S1PL*CpF5F6S3PL*(CpcF6F4S2PR*mLF4 + CpcF6F4S2PL*mLF6) +
      CpcF4cF5S1PR*CpF5F6S3PR*(CpcF6F4S2PL*mLF4 + CpcF6F4S2PR*mLF6)) + (
      CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6S3PL + CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR
      )*(Sqr(mext1) + Sqr(mext2) - Sqr(mext3)) + 2*(CpcF4cF5S1PR*CpcF6F4S2PR*
      CpF5F6S3PL + CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6S3PR)*Sqr(mLF4))));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS3S5S6 coupling Cp[S[3], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g2n2_SSS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS3S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*c0tmp1*
      CpS1cS4cS5*CpS2S4cS6*CpS3S5S6);

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLU4 mass of internal field U[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU4cU5 coupling Cp[S[1], -U[4], -U[5]]
 * @param[in] CpS2cU6U4 coupling Cp[S[2], -U[6], U[4]]
 * @param[in] CpS3U5U6 coupling Cp[S[3], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g3n3_UUU(
   double mext1, double mext2, double mext3,
   double mLU4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU4cU5, const std::complex<double>&
      CpS2cU6U4, const std::complex<double>& CpS3U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*c0tmp1*
      CpS1cU4cU5*CpS2cU6U4*CpS3U5U6);

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpS3S5V6 coupling Cp[S[3], S[5], V[6]][Mom[S[3]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g4n4_SSV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpS3S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS4cS5*
      CpS2S4cV6*CpS3S5V6*(b0tmp1 + c1tmp3*Sqr(mext1) + c2tmp4*Sqr(mext1) -
      c1tmp3*Sqr(mext2) - c2tmp4*Sqr(mext2) - c1tmp3*Sqr(mext3) + c2tmp4*Sqr(
      mext3) + c0tmp2*(-Sqr(mext1) + Sqr(mext3) + Sqr(mLS4))));

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS3S6V5 coupling Cp[S[3], S[6], V[5]][Mom[S[3]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g5n5_SVS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLS6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS3S6V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS6*mLS6, mLV5*mLV5,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS4cV5*
      CpS2S4cS6*CpS3S6V5*(b0tmp1 - c2tmp4*Sqr(mext1) - c0tmp2*Sqr(mext2) +
      c2tmp4*Sqr(mext2) + c0tmp2*Sqr(mext3) - c2tmp4*Sqr(mext3) + c1tmp3*(-Sqr(
      mext1) + Sqr(mext2) + Sqr(mext3)) + c0tmp2*Sqr(mLS4)));

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS3S5S6 coupling Cp[S[3], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g6n6_VSS(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS3S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS5cV4*
      CpS2cS6V4*CpS3S5S6*(b0tmp1 + c1tmp3*Sqr(mext1) + 3*c2tmp4*Sqr(mext1) + 3*
      c1tmp3*Sqr(mext2) + c2tmp4*Sqr(mext2) - c1tmp3*Sqr(mext3) - c2tmp4*Sqr(
      mext3) + c0tmp2*(2*(Sqr(mext1) + Sqr(mext2) - Sqr(mext3)) + Sqr(mLV4))));

   return result;
}

/**
 * @brief Evaluates T1G7N7 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpS3V5V6 coupling Cp[S[3], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g7n7_SVV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpS3V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*CpS1cS4cV5*
      CpS2S4cV6*CpS3V5V6*(2*b0tmp1 - c1tmp3*Sqr(mext1) - 3*c2tmp4*Sqr(mext1) -
      3*c1tmp3*Sqr(mext2) - c2tmp4*Sqr(mext2) + c1tmp3*Sqr(mext3) + c2tmp4*Sqr(
      mext3) + c0tmp2*(Sqr(mext1) + Sqr(mext2) - Sqr(mext3) + 2*Sqr(mLS4))));

   return result;
}

/**
 * @brief Evaluates T1G8N8 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpS3S5V6 coupling Cp[S[3], S[5], V[6]][Mom[S[3]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g8n8_VSV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpS3S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*CpS1cS5cV4*
      CpS2V4cV6*CpS3S5V6*(2*b0tmp1 + 4*c1tmp3*Sqr(mext1) + 7*c2tmp4*Sqr(mext1)
      + 2*c1tmp3*Sqr(mext2) - c2tmp4*Sqr(mext2) - 4*c1tmp3*Sqr(mext3) + c2tmp4*
      Sqr(mext3) + 2*c0tmp2*(3*Sqr(mext1) - Sqr(mext2) + Sqr(mext3) + Sqr(mLV4)
      )));

   return result;
}

/**
 * @brief Evaluates T1G9N9 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS3S6V5 coupling Cp[S[3], S[6], V[5]][Mom[S[3]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g9n9_VVS(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLS6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS3S6V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS6*mLS6, mLV5*mLV5,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*CpS1cV4cV5*
      CpS2cS6V4*CpS3S6V5*(2*b0tmp1 - c1tmp3*Sqr(mext1) + 2*c2tmp4*Sqr(mext1) +
      7*c1tmp3*Sqr(mext2) + 4*c2tmp4*Sqr(mext2) + c1tmp3*Sqr(mext3) - 4*c2tmp4*
      Sqr(mext3) + 2*c0tmp2*(-Sqr(mext1) + 3*Sqr(mext2) + Sqr(mext3) + Sqr(mLV4
      ))));

   return result;
}

/**
 * @brief Evaluates T1G10N10 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpS3V5V6 coupling Cp[S[3], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t1g10n10_VVV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpS3V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,4)*c0tmp1*
      CpS1cV4cV5*CpS2V4cV6*CpS3V5V6);

   return result;
}

/**
 * @brief Evaluates T2G1N11 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 * @param[in] CpS3S4cS5S5 coupling Cp[S[3], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t2g1n11_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpS1S2cS4, const std::complex<double>&
      CpS3S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*a0tmp1*CpS1S2cS4*CpS3S4cS5S5*Den(
      Sqr(mext3),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T2G2N12 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 * @param[in] CpS3S4cV5V5 coupling Cp[S[3], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t2g2n12_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpS1S2cS4, const std::complex<double>&
      CpS3S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S2cS4*CpS3S4cV5V5*Den(Sqr(mext3),
      Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T3G1N13 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 * @param[in] CpS2S4cS5S5 coupling Cp[S[2], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t3g1n13_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpS1S3cS4, const std::complex<double>&
      CpS2S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*a0tmp1*CpS1S3cS4*CpS2S4cS5S5*Den(
      Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T3G2N14 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 * @param[in] CpS2S4cV5V5 coupling Cp[S[2], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t3g2n14_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpS1S3cS4, const std::complex<double>&
      CpS2S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S3cS4*CpS2S4cV5V5*Den(Sqr(mext2),
      Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T4G1N15 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S3S4S5 coupling Cp[S[2], S[3], S[4], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t4g1n15_SS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S3S4S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS4*mLS4, mLS5*mLS5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*b0tmp1*CpS1cS4cS5*CpS2S3S4S5);

   return result;
}

/**
 * @brief Evaluates T4G2N16 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2S3V4V5 coupling Cp[S[2], S[3], V[4], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t4g2n16_VV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2S3V4V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV4*mLV4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1cV4cV5*CpS2S3V4V5*(-2*b0tmp1 +
      finite));

   return result;
}

/**
 * @brief Evaluates T5G1N17 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 * @param[in] CpS1S4cS5S5 coupling Cp[S[1], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t5g1n17_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpS2S3cS4, const std::complex<double>&
      CpS1S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*a0tmp1*CpS1S4cS5S5*CpS2S3cS4*Den(
      Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T5G2N18 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 * @param[in] CpS1S4cV5V5 coupling Cp[S[1], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t5g2n18_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpS2S3cS4, const std::complex<double>&
      CpS1S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S4cV5V5*CpS2S3cS4*Den(Sqr(mext1),
      Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T6G1N19 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS2cS4cS5 coupling Cp[S[2], -S[4], -S[5]]
 * @param[in] CpS1S3S4S5 coupling Cp[S[1], S[3], S[4], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t6g1n19_SS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5,
   const std::complex<double>& CpS2cS4cS5, const std::complex<double>&
      CpS1S3S4S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS4*mLS4, mLS5*mLS5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*b0tmp1*CpS1S3S4S5*CpS2cS4cS5);

   return result;
}

/**
 * @brief Evaluates T6G2N20 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS2cV4cV5 coupling Cp[S[2], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S3V4V5 coupling Cp[S[1], S[3], V[4], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t6g2n20_VV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5,
   const std::complex<double>& CpS2cV4cV5, const std::complex<double>&
      CpS1S3V4V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLV4*mLV4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S3V4V5*CpS2cV4cV5*(-2*b0tmp1 +
      finite));

   return result;
}

/**
 * @brief Evaluates T7G1N21 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS3cS4cS5 coupling Cp[S[3], -S[4], -S[5]]
 * @param[in] CpS1S2S4S5 coupling Cp[S[1], S[2], S[4], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t7g1n21_SS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5,
   const std::complex<double>& CpS3cS4cS5, const std::complex<double>&
      CpS1S2S4S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS4*mLS4, mLS5*mLS5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(-0.5*b0tmp1*CpS1S2S4S5*CpS3cS4cS5);

   return result;
}

/**
 * @brief Evaluates T7G2N22 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS3cV4cV5 coupling Cp[S[3], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S2V4V5 coupling Cp[S[1], S[2], V[4], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t7g2n22_VV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5,
   const std::complex<double>& CpS3cV4cV5, const std::complex<double>&
      CpS1S2V4V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV4*mLV4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(CpS1S2V4V5*CpS3cV4cV5*(-2*b0tmp1 +
      finite));

   return result;
}

/**
 * @brief Evaluates T8G1N23 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S3PL coupling Cp[-F[6], -F[5], S[3]][PL]
 * @param[in] CpcF6cF5S3PR coupling Cp[-F[6], -F[5], S[3]][PR]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g1n23_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6S4PL, const std::complex<double>&
      CpF5F6S4PR, const std::complex<double>& CpcF6cF5S3PL, const std::complex<
      double>& CpcF6cF5S3PR, const std::complex<double>& CpS1S2cS4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1S2cS4*Den
      (Sqr(mext3),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S3PR*CpF5F6S4PL*mLF5 +
      CpcF6cF5S3PL*CpF5F6S4PR*mLF5 + CpcF6cF5S3PL*CpF5F6S4PL*mLF6 +
      CpcF6cF5S3PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S3PR*CpF5F6S4PL + CpcF6cF5S3PL*
      CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext3))));

   return result;
}

/**
 * @brief Evaluates T8G2N24 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 * @param[in] CpS3cS5cS6 coupling Cp[S[3], -S[5], -S[6]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g2n24_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1S2cS4, const std::complex<double>&
      CpS3cS5cS6, const std::complex<double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b0tmp1*
      CpS1S2cS4*CpS3cS5cS6*CpS4S5S6*Den(Sqr(mext3),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T8G3N25 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 * @param[in] CpS3cU6cU5 coupling Cp[S[3], -U[6], -U[5]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g3n25_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpS1S2cS4, const std::complex<double>&
      CpS3cU6cU5, const std::complex<double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*b0tmp1*
      CpS1S2cS4*CpS3cU6cU5*CpS4U5U6*Den(Sqr(mext3),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T8G4N26 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 * @param[in] CpS3cS5cV6 coupling Cp[S[3], -S[5], -V[6]][Mom[S[3]] - Mom[-S[5]]
    ]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g4n26_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1S2cS4, const std::complex<double>&
      CpS3cS5cV6, const std::complex<double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1S2cS4*
      CpS3cS5cV6*CpS4S5V6*(b0tmp2 + Den(Sqr(mext3),Sqr(mIS4))*(a0tmp1 - 2*
      b1tmp3*Sqr(mext3) + b0tmp2*(Sqr(mIS4) + Sqr(mLS5)))));

   return result;
}

/**
 * @brief Evaluates T8G5N27 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S2cS4 coupling Cp[S[1], S[2], -S[4]]
 * @param[in] CpS3cV5cV6 coupling Cp[S[3], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g5n27_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1S2cS4, const std::complex<double>&
      CpS3cV5cV6, const std::complex<double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1S2cS4*
      CpS3cV5cV6*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext3),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T8G6N28 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S3PL coupling Cp[-F[6], -F[5], S[3]][PL]
 * @param[in] CpcF6cF5S3PR coupling Cp[-F[6], -F[5], S[3]][PR]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g6n28_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PR, const std::complex<double>&
      CpF5F6V4PL, const std::complex<double>& CpcF6cF5S3PL, const std::complex<
      double>& CpcF6cF5S3PR, const std::complex<double>& CpS1S2cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1S2cV4*(
      b0tmp1*(CpcF6cF5S3PR*CpF5F6V4PL + CpcF6cF5S3PL*CpF5F6V4PR)*mLF5 + b1tmp2*
      (CpcF6cF5S3PR*CpF5F6V4PL*mLF5 + CpcF6cF5S3PL*CpF5F6V4PR*mLF5 +
      CpcF6cF5S3PL*CpF5F6V4PL*mLF6 + CpcF6cF5S3PR*CpF5F6V4PR*mLF6))*Den(Sqr(
      mext3),Sqr(mIV4))*(Sqr(mext1) - Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T8G7N29 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpS3cS5cS6 coupling Cp[S[3], -S[5], -S[6]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g7n29_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpS3cS5cS6, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*(b0tmp1 + 2
      *b1tmp2)*CpS1S2cV4*CpS3cS5cS6*CpS5S6V4*Den(Sqr(mext3),Sqr(mIV4))*(Sqr(
      mext1) - Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T8G8N30 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpS3cU6cU5 coupling Cp[S[3], -U[6], -U[5]]
 * @param[in] CpU5U6V41 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 * @param[in] CpU5U6V42 coupling Cp[U[5], U[6], V[4]][Mom[U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g8n30_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpS3cU6cU5, const std::complex<double>& CpU5U6V41, const std::complex<
      double>& CpU5U6V42,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLU5*mLU5, mLU6*mLU6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*CpS1S2cV4*
      CpS3cU6cU5*(b1tmp2*(CpU5U6V41 - CpU5U6V42) - b0tmp1*CpU5U6V42)*Den(Sqr(
      mext3),Sqr(mIV4))*(Sqr(mext1) - Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T8G9N31 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpS3cS5cV6 coupling Cp[S[3], -S[5], -V[6]][Mom[S[3]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g9n31_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpS3cS5cV6, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1 -
      b1tmp2)*CpS1S2cV4*CpS3cS5cV6*CpS5V4V6*Den(Sqr(mext3),Sqr(mIV4))*(Sqr(
      mext1) - Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T8G10N32 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpS3cV5cV6 coupling Cp[S[3], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t8g10n32_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpS3cV5cV6, const std::complex<double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1.5)*(b0tmp1 +
      2*b1tmp2)*CpS1S2cV4*CpS3cV5cV6*CpV4V5V6*Den(Sqr(mext3),Sqr(mIV4))*(Sqr(
      mext1) - Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T9G1N33 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S2PL coupling Cp[-F[6], -F[5], S[2]][PL]
 * @param[in] CpcF6cF5S2PR coupling Cp[-F[6], -F[5], S[2]][PR]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g1n33_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6S4PL, const std::complex<double>&
      CpF5F6S4PR, const std::complex<double>& CpcF6cF5S2PL, const std::complex<
      double>& CpcF6cF5S2PR, const std::complex<double>& CpS1S3cS4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1S3cS4*Den
      (Sqr(mext2),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S2PR*CpF5F6S4PL*mLF5 +
      CpcF6cF5S2PL*CpF5F6S4PR*mLF5 + CpcF6cF5S2PL*CpF5F6S4PL*mLF6 +
      CpcF6cF5S2PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S2PR*CpF5F6S4PL + CpcF6cF5S2PL*
      CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext2))));

   return result;
}

/**
 * @brief Evaluates T9G2N34 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 * @param[in] CpS2cS5cS6 coupling Cp[S[2], -S[5], -S[6]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g2n34_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1S3cS4, const std::complex<double>&
      CpS2cS5cS6, const std::complex<double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b0tmp1*
      CpS1S3cS4*CpS2cS5cS6*CpS4S5S6*Den(Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T9G3N35 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 * @param[in] CpS2cU6cU5 coupling Cp[S[2], -U[6], -U[5]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g3n35_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpS1S3cS4, const std::complex<double>&
      CpS2cU6cU5, const std::complex<double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*b0tmp1*
      CpS1S3cS4*CpS2cU6cU5*CpS4U5U6*Den(Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T9G4N36 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 * @param[in] CpS2cS5cV6 coupling Cp[S[2], -S[5], -V[6]][Mom[S[2]] - Mom[-S[5]]
    ]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g4n36_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1S3cS4, const std::complex<double>&
      CpS2cS5cV6, const std::complex<double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1S3cS4*
      CpS2cS5cV6*CpS4S5V6*(b0tmp2 + Den(Sqr(mext2),Sqr(mIS4))*(a0tmp1 - 2*
      b1tmp3*Sqr(mext2) + b0tmp2*(Sqr(mIS4) + Sqr(mLS5)))));

   return result;
}

/**
 * @brief Evaluates T9G5N37 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S3cS4 coupling Cp[S[1], S[3], -S[4]]
 * @param[in] CpS2cV5cV6 coupling Cp[S[2], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g5n37_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1S3cS4, const std::complex<double>&
      CpS2cV5cV6, const std::complex<double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1S3cS4*
      CpS2cV5cV6*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T9G6N38 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S2PL coupling Cp[-F[6], -F[5], S[2]][PL]
 * @param[in] CpcF6cF5S2PR coupling Cp[-F[6], -F[5], S[2]][PR]
 * @param[in] CpS1S3cV4 coupling Cp[S[1], S[3], -V[4]][Mom[S[1]] - Mom[S[3]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g6n38_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PR, const std::complex<double>&
      CpF5F6V4PL, const std::complex<double>& CpcF6cF5S2PL, const std::complex<
      double>& CpcF6cF5S2PR, const std::complex<double>& CpS1S3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1S3cV4*(
      b0tmp1*(CpcF6cF5S2PR*CpF5F6V4PL + CpcF6cF5S2PL*CpF5F6V4PR)*mLF5 + b1tmp2*
      (CpcF6cF5S2PR*CpF5F6V4PL*mLF5 + CpcF6cF5S2PL*CpF5F6V4PR*mLF5 +
      CpcF6cF5S2PL*CpF5F6V4PL*mLF6 + CpcF6cF5S2PR*CpF5F6V4PR*mLF6))*Den(Sqr(
      mext2),Sqr(mIV4))*(Sqr(mext1) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T9G7N39 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1S3cV4 coupling Cp[S[1], S[3], -V[4]][Mom[S[1]] - Mom[S[3]]]
 * @param[in] CpS2cS5cS6 coupling Cp[S[2], -S[5], -S[6]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g7n39_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1S3cV4, const std::complex<double>&
      CpS2cS5cS6, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*(b0tmp1 + 2
      *b1tmp2)*CpS1S3cV4*CpS2cS5cS6*CpS5S6V4*Den(Sqr(mext2),Sqr(mIV4))*(Sqr(
      mext1) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T9G8N40 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1S3cV4 coupling Cp[S[1], S[3], -V[4]][Mom[S[1]] - Mom[S[3]]]
 * @param[in] CpS2cU6cU5 coupling Cp[S[2], -U[6], -U[5]]
 * @param[in] CpU5U6V41 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 * @param[in] CpU5U6V42 coupling Cp[U[5], U[6], V[4]][Mom[U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g8n40_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1S3cV4, const std::complex<double>&
      CpS2cU6cU5, const std::complex<double>& CpU5U6V41, const std::complex<
      double>& CpU5U6V42,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLU5*mLU5, mLU6*mLU6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*CpS1S3cV4*
      CpS2cU6cU5*(b1tmp2*(CpU5U6V41 - CpU5U6V42) - b0tmp1*CpU5U6V42)*Den(Sqr(
      mext2),Sqr(mIV4))*(Sqr(mext1) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T9G9N41 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S3cV4 coupling Cp[S[1], S[3], -V[4]][Mom[S[1]] - Mom[S[3]]]
 * @param[in] CpS2cS5cV6 coupling Cp[S[2], -S[5], -V[6]][Mom[S[2]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g9n41_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1S3cV4, const std::complex<double>&
      CpS2cS5cV6, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1 -
      b1tmp2)*CpS1S3cV4*CpS2cS5cV6*CpS5V4V6*Den(Sqr(mext2),Sqr(mIV4))*(Sqr(
      mext1) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T9G10N42 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S3cV4 coupling Cp[S[1], S[3], -V[4]][Mom[S[1]] - Mom[S[3]]]
 * @param[in] CpS2cV5cV6 coupling Cp[S[2], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t9g10n42_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1S3cV4, const std::complex<double>&
      CpS2cV5cV6, const std::complex<double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1.5)*(b0tmp1 +
      2*b1tmp2)*CpS1S3cV4*CpS2cV5cV6*CpV4V5V6*Den(Sqr(mext2),Sqr(mIV4))*(Sqr(
      mext1) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G1N43 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g1n43_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6S4PL, const std::complex<double>&
      CpF5F6S4PR, const std::complex<double>& CpcF6cF5S1PL, const std::complex<
      double>& CpcF6cF5S1PR, const std::complex<double>& CpS2S3cS4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS2S3cS4*Den
      (Sqr(mext1),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S1PR*CpF5F6S4PL*mLF5 +
      CpcF6cF5S1PL*CpF5F6S4PR*mLF5 + CpcF6cF5S1PL*CpF5F6S4PL*mLF6 +
      CpcF6cF5S1PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S1PR*CpF5F6S4PL + CpcF6cF5S1PL*
      CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext1))));

   return result;
}

/**
 * @brief Evaluates T10G2N44 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g2n44_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cS6, const std::complex<double>&
      CpS2S3cS4, const std::complex<double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b0tmp1*
      CpS1cS5cS6*CpS2S3cS4*CpS4S5S6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G3N45 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g3n45_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU6cU5, const std::complex<double>&
      CpS2S3cS4, const std::complex<double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*b0tmp1*
      CpS1cU6cU5*CpS2S3cS4*CpS4U5U6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G4N46 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g4n46_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV6, const std::complex<double>&
      CpS2S3cS4, const std::complex<double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS5cV6*
      CpS2S3cS4*CpS4S5V6*Den(Sqr(mext1),Sqr(mIS4))*(a0tmp1 - 2*b1tmp3*Sqr(mext1
      ) + b0tmp2*(Sqr(mext1) + Sqr(mLS5))));

   return result;
}

/**
 * @brief Evaluates T10G5N47 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS2S3cS4 coupling Cp[S[2], S[3], -S[4]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g5n47_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV5cV6, const std::complex<double>&
      CpS2S3cS4, const std::complex<double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1cV5cV6*
      CpS2S3cS4*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G6N48 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 * @param[in] CpS2S3cV4 coupling Cp[S[2], S[3], -V[4]][Mom[S[2]] - Mom[S[3]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g6n48_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PR, const std::complex<double>&
      CpF5F6V4PL, const std::complex<double>& CpcF6cF5S1PL, const std::complex<
      double>& CpcF6cF5S1PR, const std::complex<double>& CpS2S3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS2S3cV4*(
      b0tmp1*(CpcF6cF5S1PR*CpF5F6V4PL + CpcF6cF5S1PL*CpF5F6V4PR)*mLF5 + b1tmp2*
      (CpcF6cF5S1PR*CpF5F6V4PL*mLF5 + CpcF6cF5S1PL*CpF5F6V4PR*mLF5 +
      CpcF6cF5S1PL*CpF5F6V4PL*mLF6 + CpcF6cF5S1PR*CpF5F6V4PR*mLF6))*Den(Sqr(
      mext1),Sqr(mIV4))*(Sqr(mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G7N49 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS2S3cV4 coupling Cp[S[2], S[3], -V[4]][Mom[S[2]] - Mom[S[3]]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g7n49_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cS6, const std::complex<double>&
      CpS2S3cV4, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*(b0tmp1 + 2
      *b1tmp2)*CpS1cS5cS6*CpS2S3cV4*CpS5S6V4*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(
      mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G8N50 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpS2S3cV4 coupling Cp[S[2], S[3], -V[4]][Mom[S[2]] - Mom[S[3]]]
 * @param[in] CpU5U6V41 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 * @param[in] CpU5U6V42 coupling Cp[U[5], U[6], V[4]][Mom[U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g8n50_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU6cU5, const std::complex<double>&
      CpS2S3cV4, const std::complex<double>& CpU5U6V41, const std::complex<
      double>& CpU5U6V42,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*CpS1cU6cU5
      *CpS2S3cV4*(b1tmp2*(CpU5U6V41 - CpU5U6V42) - b0tmp1*CpU5U6V42)*Den(Sqr(
      mext1),Sqr(mIV4))*(Sqr(mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G9N51 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2S3cV4 coupling Cp[S[2], S[3], -V[4]][Mom[S[2]] - Mom[S[3]]]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g9n51_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV6, const std::complex<double>&
      CpS2S3cV4, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1 -
      b1tmp2)*CpS1cS5cV6*CpS2S3cV4*CpS5V4V6*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(
      mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G10N52 diagram for process S -> SS
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field S[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS2S3cV4 coupling Cp[S[2], S[3], -V[4]][Mom[S[2]] - Mom[S[3]]]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSS calculate_diagram_SSS_t10g10n52_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV5cV6, const std::complex<double>&
      CpS2S3cV4, const std::complex<double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSS result;

   result.m_decay = mext1;
   result.m_scalar_1 = mext2;
   result.m_scalar_2 = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1.5)*(b0tmp1 +
      2*b1tmp2)*CpS1cV5cV6*CpS2S3cV4*CpV4V5V6*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(
      mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T1G1N1 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F6V3PL coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V3PR coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6F4V2PL coupling Cp[-F[6], F[4], V[2]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpcF6F4V2PR coupling Cp[-F[6], F[4], V[2]][LorentzProduct[gamma[
    lt3], PR]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g1n1_FFF(
   double mext1, double mext2, double mext3,
   double mLF4, double mLF5, double mLF6,
   const std::complex<double>& CpcF4cF5S1PL, const std::complex<double>&
      CpcF4cF5S1PR, const std::complex<double>& CpF5F6V3PL, const std::complex<
      double>& CpF5F6V3PR, const std::complex<double>& CpcF6F4V2PL, const std::
      complex<double>& CpcF6F4V2PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c00tmp3 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c1tmp4 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c12tmp5 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c2tmp6 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c22tmp7 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_eps = oneOver16PiSqr*(-2*(c0tmp2*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL - CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mLF4 +
      c2tmp6*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4 - CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4 + CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5 -
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5) + c1tmp4*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 - CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 -
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6 + CpcF4cF5S1PR*CpcF6F4V2PL*
      CpF5F6V3PR*mLF6)));
   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*(-4*
      c00tmp3*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4 + CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4 + CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5 +
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5) - 2*c0tmp2*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PL*mLF4*mLF5*mLF6 - 2*c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL
      *CpF5F6V3PR*mLF4*mLF5*mLF6 + 2*b0tmp1*(CpcF4cF5S1PL*(CpcF6F4V2PL*
      CpF5F6V3PL*mLF4 + CpcF6F4V2PR*CpF5F6V3PR*mLF5 - CpcF6F4V2PR*CpF5F6V3PL*
      mLF6) + CpcF4cF5S1PR*(CpcF6F4V2PR*CpF5F6V3PR*mLF4 + CpcF6F4V2PL*
      CpF5F6V3PL*mLF5 - CpcF6F4V2PL*CpF5F6V3PR*mLF6)) + 2*c0tmp2*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*Cube(mLF4) + 2*c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PR*
      CpF5F6V3PR*Cube(mLF4) + c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4*
      Sqr(mext1) + c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(mext1) +
      3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(mext1) + c0tmp2*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(mext1) + c1tmp4*CpcF4cF5S1PR
      *CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(mext1) + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(mext1) + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*
      CpF5F6V3PL*mLF5*Sqr(mext1) + c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*
      mLF5*Sqr(mext1) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6*Sqr(
      mext1) - 2*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6*Sqr(mext1) -
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6*Sqr(mext1) - 2*c2tmp6*
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6*Sqr(mext1) + c0tmp2*CpcF4cF5S1PL
      *CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(mext2) + 3*c1tmp4*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(mext2) + c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*mLF4*Sqr(mext2) + c0tmp2*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*
      mLF4*Sqr(mext2) + 3*c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(
      mext2) + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(mext2) + 2*
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5*Sqr(mext2) + c2tmp6*
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5*Sqr(mext2) + 2*c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5*Sqr(mext2) + c2tmp6*CpcF4cF5S1PL
      *CpcF6F4V2PR*CpF5F6V3PR*mLF5*Sqr(mext2) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR
      *CpF5F6V3PL*mLF6*Sqr(mext2) - c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*
      mLF6*Sqr(mext2) - c0tmp2*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(
      mext3) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(mext3) -
      c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4*Sqr(mext3) - c0tmp2*
      CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(mext3) - c1tmp4*CpcF4cF5S1PR
      *CpcF6F4V2PR*CpF5F6V3PR*mLF4*Sqr(mext3) - c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PR
      *CpF5F6V3PR*mLF4*Sqr(mext3) - c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mLF5*Sqr(mext3) - c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5*Sqr(
      mext3) + c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6*Sqr(mext3) +
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6*Sqr(mext3) + 2*c0tmp2*
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5*Sqr(mLF4) + 2*c0tmp2*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5*Sqr(mLF4) - 2*c0tmp2*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6*Sqr(mLF4) - 2*c0tmp2*
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6*Sqr(mLF4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*mLF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 2*
      c22tmp7*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
       + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mLF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mLF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mLF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mLF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6 -
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6));
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*mLF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 2*
      c22tmp7*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
       + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mLF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mLF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mLF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mLF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6 -
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6));
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*mLF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 2*
      c22tmp7*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
       + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mLF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mLF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mLF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mLF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6 -
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6));
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,2)*(c1tmp4*
      CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 2*c22tmp7*CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + 3*c2tmp6*CpcF4cF5S1PL*CpcF6F4V2PL*
      CpF5F6V3PL*mLF4 + c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 2*
      c22tmp7*CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 + 3*c2tmp6*CpcF4cF5S1PR*
      CpcF6F4V2PR*CpF5F6V3PR*mLF4 + c0tmp2*(CpcF4cF5S1PL*CpcF6F4V2PL*CpF5F6V3PL
       + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR)*mLF4 + 2*c22tmp7*CpcF4cF5S1PR*
      CpcF6F4V2PL*CpF5F6V3PL*mLF5 + c2tmp6*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*
      mLF5 + 2*c22tmp7*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + c2tmp6*
      CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PR*mLF5 + 2*c12tmp5*(CpcF4cF5S1PL*
      CpcF6F4V2PL*CpF5F6V3PL*mLF4 + CpcF4cF5S1PR*CpcF6F4V2PR*CpF5F6V3PR*mLF4 +
      CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PL*mLF5 + CpcF4cF5S1PL*CpcF6F4V2PR*
      CpF5F6V3PR*mLF5) - c1tmp4*CpcF4cF5S1PL*CpcF6F4V2PR*CpF5F6V3PL*mLF6 -
      c1tmp4*CpcF4cF5S1PR*CpcF6F4V2PL*CpF5F6V3PR*mLF6));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS4cS6V2 coupling Cp[S[4], -S[6], V[2]][Mom[S[4]] - Mom[-S[6]]]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g2n2_SSS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS4cS6V2, const std::complex<double>& CpS5S6V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c00tmp1 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c12tmp2 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c22tmp4 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,4)*c00tmp1*
      CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,4)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpS1cS4cS5*CpS4cS6V2*CpS5S6V3);

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLU4 mass of internal field U[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU4cU5 coupling Cp[S[1], -U[4], -U[5]]
 * @param[in] CpU5U6V3 coupling Cp[U[5], U[6], V[3]][Mom[U[5]]]
 * @param[in] CpcU6U4V2 coupling Cp[-U[6], U[4], V[2]][Mom[-U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g3n3_UUU(
   double mext1, double mext2, double mext3,
   double mLU4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU4cU5, const std::complex<double>& CpU5U6V3
      , const std::complex<double>& CpcU6U4V2,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c00tmp1 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);
   const auto c12tmp2 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);
   const auto c22tmp4 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*c00tmp1*
      CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,1)*(c12tmp2 +
      c22tmp4 + c2tmp3)*CpcU6U4V2*CpS1cU4cU5*CpU5U6V3);

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS4V2cV6 coupling Cp[S[4], V[2], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g4n4_SSV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS4V2cV6, const std::complex<double>& CpS5V3V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*c0tmp1*
      CpS1cS4cS5*CpS4V2cV6*CpS5V3V6);

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS4cS6V2 coupling Cp[S[4], -S[6], V[2]][Mom[S[4]] - Mom[-S[6]]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g5n5_SVS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLS6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS4cS6V2, const std::complex<double>& CpS6V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c00tmp1 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c12tmp2 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c22tmp4 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,2)*c00tmp1*
      CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,2)*(c12tmp2 +
      c22tmp4 - c2tmp3)*CpS1cS4cV5*CpS4cS6V2*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 * @param[in] CpcS6V2V4 coupling Cp[-S[6], V[2], V[4]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g6n6_VSS(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>& CpS5S6V3
      , const std::complex<double>& CpcS6V2V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c00tmp2 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c12tmp4 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c2tmp5 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c22tmp6 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,2)*c00tmp2*
      CpcS6V2V4*CpS1cS5cV4*CpS5S6V3);
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,2)*(2*c0tmp1
      + c12tmp4 + 2*c1tmp3 + c22tmp6 + 3*c2tmp5)*CpcS6V2V4*CpS1cS5cV4*CpS5S6V3)
      ;

   return result;
}

/**
 * @brief Evaluates T1G7N7 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS4V2cV6 coupling Cp[S[4], V[2], -V[6]][g[lt2, lt3]]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g7n7_SVV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS4V2cV6, const std::complex<double>& CpV3V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c00tmp3 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c1tmp4 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c12tmp5 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c2tmp6 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c22tmp7 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1cS4cV5*
      CpS4V2cV6*CpV3V5V6*(b0tmp1 - c00tmp3 - c1tmp4*Sqr(mext1) - c2tmp6*Sqr(
      mext1) - c0tmp2*Sqr(mext2) + c1tmp4*Sqr(mext2) + c2tmp6*Sqr(mext2) +
      c0tmp2*Sqr(mext3) + c1tmp4*Sqr(mext3) - c2tmp6*Sqr(mext3) + c0tmp2*Sqr(
      mLS4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,-1)*(c12tmp5
      - 4*c1tmp4 + c22tmp7 - c2tmp6)*CpS1cS4cV5*CpS4V2cV6*CpV3V5V6);

   return result;
}

/**
 * @brief Evaluates T1G8N8 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 * @param[in] CpV2V4cV6 coupling Cp[V[2], V[4], -V[6]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[4]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[V[4
    ]] + Mom[-V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g8n8_VSV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>& CpS5V3V6
      , const std::complex<double>& CpV2V4cV6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c00tmp3 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c1tmp4 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c12tmp5 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c2tmp6 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c22tmp7 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1cS5cV4*
      CpS5V3V6*CpV2V4cV6*(b0tmp1 - c00tmp3 + 2*c0tmp2*Sqr(mext1) + c1tmp4*Sqr(
      mext1) + 3*c2tmp6*Sqr(mext1) + 2*c0tmp2*Sqr(mext2) + 3*c1tmp4*Sqr(mext2)
      + c2tmp6*Sqr(mext2) - 2*c0tmp2*Sqr(mext3) - c1tmp4*Sqr(mext3) - c2tmp6*
      Sqr(mext3) + c0tmp2*Sqr(mLV4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp2
       + c12tmp5 - 2*c1tmp4 + c22tmp7 + 3*c2tmp6)*CpS1cS5cV4*CpS5V3V6*CpV2V4cV6
      );

   return result;
}

/**
 * @brief Evaluates T1G9N9 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpcS6V2V4 coupling Cp[-S[6], V[2], V[4]][g[lt2, lt3]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g9n9_VVS(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLS6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpcS6V2V4, const std::complex<double>& CpS6V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*c0tmp1*
      CpcS6V2V4*CpS1cV4cV5*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G10N10 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpV2V4cV6 coupling Cp[V[2], V[4], -V[6]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[4]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[V[4
    ]] + Mom[-V[6]])]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t1g10n10_VVV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpV2V4cV6, const std::complex<double>& CpV3V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c00tmp3 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c1tmp4 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c12tmp5 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c2tmp6 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c22tmp7 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-0.5)*
      CpS1cV4cV5*CpV2V4cV6*CpV3V5V6*(4*b0tmp1 + 20*c00tmp3 - 4*finite - 4*
      c0tmp2*Sqr(mext1) + c1tmp4*Sqr(mext1) + 4*c2tmp6*Sqr(mext1) + 6*c0tmp2*
      Sqr(mext2) + 5*c1tmp4*Sqr(mext2) + 2*c2tmp6*Sqr(mext2) + 4*c0tmp2*Sqr(
      mext3) - c1tmp4*Sqr(mext3) - 2*c2tmp6*Sqr(mext3) + 4*c0tmp2*Sqr(mLV4)));
   result.form_factor_11 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);
   result.form_factor_21 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);
   result.form_factor_12 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);
   result.form_factor_22 = oneOver16PiSqr*(std::complex<double>(0,-1)*(5*c0tmp2
       + c1tmp4 + 10*(c12tmp5 + c22tmp7 + c2tmp6))*CpS1cV4cV5*CpV2V4cV6*
      CpV3V5V6);

   return result;
}

/**
 * @brief Evaluates T2G1N11 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 * @param[in] CpcS5S5V3V4 coupling Cp[-S[5], S[5], V[3], V[4]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t2g1n11_VS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5,
   const std::complex<double>& CpS1V2cV4, const std::complex<double>&
      CpcS5S5V3V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(0.5*a0tmp1*CpcS5S5V3V4*CpS1V2cV4*Den(
      Sqr(mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T2G2N12 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 * @param[in] CpV3V4cV5V51 coupling Cp[V[3], V[4], -V[5], V[5]][g[lt1, lt2] g[
    lt3, lt4]]
 * @param[in] CpV3V4cV5V52 coupling Cp[V[3], V[4], -V[5], V[5]][g[lt1, lt4] g[
    lt2, lt3]]
 * @param[in] CpV3V4cV5V53 coupling Cp[V[3], V[4], -V[5], V[5]][g[lt1, lt3] g[
    lt2, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t2g2n12_VV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5,
   const std::complex<double>& CpS1V2cV4, const std::complex<double>&
      CpV3V4cV5V51, const std::complex<double>& CpV3V4cV5V52, const std::
      complex<double>& CpV3V4cV5V53,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*CpS1V2cV4*Den(Sqr(mext3),Sqr(
      mIV4))*(a0tmp1*(4*CpV3V4cV5V51 + CpV3V4cV5V52 + CpV3V4cV5V53) - 2*
      CpV3V4cV5V51*finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T3G1N13 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpcS5S5V2V4 coupling Cp[-S[5], S[5], V[2], V[4]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t3g1n13_VS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpcS5S5V2V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(0.5*a0tmp1*CpcS5S5V2V4*CpS1V3cV4*Den(
      Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T3G2N14 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpV2V4cV5V51 coupling Cp[V[2], V[4], -V[5], V[5]][g[lt1, lt2] g[
    lt3, lt4]]
 * @param[in] CpV2V4cV5V52 coupling Cp[V[2], V[4], -V[5], V[5]][g[lt1, lt4] g[
    lt2, lt3]]
 * @param[in] CpV2V4cV5V53 coupling Cp[V[2], V[4], -V[5], V[5]][g[lt1, lt3] g[
    lt2, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t3g2n14_VV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpV2V4cV5V51, const std::complex<double>& CpV2V4cV5V52, const std::
      complex<double>& CpV2V4cV5V53,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*CpS1V3cV4*Den(Sqr(mext2),Sqr(
      mIV4))*(a0tmp1*(4*CpV2V4cV5V51 + CpV2V4cV5V52 + CpV2V4cV5V53) - 2*
      CpV2V4cV5V51*finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T4G1N15 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS4S5V2V3 coupling Cp[S[4], S[5], V[2], V[3]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t4g1n15_SS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS4S5V2V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS4*mLS4, mLS5*mLS5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*b0tmp1*CpS1cS4cS5*CpS4S5V2V3);

   return result;
}

/**
 * @brief Evaluates T4G2N16 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpV2V3V4V51 coupling Cp[V[2], V[3], V[4], V[5]][g[lt1, lt2] g[lt3
    , lt4]]
 * @param[in] CpV2V3V4V52 coupling Cp[V[2], V[3], V[4], V[5]][g[lt1, lt4] g[lt2
    , lt3]]
 * @param[in] CpV2V3V4V53 coupling Cp[V[2], V[3], V[4], V[5]][g[lt1, lt3] g[lt2
    , lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t4g2n16_VV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpV2V3V4V51, const std::complex<double>& CpV2V3V4V52, const std::complex<
      double>& CpV2V3V4V53,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV4*mLV4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*b0tmp1*CpS1cV4cV5*(4*CpV2V3V4V51
       + CpV2V3V4V52 + CpV2V3V4V53) + CpS1cV4cV5*CpV2V3V4V51*finite);

   return result;
}

/**
 * @brief Evaluates T5G1N17 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 * @param[in] CpS1S4cS5S5 coupling Cp[S[1], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t5g1n17_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpcS4V2V3, const std::complex<double>&
      CpS1S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(-0.5*a0tmp1*CpcS4V2V3*CpS1S4cS5S5*Den
      (Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T5G2N18 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 * @param[in] CpS1S4cV5V5 coupling Cp[S[1], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t5g2n18_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpcS4V2V3, const std::complex<double>&
      CpS1S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(CpcS4V2V3*CpS1S4cV5V5*Den(Sqr(mext1),
      Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T6G1N19 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpcS4V2cV5 coupling Cp[-S[4], V[2], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S4V3V5 coupling Cp[S[1], S[4], V[3], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t6g1n19_SV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5,
   const std::complex<double>& CpcS4V2cV5, const std::complex<double>&
      CpS1S4V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS4*mLS4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(b0tmp1*CpcS4V2cV5*CpS1S4V3V5);

   return result;
}

/**
 * @brief Evaluates T7G1N20 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpcS4V3cV5 coupling Cp[-S[4], V[3], -V[5]][g[lt2, lt3]]
 * @param[in] CpS1S4V2V5 coupling Cp[S[1], S[4], V[2], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t7g1n20_SV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5,
   const std::complex<double>& CpcS4V3cV5, const std::complex<double>&
      CpS1S4V2V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS4*mLS4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(b0tmp1*CpcS4V3cV5*CpS1S4V2V5);

   return result;
}

/**
 * @brief Evaluates T8G6N26 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6cF5V3PL coupling Cp[-F[6], -F[5], V[3]][LorentzProduct[gamma
    [lt3], PL]]
 * @param[in] CpcF6cF5V3PR coupling Cp[-F[6], -F[5], V[3]][LorentzProduct[gamma
    [lt3], PR]]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t8g6n26_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PL, const std::complex<double>&
      CpF5F6V4PR, const std::complex<double>& CpcF6cF5V3PL, const std::complex<
      double>& CpcF6cF5V3PR, const std::complex<double>& CpS1V2cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b00tmp3 = lib.B00(mext3*mext3, mLF5*mLF5, mLF6*
      mLF6, scale*scale);
   const auto b1tmp4 = lib.B1(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1V2cV4*
      Den(Sqr(mext3),Sqr(mIV4))*(a0tmp1*(CpcF6cF5V3PL*CpF5F6V4PL + CpcF6cF5V3PR
      *CpF5F6V4PR) - 2*b00tmp3*(CpcF6cF5V3PL*CpF5F6V4PL + CpcF6cF5V3PR*
      CpF5F6V4PR) - b0tmp2*CpcF6cF5V3PR*CpF5F6V4PL*mLF5*mLF6 - b0tmp2*
      CpcF6cF5V3PL*CpF5F6V4PR*mLF5*mLF6 + b1tmp4*CpcF6cF5V3PL*CpF5F6V4PL*Sqr(
      mext3) + b1tmp4*CpcF6cF5V3PR*CpF5F6V4PR*Sqr(mext3) + b0tmp2*CpcF6cF5V3PL*
      CpF5F6V4PL*Sqr(mLF5) + b0tmp2*CpcF6cF5V3PR*CpF5F6V4PR*Sqr(mLF5)));

   return result;
}

/**
 * @brief Evaluates T8G7N27 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 * @param[in] CpcS5cS6V3 coupling Cp[-S[5], -S[6], V[3]][Mom[-S[5]] - Mom[-S[6]
    ]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t8g7n27_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1V2cV4, const std::complex<double>&
      CpcS5cS6V3, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b00tmp1 = lib.B00(mext3*mext3, mLS5*mLS5, mLS6*
      mLS6, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-2)*b00tmp1*
      CpcS5cS6V3*CpS1V2cV4*CpS5S6V4*Den(Sqr(mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T8G8N28 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 * @param[in] CpU5U6V4 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 * @param[in] CpcU6cU5V3 coupling Cp[-U[6], -U[5], V[3]][Mom[-U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t8g8n28_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1V2cV4, const std::complex<double>& CpU5U6V4,
      const std::complex<double>& CpcU6cU5V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b00tmp1 = lib.B00(mext3*mext3, mLU5*mLU5, mLU6*
      mLU6, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b00tmp1*
      CpcU6cU5V3*CpS1V2cV4*CpU5U6V4*Den(Sqr(mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T8G9N29 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 * @param[in] CpcS5V3cV6 coupling Cp[-S[5], V[3], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t8g9n29_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1V2cV4, const std::complex<double>&
      CpcS5V3cV6, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*b0tmp1*
      CpcS5V3cV6*CpS1V2cV4*CpS5V4V6*Den(Sqr(mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T8G10N30 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1V2cV4 coupling Cp[S[1], V[2], -V[4]][g[lt2, lt3]]
 * @param[in] CpV3cV5cV6 coupling Cp[V[3], -V[5], -V[6]][g[lt1, lt2] (-Mom[V[3]
    ] + Mom[-V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[
    -V[5]] + Mom[-V[6]])]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t8g10n30_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1V2cV4, const std::complex<double>&
      CpV3cV5cV6, const std::complex<double>& CpV4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b00tmp3 = lib.B00(mext3*mext3, mLV5*mLV5, mLV6*
      mLV6, scale*scale);
   const auto b1tmp4 = lib.B1(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-
      0.16666666666666666)*CpS1V2cV4*CpV3cV5cV6*CpV4V5V6*(15*b0tmp2 + Den(Sqr(
      mext3),Sqr(mIV4))*(6*a0tmp1 + 30*b00tmp3 + 6*b1tmp4*Sqr(mext3) + 2*finite
      *Sqr(mext3) + 15*b0tmp2*Sqr(mIV4) + 6*b0tmp2*Sqr(mLV5) - 6*finite*Sqr(
      mLV5) - 6*finite*Sqr(mLV6))));

   return result;
}

/**
 * @brief Evaluates T9G6N36 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6cF5V2PL coupling Cp[-F[6], -F[5], V[2]][LorentzProduct[gamma
    [lt3], PL]]
 * @param[in] CpcF6cF5V2PR coupling Cp[-F[6], -F[5], V[2]][LorentzProduct[gamma
    [lt3], PR]]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t9g6n36_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PL, const std::complex<double>&
      CpF5F6V4PR, const std::complex<double>& CpcF6cF5V2PL, const std::complex<
      double>& CpcF6cF5V2PR, const std::complex<double>& CpS1V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b00tmp3 = lib.B00(mext2*mext2, mLF5*mLF5, mLF6*
      mLF6, scale*scale);
   const auto b1tmp4 = lib.B1(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1V3cV4*
      Den(Sqr(mext2),Sqr(mIV4))*(a0tmp1*(CpcF6cF5V2PL*CpF5F6V4PL + CpcF6cF5V2PR
      *CpF5F6V4PR) - 2*b00tmp3*(CpcF6cF5V2PL*CpF5F6V4PL + CpcF6cF5V2PR*
      CpF5F6V4PR) - b0tmp2*CpcF6cF5V2PR*CpF5F6V4PL*mLF5*mLF6 - b0tmp2*
      CpcF6cF5V2PL*CpF5F6V4PR*mLF5*mLF6 + b1tmp4*CpcF6cF5V2PL*CpF5F6V4PL*Sqr(
      mext2) + b1tmp4*CpcF6cF5V2PR*CpF5F6V4PR*Sqr(mext2) + b0tmp2*CpcF6cF5V2PL*
      CpF5F6V4PL*Sqr(mLF5) + b0tmp2*CpcF6cF5V2PR*CpF5F6V4PR*Sqr(mLF5)));

   return result;
}

/**
 * @brief Evaluates T9G7N37 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpcS5cS6V2 coupling Cp[-S[5], -S[6], V[2]][Mom[-S[5]] - Mom[-S[6]
    ]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t9g7n37_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpcS5cS6V2, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b00tmp1 = lib.B00(mext2*mext2, mLS5*mLS5, mLS6*
      mLS6, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-2)*b00tmp1*
      CpcS5cS6V2*CpS1V3cV4*CpS5S6V4*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G8N38 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpU5U6V4 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 * @param[in] CpcU6cU5V2 coupling Cp[-U[6], -U[5], V[2]][Mom[-U[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t9g8n38_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>& CpU5U6V4,
      const std::complex<double>& CpcU6cU5V2,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b00tmp1 = lib.B00(mext2*mext2, mLU5*mLU5, mLU6*
      mLU6, scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b00tmp1*
      CpcU6cU5V2*CpS1V3cV4*CpU5U6V4*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G9N39 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpcS5V2cV6 coupling Cp[-S[5], V[2], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t9g9n39_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpcS5V2cV6, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*b0tmp1*
      CpcS5V2cV6*CpS1V3cV4*CpS5V4V6*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G10N40 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpV2cV5cV6 coupling Cp[V[2], -V[5], -V[6]][g[lt1, lt2] (-Mom[V[2]
    ] + Mom[-V[5]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[
    -V[5]] + Mom[-V[6]])]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t9g10n40_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpV2cV5cV6, const std::complex<double>& CpV4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b00tmp3 = lib.B00(mext2*mext2, mLV5*mLV5, mLV6*
      mLV6, scale*scale);
   const auto b1tmp4 = lib.B1(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-
      0.16666666666666666)*CpS1V3cV4*CpV2cV5cV6*CpV4V5V6*(15*b0tmp2 + Den(Sqr(
      mext2),Sqr(mIV4))*(6*a0tmp1 + 30*b00tmp3 + 6*b1tmp4*Sqr(mext2) + 2*finite
      *Sqr(mext2) + 15*b0tmp2*Sqr(mIV4) + 6*b0tmp2*Sqr(mLV5) - 6*finite*Sqr(
      mLV5) - 6*finite*Sqr(mLV6))));

   return result;
}

/**
 * @brief Evaluates T10G1N41 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g1n41_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6S4PL, const std::complex<double>&
      CpF5F6S4PR, const std::complex<double>& CpcF6cF5S1PL, const std::complex<
      double>& CpcF6cF5S1PR, const std::complex<double>& CpcS4V2V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpcS4V2V3*
      Den(Sqr(mext1),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S1PR*CpF5F6S4PL*mLF5 +
      CpcF6cF5S1PL*CpF5F6S4PR*mLF5 + CpcF6cF5S1PL*CpF5F6S4PL*mLF6 +
      CpcF6cF5S1PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S1PR*CpF5F6S4PL + CpcF6cF5S1PL*
      CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext1))));

   return result;
}

/**
 * @brief Evaluates T10G2N42 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g2n42_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cS6, const std::complex<double>&
      CpcS4V2V3, const std::complex<double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b0tmp1*
      CpcS4V2V3*CpS1cS5cS6*CpS4S5S6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G3N43 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g3n43_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU6cU5, const std::complex<double>&
      CpcS4V2V3, const std::complex<double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,0.5)*b0tmp1*
      CpcS4V2V3*CpS1cU6cU5*CpS4U5U6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G4N44 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g4n44_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV6, const std::complex<double>&
      CpcS4V2V3, const std::complex<double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*CpcS4V2V3*
      CpS1cS5cV6*CpS4S5V6*Den(Sqr(mext1),Sqr(mIS4))*(a0tmp1 - 2*b1tmp3*Sqr(
      mext1) + b0tmp2*(Sqr(mext1) + Sqr(mLS5))));

   return result;
}

/**
 * @brief Evaluates T10G5N45 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpcS4V2V3 coupling Cp[-S[4], V[2], V[3]][g[lt2, lt3]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g5n45_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV5cV6, const std::complex<double>&
      CpcS4V2V3, const std::complex<double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1)*CpcS4V2V3*
      CpS1cV5cV6*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G6N46 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 * @param[in] CpV2V3cV4 coupling Cp[V[2], V[3], -V[4]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[3]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[4]]) + g[lt2, lt3] (-Mom[V[3
    ]] + Mom[-V[4]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g6n46_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PR, const std::complex<double>&
      CpF5F6V4PL, const std::complex<double>& CpcF6cF5S1PL, const std::complex<
      double>& CpcF6cF5S1PR, const std::complex<double>& CpV2V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*CpV2V3cV4*
      (b0tmp1*(CpcF6cF5S1PR*CpF5F6V4PL + CpcF6cF5S1PL*CpF5F6V4PR)*mLF5 + b1tmp2
      *(CpcF6cF5S1PR*CpF5F6V4PL*mLF5 + CpcF6cF5S1PL*CpF5F6V4PR*mLF5 +
      CpcF6cF5S1PL*CpF5F6V4PL*mLF6 + CpcF6cF5S1PR*CpF5F6V4PR*mLF6))*Den(Sqr(
      mext1),Sqr(mIV4))*(Sqr(mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G7N47 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 * @param[in] CpV2V3cV4 coupling Cp[V[2], V[3], -V[4]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[3]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[4]]) + g[lt2, lt3] (-Mom[V[3
    ]] + Mom[-V[4]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g7n47_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cS6, const std::complex<double>& CpS5S6V4
      , const std::complex<double>& CpV2V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-0.5)*(b0tmp1
      + 2*b1tmp2)*CpS1cS5cS6*CpS5S6V4*CpV2V3cV4*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(
      mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G8N48 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpU5U6V4 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 * @param[in] CpV2V3cV4 coupling Cp[V[2], V[3], -V[4]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[3]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[4]]) + g[lt2, lt3] (-Mom[V[3
    ]] + Mom[-V[4]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g8n48_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU6cU5, const std::complex<double>& CpU5U6V4
      , const std::complex<double>& CpV2V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b1tmp1 = lib.B1(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,0.5)*b1tmp1*
      CpS1cU6cU5*CpU5U6V4*CpV2V3cV4*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(mext2) - Sqr
      (mext3)));

   return result;
}

/**
 * @brief Evaluates T10G9N49 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 * @param[in] CpV2V3cV4 coupling Cp[V[2], V[3], -V[4]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[3]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[4]]) + g[lt2, lt3] (-Mom[V[3
    ]] + Mom[-V[4]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g9n49_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV6, const std::complex<double>& CpS5V4V6
      , const std::complex<double>& CpV2V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,-1)*(b0tmp1 -
      b1tmp2)*CpS1cS5cV6*CpS5V4V6*CpV2V3cV4*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(
      mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T10G10N50 diagram for process S -> VV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field V[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpV2V3cV4 coupling Cp[V[2], V[3], -V[4]][g[lt1, lt2] (-Mom[V[2]]
    + Mom[V[3]]) + g[lt1, lt3] (Mom[V[2]] - Mom[-V[4]]) + g[lt2, lt3] (-Mom[V[3
    ]] + Mom[-V[4]])]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SVV calculate_diagram_SVV_t10g10n50_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV5cV6, const std::complex<double>&
      CpV2V3cV4, const std::complex<double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SVV result;

   result.m_decay = mext1;
   result.m_vector_1 = mext2;
   result.m_vector_2 = mext3;

   result.form_factor_g = oneOver16PiSqr*(std::complex<double>(0,1.5)*(b0tmp1 +
      2*b1tmp2)*CpS1cV5cV6*CpV2V3cV4*CpV4V5V6*Den(Sqr(mext1),Sqr(mIV4))*(Sqr(
      mext2) - Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T1G1N1 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F6V3PL coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V3PR coupling Cp[F[5], F[6], V[3]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6F4S2PL coupling Cp[-F[6], F[4], S[2]][PL]
 * @param[in] CpcF6F4S2PR coupling Cp[-F[6], F[4], S[2]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g1n1_FFF(
   double mext1, double mext2, double mext3,
   double mLF4, double mLF5, double mLF6,
   const std::complex<double>& CpcF4cF5S1PL, const std::complex<double>&
      CpcF4cF5S1PR, const std::complex<double>& CpF5F6V3PL, const std::complex<
      double>& CpF5F6V3PR, const std::complex<double>& CpcF6F4S2PL, const std::
      complex<double>& CpcF6F4S2PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLF4*mLF4, mLF6*mLF6, mLF5*mLF5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*(b0tmp1*
      CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6V3PL + b0tmp1*CpcF4cF5S1PR*CpcF6F4S2PL*
      CpF5F6V3PR + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PR*CpF5F6V3PL*mLF4*mLF5 +
      c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6V3PR*mLF4*mLF5 + c2tmp4*
      CpcF4cF5S1PL*CpcF6F4S2PL*CpF5F6V3PL*mLF4*mLF6 + c2tmp4*CpcF4cF5S1PR*
      CpcF6F4S2PR*CpF5F6V3PR*mLF4*mLF6 + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*
      CpF5F6V3PL*mLF5*mLF6 + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6V3PR*mLF5*
      mLF6 + c0tmp2*mLF4*(CpcF4cF5S1PL*(2*CpcF6F4S2PR*CpF5F6V3PL*mLF4 +
      CpcF6F4S2PL*CpF5F6V3PR*mLF5 + CpcF6F4S2PL*CpF5F6V3PL*mLF6) + CpcF4cF5S1PR
      *(2*CpcF6F4S2PL*CpF5F6V3PR*mLF4 + CpcF6F4S2PR*CpF5F6V3PL*mLF5 +
      CpcF6F4S2PR*CpF5F6V3PR*mLF6)) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*
      CpF5F6V3PL*Sqr(mext1) + c2tmp4*CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6V3PR*Sqr(
      mext1) + c2tmp4*CpcF4cF5S1PL*CpcF6F4S2PR*CpF5F6V3PL*Sqr(mLF4) + c2tmp4*
      CpcF4cF5S1PR*CpcF6F4S2PL*CpF5F6V3PR*Sqr(mLF4) + c1tmp3*(CpcF4cF5S1PL*
      CpcF6F4S2PL*mLF4*(CpF5F6V3PR*mLF5 + CpF5F6V3PL*mLF6) + CpcF4cF5S1PR*
      CpcF6F4S2PR*mLF4*(CpF5F6V3PL*mLF5 + CpF5F6V3PR*mLF6) + CpcF4cF5S1PL*
      CpcF6F4S2PR*(CpF5F6V3PR*mLF5*mLF6 + CpF5F6V3PL*(Sqr(mext2) + Sqr(mLF4)))
      + CpcF4cF5S1PR*CpcF6F4S2PL*(CpF5F6V3PL*mLF5*mLF6 + CpF5F6V3PR*(Sqr(mext2)
      + Sqr(mLF4))))));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g2n2_SSS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS5S6V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*(c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cS4cS5*CpS2S4cS6*CpS5S6V3);

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLU4 mass of internal field U[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU4cU5 coupling Cp[S[1], -U[4], -U[5]]
 * @param[in] CpS2cU6U4 coupling Cp[S[2], -U[6], U[4]]
 * @param[in] CpU5U6V3 coupling Cp[U[5], U[6], V[3]][Mom[U[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g3n3_UUU(
   double mext1, double mext2, double mext3,
   double mLU4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU4cU5, const std::complex<double>&
      CpS2cU6U4, const std::complex<double>& CpU5U6V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLU4*mLU4, mLU6*mLU6, mLU5*mLU5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cU4cU5*CpS2cU6U4*CpU5U6V3);

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g4n4_SSV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS4cS5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpS5V3V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*(c0tmp1 -
      c1tmp2 - c2tmp3)*CpS1cS4cS5*CpS2S4cV6*CpS5V3V6);

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cS6 coupling Cp[S[2], S[4], -S[6]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g5n5_SVS(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLS6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cS6, const std::complex<double>& CpS6V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLS6*mLS6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(c0tmp1 -
      c1tmp2 - c2tmp3)*CpS1cS4cV5*CpS2S4cS6*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS5S6V3 coupling Cp[S[5], S[6], V[3]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g6n6_VSS(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS5S6V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c00tmp2 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c11tmp4 = lib.C11(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c12tmp5 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c2tmp6 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);
   const auto c22tmp7 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*CpS1cS5cV4*
      CpS2cS6V4*CpS5S6V3*(4*c00tmp2 + c11tmp4*Sqr(mext1) + 4*c12tmp5*Sqr(mext1)
      + 3*c1tmp3*Sqr(mext1) + 3*c22tmp7*Sqr(mext1) + 5*c2tmp6*Sqr(mext1) + 3*
      c11tmp4*Sqr(mext2) + 4*c12tmp5*Sqr(mext2) + 5*c1tmp3*Sqr(mext2) + c22tmp7
      *Sqr(mext2) + 3*c2tmp6*Sqr(mext2) - c11tmp4*Sqr(mext3) - 2*c12tmp5*Sqr(
      mext3) - 3*c1tmp3*Sqr(mext3) - c22tmp7*Sqr(mext3) - 3*c2tmp6*Sqr(mext3) +
      c1tmp3*Sqr(mLV4) + c2tmp6*Sqr(mLV4) + c0tmp1*(2*(Sqr(mext1) + Sqr(mext2)
      - Sqr(mext3)) + Sqr(mLV4))));

   return result;
}

/**
 * @brief Evaluates T1G7N7 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4cV6 coupling Cp[S[2], S[4], -V[6]][Mom[S[2]] - Mom[S[4]]]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g7n7_SVV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4cV6, const std::complex<double>& CpV3V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c00tmp3 = lib.C00(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c1tmp4 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c11tmp5 = lib.C11(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c12tmp6 = lib.C12(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c2tmp7 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c22tmp8 = lib.C22(mext2*mext2, mext3*mext3, mext1*
      mext1, mLS4*mLS4, mLV6*mLV6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS1cS4cV5*
      CpS2S4cV6*CpV3V5V6*(4*b0tmp1 - 4*c00tmp3 - c11tmp5*Sqr(mext1) - 4*c12tmp6
      *Sqr(mext1) + c1tmp4*Sqr(mext1) - 3*c22tmp8*Sqr(mext1) - c2tmp7*Sqr(mext1
      ) - 3*c11tmp5*Sqr(mext2) - 4*c12tmp6*Sqr(mext2) - c1tmp4*Sqr(mext2) -
      c22tmp8*Sqr(mext2) + c2tmp7*Sqr(mext2) + c0tmp2*Sqr(mext3) + c11tmp5*Sqr(
      mext3) + 2*c12tmp6*Sqr(mext3) - 2*c1tmp4*Sqr(mext3) + c22tmp8*Sqr(mext3)
      - 2*c2tmp7*Sqr(mext3) + 4*c0tmp2*Sqr(mLS4)));

   return result;
}

/**
 * @brief Evaluates T1G8N8 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V3V6 coupling Cp[S[5], V[3], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g8n8_VSV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV4, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpS5V3V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*(2*c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cS5cV4*CpS2V4cV6*CpS5V3V6);

   return result;
}

/**
 * @brief Evaluates T1G9N9 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2cS6V4 coupling Cp[S[2], -S[6], V[4]][Mom[S[2]] - Mom[-S[6]]]
 * @param[in] CpS6V3V5 coupling Cp[S[6], V[3], V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g9n9_VVS(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLS6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2cS6V4, const std::complex<double>& CpS6V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLS6*mLS6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(2*c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cV4cV5*CpS2cS6V4*CpS6V3V5);

   return result;
}

/**
 * @brief Evaluates T1G10N10 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 * @param[in] CpS2V4cV6 coupling Cp[S[2], V[4], -V[6]][g[lt2, lt3]]
 * @param[in] CpV3V5V6 coupling Cp[V[3], V[5], V[6]][g[lt1, lt2] (-Mom[V[3]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t1g10n10_VVV(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV4cV5, const std::complex<double>&
      CpS2V4cV6, const std::complex<double>& CpV3V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext3*mext3, mext1*
      mext1, mLV4*mLV4, mLV6*mLV6, mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-6)*(c0tmp1 +
      c1tmp2 + c2tmp3)*CpS1cV4cV5*CpS2V4cV6*CpV3V5V6);

   return result;
}

/**
 * @brief Evaluates T2G1N11 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpcS5S5V3V4 coupling Cp[-S[5], S[5], V[3], V[4]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t2g1n11_VS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpcS5S5V3V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(a0tmp1*CpcS5S5V3V4*CpS1S2cV4*Den(Sqr(
      mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T2G2N12 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpV3V4cV5V51 coupling Cp[V[3], V[4], -V[5], V[5]][g[lt1, lt2] g[
    lt3, lt4]]
 * @param[in] CpV3V4cV5V52 coupling Cp[V[3], V[4], -V[5], V[5]][g[lt1, lt4] g[
    lt2, lt3]]
 * @param[in] CpV3V4cV5V53 coupling Cp[V[3], V[4], -V[5], V[5]][g[lt1, lt3] g[
    lt2, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t2g2n12_VV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpV3V4cV5V51, const std::complex<double>& CpV3V4cV5V52, const std::
      complex<double>& CpV3V4cV5V53,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(-(CpS1S2cV4*Den(Sqr(mext3),Sqr(mIV4))*(
      a0tmp1*(4*CpV3V4cV5V51 + CpV3V4cV5V52 + CpV3V4cV5V53) - 2*CpV3V4cV5V51*
      finite*Sqr(mLV5))));

   return result;
}

/**
 * @brief Evaluates T3G1N13 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 * @param[in] CpS2S4cS5S5 coupling Cp[S[2], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t3g1n13_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpS1cS4V3, const std::complex<double>&
      CpS2S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(-(a0tmp1*CpS1cS4V3*CpS2S4cS5S5*Den(Sqr(
      mext2),Sqr(mIS4))));

   return result;
}

/**
 * @brief Evaluates T3G2N14 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 * @param[in] CpS2S4cV5V5 coupling Cp[S[2], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t3g2n14_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpS1cS4V3, const std::complex<double>&
      CpS2S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(2*CpS1cS4V3*CpS2S4cV5V5*Den(Sqr(mext2),
      Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T4G1N15 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 * @param[in] CpS2S4V3V5 coupling Cp[S[2], S[4], V[3], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t4g1n15_SV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5,
   const std::complex<double>& CpS1cS4cV5, const std::complex<double>&
      CpS2S4V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS4*mLS4, mLV5*mLV5,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS4*mLS4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*((b0tmp1 - b1tmp2)*CpS1cS4cV5*CpS2S4V3V5
      );

   return result;
}

/**
 * @brief Evaluates T5G1N16 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 * @param[in] CpS1S4cS5S5 coupling Cp[S[1], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t5g1n16_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpS2cS4V3, const std::complex<double>&
      CpS1S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(a0tmp1*CpS1S4cS5S5*CpS2cS4V3*Den(Sqr(
      mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T5G2N17 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 * @param[in] CpS1S4cV5V5 coupling Cp[S[1], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t5g2n17_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpS2cS4V3, const std::complex<double>&
      CpS1S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(-2*CpS1S4cV5V5*CpS2cS4V3*Den(Sqr(mext1)
      ,Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5)));

   return result;
}

/**
 * @brief Evaluates T6G1N18 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpS2cS4cV5 coupling Cp[S[2], -S[4], -V[5]][Mom[S[2]] - Mom[-S[4]]
    ]
 * @param[in] CpS1S4V3V5 coupling Cp[S[1], S[4], V[3], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t6g1n18_SV(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5,
   const std::complex<double>& CpS2cS4cV5, const std::complex<double>&
      CpS1S4V3V5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS4*mLS4, mLV5*mLV5,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLS4*mLS4, mLV5*mLV5,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*((-b0tmp1 + b1tmp2)*CpS1S4V3V5*
      CpS2cS4cV5);

   return result;
}

/**
 * @brief Evaluates T8G6N26 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF6cF5V3PL coupling Cp[-F[6], -F[5], V[3]][LorentzProduct[gamma
    [lt3], PL]]
 * @param[in] CpcF6cF5V3PR coupling Cp[-F[6], -F[5], V[3]][LorentzProduct[gamma
    [lt3], PR]]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t8g6n26_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PL, const std::complex<double>&
      CpF5F6V4PR, const std::complex<double>& CpcF6cF5V3PL, const std::complex<
      double>& CpcF6cF5V3PR, const std::complex<double>& CpS1S2cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b00tmp3 = lib.B00(mext3*mext3, mLF5*mLF5, mLF6*
      mLF6, scale*scale);
   const auto b1tmp4 = lib.B1(mext3*mext3, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,2)*CpS1S2cV4*Den
      (Sqr(mext3),Sqr(mIV4))*(a0tmp1*(CpcF6cF5V3PL*CpF5F6V4PL + CpcF6cF5V3PR*
      CpF5F6V4PR) - 2*b00tmp3*(CpcF6cF5V3PL*CpF5F6V4PL + CpcF6cF5V3PR*
      CpF5F6V4PR) - b0tmp2*CpcF6cF5V3PR*CpF5F6V4PL*mLF5*mLF6 - b0tmp2*
      CpcF6cF5V3PL*CpF5F6V4PR*mLF5*mLF6 + b1tmp4*CpcF6cF5V3PL*CpF5F6V4PL*Sqr(
      mext3) + b1tmp4*CpcF6cF5V3PR*CpF5F6V4PR*Sqr(mext3) + b0tmp2*CpcF6cF5V3PL*
      CpF5F6V4PL*Sqr(mLF5) + b0tmp2*CpcF6cF5V3PR*CpF5F6V4PR*Sqr(mLF5)));

   return result;
}

/**
 * @brief Evaluates T8G7N27 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpcS5cS6V3 coupling Cp[-S[5], -S[6], V[3]][Mom[-S[5]] - Mom[-S[6]
    ]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t8g7n27_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpcS5cS6V3, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b00tmp1 = lib.B00(mext3*mext3, mLS5*mLS5, mLS6*
      mLS6, scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-4)*b00tmp1*
      CpcS5cS6V3*CpS1S2cV4*CpS5S6V4*Den(Sqr(mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T8G9N29 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpcS5V3cV6 coupling Cp[-S[5], V[3], -V[6]][g[lt2, lt3]]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t8g9n29_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpcS5V3cV6, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*b0tmp1*
      CpcS5V3cV6*CpS1S2cV4*CpS5V4V6*Den(Sqr(mext3),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T8G10N30 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1S2cV4 coupling Cp[S[1], S[2], -V[4]][Mom[S[1]] - Mom[S[2]]]
 * @param[in] CpV3cV5cV6 coupling Cp[V[3], -V[5], -V[6]][g[lt1, lt2] (-Mom[V[3]
    ] + Mom[-V[5]]) + g[lt1, lt3] (Mom[V[3]] - Mom[-V[6]]) + g[lt2, lt3] (-Mom[
    -V[5]] + Mom[-V[6]])]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t8g10n30_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1S2cV4, const std::complex<double>&
      CpV3cV5cV6, const std::complex<double>& CpV4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b00tmp3 = lib.B00(mext3*mext3, mLV5*mLV5, mLV6*
      mLV6, scale*scale);
   const auto b1tmp4 = lib.B1(mext3*mext3, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-
      0.3333333333333333)*CpS1S2cV4*CpV3cV5cV6*CpV4V5V6*(15*b0tmp2 + Den(Sqr(
      mext3),Sqr(mIV4))*(6*a0tmp1 + 30*b00tmp3 + 6*b1tmp4*Sqr(mext3) + 2*finite
      *Sqr(mext3) + 15*b0tmp2*Sqr(mIV4) + 6*b0tmp2*Sqr(mLV5) - 6*finite*Sqr(
      mLV5) - 6*finite*Sqr(mLV6))));

   return result;
}

/**
 * @brief Evaluates T9G1N31 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S2PL coupling Cp[-F[6], -F[5], S[2]][PL]
 * @param[in] CpcF6cF5S2PR coupling Cp[-F[6], -F[5], S[2]][PR]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g1n31_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6S4PL, const std::complex<double>&
      CpF5F6S4PR, const std::complex<double>& CpcF6cF5S2PL, const std::complex<
      double>& CpcF6cF5S2PR, const std::complex<double>& CpS1cS4V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,2)*CpS1cS4V3*Den
      (Sqr(mext2),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S2PR*CpF5F6S4PL*mLF5 +
      CpcF6cF5S2PL*CpF5F6S4PR*mLF5 + CpcF6cF5S2PL*CpF5F6S4PL*mLF6 +
      CpcF6cF5S2PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S2PR*CpF5F6S4PL + CpcF6cF5S2PL*
      CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext2))));

   return result;
}

/**
 * @brief Evaluates T9G2N32 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 * @param[in] CpS2cS5cS6 coupling Cp[S[2], -S[5], -S[6]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g2n32_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS4V3, const std::complex<double>&
      CpS2cS5cS6, const std::complex<double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*b0tmp1*
      CpS1cS4V3*CpS2cS5cS6*CpS4S5S6*Den(Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T9G3N33 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 * @param[in] CpS2cU6cU5 coupling Cp[S[2], -U[6], -U[5]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g3n33_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cS4V3, const std::complex<double>&
      CpS2cU6cU5, const std::complex<double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*b0tmp1*
      CpS1cS4V3*CpS2cU6cU5*CpS4U5U6*Den(Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T9G4N34 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 * @param[in] CpS2cS5cV6 coupling Cp[S[2], -S[5], -V[6]][Mom[S[2]] - Mom[-S[5]]
    ]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g4n34_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS4V3, const std::complex<double>&
      CpS2cS5cV6, const std::complex<double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*CpS1cS4V3*
      CpS2cS5cV6*CpS4S5V6*(b0tmp2 + Den(Sqr(mext2),Sqr(mIS4))*(a0tmp1 - 2*
      b1tmp3*Sqr(mext2) + b0tmp2*(Sqr(mIS4) + Sqr(mLS5)))));

   return result;
}

/**
 * @brief Evaluates T9G5N35 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS4V3 coupling Cp[S[1], -S[4], V[3]][Mom[S[1]] - Mom[-S[4]]]
 * @param[in] CpS2cV5cV6 coupling Cp[S[2], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g5n35_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cS4V3, const std::complex<double>&
      CpS2cV5cV6, const std::complex<double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,2)*CpS1cS4V3*
      CpS2cV5cV6*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext2),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T9G6N36 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S2PL coupling Cp[-F[6], -F[5], S[2]][PL]
 * @param[in] CpcF6cF5S2PR coupling Cp[-F[6], -F[5], S[2]][PR]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g6n36_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PR, const std::complex<double>&
      CpF5F6V4PL, const std::complex<double>& CpcF6cF5S2PL, const std::complex<
      double>& CpcF6cF5S2PR, const std::complex<double>& CpS1V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*CpS1V3cV4*(
      b0tmp1*(CpcF6cF5S2PR*CpF5F6V4PL + CpcF6cF5S2PL*CpF5F6V4PR)*mLF5 + b1tmp2*
      (CpcF6cF5S2PR*CpF5F6V4PL*mLF5 + CpcF6cF5S2PL*CpF5F6V4PR*mLF5 +
      CpcF6cF5S2PL*CpF5F6V4PL*mLF6 + CpcF6cF5S2PR*CpF5F6V4PR*mLF6))*Den(Sqr(
      mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G7N37 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpS2cS5cS6 coupling Cp[S[2], -S[5], -S[6]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g7n37_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpS2cS5cS6, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*(b0tmp1 + 2
      *b1tmp2)*CpS1V3cV4*CpS2cS5cS6*CpS5S6V4*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G8N38 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpS2cU6cU5 coupling Cp[S[2], -U[6], -U[5]]
 * @param[in] CpU5U6V4 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g8n38_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpS2cU6cU5, const std::complex<double>& CpU5U6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b1tmp1 = lib.B1(mext2*mext2, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*b1tmp1*
      CpS1V3cV4*CpS2cU6cU5*CpU5U6V4*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G9N39 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpS2cS5cV6 coupling Cp[S[2], -S[5], -V[6]][Mom[S[2]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g9n39_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpS2cS5cV6, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1 -
      b1tmp2)*CpS1V3cV4*CpS2cS5cV6*CpS5V4V6*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T9G10N40 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1V3cV4 coupling Cp[S[1], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpS2cV5cV6 coupling Cp[S[2], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t9g10n40_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1V3cV4, const std::complex<double>&
      CpS2cV5cV6, const std::complex<double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1.5)*(b0tmp1 +
      2*b1tmp2)*CpS1V3cV4*CpS2cV5cV6*CpV4V5V6*Den(Sqr(mext2),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T10G1N41 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g1n41_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6S4PL, const std::complex<double>&
      CpF5F6S4PR, const std::complex<double>& CpcF6cF5S1PL, const std::complex<
      double>& CpcF6cF5S1PR, const std::complex<double>& CpS2cS4V3,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*CpS2cS4V3*
      Den(Sqr(mext1),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S1PR*CpF5F6S4PL*mLF5 +
      CpcF6cF5S1PL*CpF5F6S4PR*mLF5 + CpcF6cF5S1PL*CpF5F6S4PL*mLF6 +
      CpcF6cF5S1PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S1PR*CpF5F6S4PL + CpcF6cF5S1PL*
      CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext1))));

   return result;
}

/**
 * @brief Evaluates T10G2N42 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g2n42_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cS6, const std::complex<double>&
      CpS2cS4V3, const std::complex<double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1)*b0tmp1*
      CpS1cS5cS6*CpS2cS4V3*CpS4S5S6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G3N43 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g3n43_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU6cU5, const std::complex<double>&
      CpS2cS4V3, const std::complex<double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*b0tmp1*
      CpS1cU6cU5*CpS2cS4V3*CpS4U5U6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G4N44 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g4n44_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV6, const std::complex<double>&
      CpS2cS4V3, const std::complex<double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,2)*CpS1cS5cV6*
      CpS2cS4V3*CpS4S5V6*Den(Sqr(mext1),Sqr(mIS4))*(a0tmp1 - 2*b1tmp3*Sqr(mext1
      ) + b0tmp2*(Sqr(mext1) + Sqr(mLS5))));

   return result;
}

/**
 * @brief Evaluates T10G5N45 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS2cS4V3 coupling Cp[S[2], -S[4], V[3]][Mom[S[2]] - Mom[-S[4]]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g5n45_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV5cV6, const std::complex<double>&
      CpS2cS4V3, const std::complex<double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-2)*CpS1cV5cV6*
      CpS2cS4V3*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T10G6N46 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 * @param[in] CpS2V3cV4 coupling Cp[S[2], V[3], -V[4]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g6n46_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF5F6V4PR, const std::complex<double>&
      CpF5F6V4PL, const std::complex<double>& CpcF6cF5S1PL, const std::complex<
      double>& CpcF6cF5S1PR, const std::complex<double>& CpS2V3cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*CpS2V3cV4*(
      b0tmp1*(CpcF6cF5S1PR*CpF5F6V4PL + CpcF6cF5S1PL*CpF5F6V4PR)*mLF5 + b1tmp2*
      (CpcF6cF5S1PR*CpF5F6V4PL*mLF5 + CpcF6cF5S1PL*CpF5F6V4PR*mLF5 +
      CpcF6cF5S1PL*CpF5F6V4PL*mLF6 + CpcF6cF5S1PR*CpF5F6V4PR*mLF6))*Den(Sqr(
      mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T10G7N47 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS2V3cV4 coupling Cp[S[2], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g7n47_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpS1cS5cS6, const std::complex<double>&
      CpS2V3cV4, const std::complex<double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-0.5)*(b0tmp1 +
      2*b1tmp2)*CpS1cS5cS6*CpS2V3cV4*CpS5S6V4*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T10G8N48 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpS2V3cV4 coupling Cp[S[2], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpU5U6V4 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g8n48_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpS1cU6cU5, const std::complex<double>&
      CpS2V3cV4, const std::complex<double>& CpU5U6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b1tmp1 = lib.B1(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,0.5)*b1tmp1*
      CpS1cU6cU5*CpS2V3cV4*CpU5U6V4*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T10G9N49 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS2V3cV4 coupling Cp[S[2], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g9n49_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpS1cS5cV6, const std::complex<double>&
      CpS2V3cV4, const std::complex<double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,-1)*(b0tmp1 -
      b1tmp2)*CpS1cS5cV6*CpS2V3cV4*CpS5V4V6*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T10G10N50 diagram for process S -> SV
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field S[2]
 * @param[in] mext3 mass of external field V[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS2V3cV4 coupling Cp[S[2], V[3], -V[4]][g[lt2, lt3]]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SSV calculate_diagram_SSV_t10g10n50_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpS1cV5cV6, const std::complex<double>&
      CpS2V3cV4, const std::complex<double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SSV result;

   result.m_decay = mext1;
   result.m_scalar = mext2;
   result.m_vector = mext3;

   result.form_factor = oneOver16PiSqr*(std::complex<double>(0,1.5)*(b0tmp1 + 2
      *b1tmp2)*CpS1cV5cV6*CpS2V3cV4*CpV4V5V6*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T1G1N1 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mLF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpF2F4cS6PL coupling Cp[F[2], F[4], -S[6]][PL]
 * @param[in] CpF2F4cS6PR coupling Cp[F[2], F[4], -S[6]][PR]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F3S6PL coupling Cp[F[5], F[3], S[6]][PL]
 * @param[in] CpF5F3S6PR coupling Cp[F[5], F[3], S[6]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g1n1_FFS(
   double mext1, double mext2, double mext3,
   double mLF4, double mLF5, double mLS6,
   const std::complex<double>& CpF2F4cS6PL, const std::complex<double>&
      CpF2F4cS6PR, const std::complex<double>& CpcF4cF5S1PL, const std::complex
      <double>& CpcF4cF5S1PR, const std::complex<double>& CpF5F3S6PL, const std
      ::complex<double>& CpF5F3S6PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLS6*mLS6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext1*mext1, mext3*mext3, mext2*
      mext2, mLF4*mLF4, mLF5*mLF5, mLS6*mLS6, scale*scale);
   const auto c1tmp3 = lib.C1(mext1*mext1, mext3*mext3, mext2*
      mext2, mLF4*mLF4, mLF5*mLF5, mLS6*mLS6, scale*scale);
   const auto c2tmp4 = lib.C2(mext1*mext1, mext3*mext3, mext2*
      mext2, mLF4*mLF4, mLF5*mLF5, mLS6*mLS6, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1*
      CpcF4cF5S1PR*CpF2F4cS6PL*CpF5F3S6PL - c2tmp4*CpcF4cF5S1PL*CpF2F4cS6PR*
      CpF5F3S6PR*mext2*mext3 - c2tmp4*CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PL*mext2
      *mLF4 - c2tmp4*CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PL*mext2*mLF5 - c1tmp3*
      CpF2F4cS6PR*CpF5F3S6PL*mext2*(CpcF4cF5S1PR*mLF4 + CpcF4cF5S1PL*mLF5) +
      c0tmp2*mLF4*(-(CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PL*mext2) + CpcF4cF5S1PL*
      CpF2F4cS6PL*CpF5F3S6PR*mext3 + CpcF4cF5S1PR*CpF2F4cS6PL*CpF5F3S6PL*mLF4 +
      CpcF4cF5S1PL*CpF2F4cS6PL*CpF5F3S6PL*mLF5) + c1tmp3*CpF2F4cS6PL*(
      CpcF4cF5S1PL*CpF5F3S6PR*mext3*mLF4 + CpcF4cF5S1PR*CpF5F3S6PR*mext3*mLF5 +
      CpcF4cF5S1PR*CpF5F3S6PL*Sqr(mext1)) + c2tmp4*CpcF4cF5S1PR*CpF2F4cS6PL*
      CpF5F3S6PL*Sqr(mext2)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1
      *CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PR - c0tmp2*CpcF4cF5S1PL*CpF2F4cS6PL*
      CpF5F3S6PR*mext2*mLF4 + c0tmp2*CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PL*mext3*
      mLF4 + c0tmp2*CpcF4cF5S1PR*CpF2F4cS6PR*CpF5F3S6PR*mLF4*mLF5 + c1tmp3*
      CpF2F4cS6PR*CpF5F3S6PL*mext3*(CpcF4cF5S1PR*mLF4 + CpcF4cF5S1PL*mLF5) -
      c2tmp4*CpF2F4cS6PL*mext2*(CpcF4cF5S1PR*CpF5F3S6PL*mext3 + CpcF4cF5S1PL*
      CpF5F3S6PR*mLF4 + CpcF4cF5S1PR*CpF5F3S6PR*mLF5) + c1tmp3*CpF5F3S6PR*(-(
      CpF2F4cS6PL*mext2*(CpcF4cF5S1PL*mLF4 + CpcF4cF5S1PR*mLF5)) + CpcF4cF5S1PL
      *CpF2F4cS6PR*Sqr(mext1)) + c2tmp4*CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PR*Sqr
      (mext2) + c0tmp2*CpcF4cF5S1PL*CpF2F4cS6PR*CpF5F3S6PR*Sqr(mLF4)));

   return result;
}

/**
 * @brief Evaluates T1G2N2 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF2cF6S4PL coupling Cp[F[2], -F[6], S[4]][PL]
 * @param[in] CpF2cF6S4PR coupling Cp[F[2], -F[6], S[4]][PR]
 * @param[in] CpF6F3S5PL coupling Cp[F[6], F[3], S[5]][PL]
 * @param[in] CpF6F3S5PR coupling Cp[F[6], F[3], S[5]][PR]
 * @param[in] CpS1cS4cS5 coupling Cp[S[1], -S[4], -S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g2n2_SSF(
   double mext1, double mext2, double mext3,
   double mLS4, double mLS5, double mLF6,
   const std::complex<double>& CpF2cF6S4PL, const std::complex<double>&
      CpF2cF6S4PR, const std::complex<double>& CpF6F3S5PL, const std::complex<
      double>& CpF6F3S5PR, const std::complex<double>& CpS1cS4cS5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLS4*mLS4, mLS5*mLS5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLS4*mLS4, mLS5*mLS5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLS4*mLS4, mLS5*mLS5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cS5*(c1tmp2*CpF2cF6S4PR*CpF6F3S5PL*mext2 + c2tmp3*CpF2cF6S4PL*
      CpF6F3S5PR*mext3 - c0tmp1*CpF2cF6S4PL*CpF6F3S5PL*mLF6));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cS5*(c1tmp2*CpF2cF6S4PL*CpF6F3S5PR*mext2 + c2tmp3*CpF2cF6S4PR*
      CpF6F3S5PL*mext3 - c0tmp1*CpF2cF6S4PR*CpF6F3S5PR*mLF6));

   return result;
}

/**
 * @brief Evaluates T1G3N3 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mLF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpF2F4cV6PR coupling Cp[F[2], F[4], -V[6]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpF2F4cV6PL coupling Cp[F[2], F[4], -V[6]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpcF4cF5S1PL coupling Cp[-F[4], -F[5], S[1]][PL]
 * @param[in] CpcF4cF5S1PR coupling Cp[-F[4], -F[5], S[1]][PR]
 * @param[in] CpF5F3V6PL coupling Cp[F[5], F[3], V[6]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F3V6PR coupling Cp[F[5], F[3], V[6]][LorentzProduct[gamma[lt3
    ], PR]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g3n3_FFV(
   double mext1, double mext2, double mext3,
   double mLF4, double mLF5, double mLV6,
   const std::complex<double>& CpF2F4cV6PR, const std::complex<double>&
      CpF2F4cV6PL, const std::complex<double>& CpcF4cF5S1PL, const std::complex
      <double>& CpcF4cF5S1PR, const std::complex<double>& CpF5F3V6PL, const std
      ::complex<double>& CpF5F3V6PR,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLV6*mLV6,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext1*mext1, mext3*mext3, mext2*
      mext2, mLF4*mLF4, mLF5*mLF5, mLV6*mLV6, scale*scale);
   const auto c1tmp3 = lib.C1(mext1*mext1, mext3*mext3, mext2*
      mext2, mLF4*mLF4, mLF5*mLF5, mLV6*mLV6, scale*scale);
   const auto c2tmp4 = lib.C2(mext1*mext1, mext3*mext3, mext2*
      mext2, mLF4*mLF4, mLF5*mLF5, mLV6*mLV6, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-2)*(2*
      b0tmp1*CpcF4cF5S1PL*CpF2F4cV6PR*CpF5F3V6PL + CpcF4cF5S1PR*(-(c1tmp3*
      CpF2F4cV6PR*CpF5F3V6PR*mext3*mLF4) + c1tmp3*CpF2F4cV6PL*CpF5F3V6PL*mext2*
      mLF5 + c2tmp4*CpF2F4cV6PL*CpF5F3V6PL*mext2*mLF5 + c0tmp2*CpF2F4cV6PR*mLF4
      *(-(CpF5F3V6PR*mext3) + 2*CpF5F3V6PL*mLF5)) + CpcF4cF5S1PL*((c0tmp2 +
      c1tmp3 + c2tmp4)*CpF2F4cV6PL*CpF5F3V6PL*mext2*mLF4 - c1tmp3*CpF2F4cV6PR*
      CpF5F3V6PR*mext3*mLF5 + CpF2F4cV6PR*CpF5F3V6PL*(-finite + 2*c1tmp3*Sqr(
      mext1) + c2tmp4*Sqr(mext1) + c2tmp4*Sqr(mext2) - c2tmp4*Sqr(mext3) + 2*
      c0tmp2*Sqr(mLF4)))));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,2)*(-2*
      b0tmp1*CpcF4cF5S1PR*CpF2F4cV6PL*CpF5F3V6PR + CpcF4cF5S1PL*(c1tmp3*
      CpF2F4cV6PL*CpF5F3V6PL*mext3*mLF4 - c1tmp3*CpF2F4cV6PR*CpF5F3V6PR*mext2*
      mLF5 - c2tmp4*CpF2F4cV6PR*CpF5F3V6PR*mext2*mLF5 + c0tmp2*CpF2F4cV6PL*mLF4
      *(CpF5F3V6PL*mext3 - 2*CpF5F3V6PR*mLF5)) + CpcF4cF5S1PR*(-((c0tmp2 +
      c1tmp3 + c2tmp4)*CpF2F4cV6PR*CpF5F3V6PR*mext2*mLF4) + c1tmp3*CpF2F4cV6PL*
      CpF5F3V6PL*mext3*mLF5 + CpF2F4cV6PL*CpF5F3V6PR*(finite - 2*c1tmp3*Sqr(
      mext1) - c2tmp4*Sqr(mext1) - c2tmp4*Sqr(mext2) + c2tmp4*Sqr(mext3) - 2*
      c0tmp2*Sqr(mLF4)))));

   return result;
}

/**
 * @brief Evaluates T1G4N4 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mLS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF2cF6S4PL coupling Cp[F[2], -F[6], S[4]][PL]
 * @param[in] CpF2cF6S4PR coupling Cp[F[2], -F[6], S[4]][PR]
 * @param[in] CpF6F3V5PL coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF6F3V5PR coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpS1cS4cV5 coupling Cp[S[1], -S[4], -V[5]][Mom[S[1]] - Mom[-S[4]]
    ]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g4n4_SVF(
   double mext1, double mext2, double mext3,
   double mLS4, double mLV5, double mLF6,
   const std::complex<double>& CpF2cF6S4PL, const std::complex<double>&
      CpF2cF6S4PR, const std::complex<double>& CpF6F3V5PL, const std::complex<
      double>& CpF6F3V5PR, const std::complex<double>& CpS1cS4cV5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF6*mLF6, mLV5*mLV5,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLS4*mLS4, mLV5*mLV5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLS4*mLS4, mLV5*mLV5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLS4*mLS4, mLV5*mLV5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpS1cS4cV5*(b0tmp1*CpF2cF6S4PL*CpF6F3V5PL - c1tmp3*CpF2cF6S4PR*CpF6F3V5PR
      *mext2*mext3 - 2*c0tmp2*CpF2cF6S4PR*CpF6F3V5PL*mext2*mLF6 - c1tmp3*
      CpF2cF6S4PR*CpF6F3V5PL*mext2*mLF6 + c0tmp2*CpF2cF6S4PL*CpF6F3V5PR*mext3*
      mLF6 - c2tmp4*(CpF6F3V5PR*mext3*(2*CpF2cF6S4PR*mext2 + CpF2cF6S4PL*mLF6)
      + CpF2cF6S4PL*CpF6F3V5PL*(Sqr(mext1) - Sqr(mext2))) - c0tmp2*CpF2cF6S4PL*
      CpF6F3V5PL*Sqr(mext2) + c0tmp2*CpF2cF6S4PL*CpF6F3V5PL*Sqr(mLS4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS4cV5*(-(b0tmp1*CpF2cF6S4PR*CpF6F3V5PR) + 2*c2tmp4*CpF2cF6S4PL*
      CpF6F3V5PL*mext2*mext3 + 2*c0tmp2*CpF2cF6S4PL*CpF6F3V5PR*mext2*mLF6 -
      c0tmp2*CpF2cF6S4PR*CpF6F3V5PL*mext3*mLF6 + c1tmp3*CpF2cF6S4PL*mext2*(
      CpF6F3V5PL*mext3 + CpF6F3V5PR*mLF6) + c2tmp4*CpF2cF6S4PR*(CpF6F3V5PL*
      mext3*mLF6 + CpF6F3V5PR*(Sqr(mext1) - Sqr(mext2))) + c0tmp2*CpF2cF6S4PR*
      CpF6F3V5PR*Sqr(mext2) - c0tmp2*CpF2cF6S4PR*CpF6F3V5PR*Sqr(mLS4)));

   return result;
}

/**
 * @brief Evaluates T1G5N5 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF2cF6V4PR coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpF2cF6V4PL coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF6F3S5PL coupling Cp[F[6], F[3], S[5]][PL]
 * @param[in] CpF6F3S5PR coupling Cp[F[6], F[3], S[5]][PR]
 * @param[in] CpS1cS5cV4 coupling Cp[S[1], -S[5], -V[4]][Mom[S[1]] - Mom[-S[5]]
    ]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g5n5_VSF(
   double mext1, double mext2, double mext3,
   double mLV4, double mLS5, double mLF6,
   const std::complex<double>& CpF2cF6V4PR, const std::complex<double>&
      CpF2cF6V4PL, const std::complex<double>& CpF6F3S5PL, const std::complex<
      double>& CpF6F3S5PR, const std::complex<double>& CpS1cS5cV4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF6*mLF6, mLS5*mLS5,
      scale*scale);
   const auto c0tmp2 = lib.C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLV4*mLV4, mLS5*mLS5, scale*scale);
   const auto c1tmp3 = lib.C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLV4*mLV4, mLS5*mLS5, scale*scale);
   const auto c2tmp4 = lib.C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLV4*mLV4, mLS5*mLS5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpS1cS5cV4*(b0tmp1*CpF2cF6V4PR*CpF6F3S5PL - c2tmp4*CpF2cF6V4PL*CpF6F3S5PR
      *mext2*mext3 - c2tmp4*CpF2cF6V4PR*CpF6F3S5PR*mext3*mLF6 + c2tmp4*
      CpF2cF6V4PR*CpF6F3S5PL*Sqr(mext1) - c2tmp4*CpF2cF6V4PR*CpF6F3S5PL*Sqr(
      mext2) - c1tmp3*(CpF2cF6V4PL*mext2*(2*CpF6F3S5PR*mext3 + CpF6F3S5PL*mLF6)
      + CpF2cF6V4PR*CpF6F3S5PL*(2*Sqr(mext1) + Sqr(mext2) - 2*Sqr(mext3))) +
      c2tmp4*CpF2cF6V4PR*CpF6F3S5PL*Sqr(mext3) + c0tmp2*(CpF2cF6V4PL*CpF6F3S5PL
      *mext2*mLF6 - 2*CpF2cF6V4PR*CpF6F3S5PR*mext3*mLF6 + CpF2cF6V4PR*
      CpF6F3S5PL*(-Sqr(mext2) + Sqr(mLV4)))));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpS1cS5cV4*(-(b0tmp1*CpF2cF6V4PL*CpF6F3S5PR) + 2*c1tmp3*CpF2cF6V4PR*
      CpF6F3S5PL*mext2*mext3 + c1tmp3*CpF2cF6V4PR*CpF6F3S5PR*mext2*mLF6 + 2*
      c0tmp2*CpF2cF6V4PL*CpF6F3S5PL*mext3*mLF6 + c2tmp4*CpF6F3S5PL*mext3*(
      CpF2cF6V4PR*mext2 + CpF2cF6V4PL*mLF6) + 2*c1tmp3*CpF2cF6V4PL*CpF6F3S5PR*
      Sqr(mext1) + c1tmp3*CpF2cF6V4PL*CpF6F3S5PR*Sqr(mext2) - 2*c1tmp3*
      CpF2cF6V4PL*CpF6F3S5PR*Sqr(mext3) - c2tmp4*CpF2cF6V4PL*CpF6F3S5PR*(Sqr(
      mext1) - Sqr(mext2) + Sqr(mext3)) - c0tmp2*CpF6F3S5PR*(CpF2cF6V4PR*mext2*
      mLF6 + CpF2cF6V4PL*(-Sqr(mext2) + Sqr(mLV4)))));

   return result;
}

/**
 * @brief Evaluates T1G6N6 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mLV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF2cF6V4PL coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2cF6V4PR coupling Cp[F[2], -F[6], V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpF6F3V5PL coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF6F3V5PR coupling Cp[F[6], F[3], V[5]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpS1cV4cV5 coupling Cp[S[1], -V[4], -V[5]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t1g6n6_VVF(
   double mext1, double mext2, double mext3,
   double mLV4, double mLV5, double mLF6,
   const std::complex<double>& CpF2cF6V4PL, const std::complex<double>&
      CpF2cF6V4PR, const std::complex<double>& CpF6F3V5PL, const std::complex<
      double>& CpF6F3V5PR, const std::complex<double>& CpS1cV4cV5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto c0tmp1 = lib.C0(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLV4*mLV4, mLV5*mLV5, scale*scale);
   const auto c1tmp2 = lib.C1(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLV4*mLV4, mLV5*mLV5, scale*scale);
   const auto c2tmp3 = lib.C2(mext2*mext2, mext1*mext1, mext3*
      mext3, mLF6*mLF6, mLV4*mLV4, mLV5*mLV5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,2)*
      CpS1cV4cV5*(c1tmp2*CpF2cF6V4PL*CpF6F3V5PL*mext2 + c2tmp3*CpF2cF6V4PR*
      CpF6F3V5PR*mext3 + 2*c0tmp1*CpF2cF6V4PR*CpF6F3V5PL*mLF6));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,2)*
      CpS1cV4cV5*(c1tmp2*CpF2cF6V4PR*CpF6F3V5PR*mext2 + c2tmp3*CpF2cF6V4PL*
      CpF6F3V5PL*mext3 + 2*c0tmp1*CpF2cF6V4PL*CpF6F3V5PR*mLF6));

   return result;
}

/**
 * @brief Evaluates T2G1N7 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpS1S4cS5S5 coupling Cp[S[1], S[4], -S[5], S[5]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t2g1n7_SS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpS1S4cS5S5,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLS5*mLS5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(0.5*a0tmp1*CpF2F3cS4PL*CpS1S4cS5S5
      *Den(Sqr(mext1),Sqr(mIS4)));
   result.form_factor_right = oneOver16PiSqr*(0.5*a0tmp1*CpF2F3cS4PR*
      CpS1S4cS5S5*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T2G2N8 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpS1S4cV5V5 coupling Cp[S[1], S[4], -V[5], V[5]][g[lt3, lt4]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t2g2n8_SV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpS1S4cV5V5,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV5*mLV5, scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(-(CpF2F3cS4PL*CpS1S4cV5V5*Den(Sqr(
      mext1),Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5))));
   result.form_factor_right = oneOver16PiSqr*(-(CpF2F3cS4PR*CpS1S4cV5V5*Den(Sqr
      (mext1),Sqr(mIS4))*(2*a0tmp1 - finite*Sqr(mLV5))));

   return result;
}

/**
 * @brief Evaluates T3G1N9 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpF2cF4S1PL coupling Cp[F[2], -F[4], S[1]][PL]
 * @param[in] CpF2cF4S1PR coupling Cp[F[2], -F[4], S[1]][PR]
 * @param[in] CpF4F5S6PL coupling Cp[F[4], F[5], S[6]][PL]
 * @param[in] CpF4F5S6PR coupling Cp[F[4], F[5], S[6]][PR]
 * @param[in] CpcF5F3cS6PL coupling Cp[-F[5], F[3], -S[6]][PL]
 * @param[in] CpcF5F3cS6PR coupling Cp[-F[5], F[3], -S[6]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t3g1n9_FFS(
   double mext1, double mext2, double mext3,
   double mIF4, double mLF5, double mLS6,
   const std::complex<double>& CpF2cF4S1PL, const std::complex<double>&
      CpF2cF4S1PR, const std::complex<double>& CpF4F5S6PL, const std::complex<
      double>& CpF4F5S6PR, const std::complex<double>& CpcF5F3cS6PL, const std
      ::complex<double>& CpcF5F3cS6PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLF5*mLF5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpF2cF4S1PL*(-(b1tmp2*mext3*(CpcF5F3cS6PL*CpF4F5S6PR*mext3 + CpcF5F3cS6PR
      *CpF4F5S6PL*mIF4)) + b0tmp1*(CpcF5F3cS6PR*CpF4F5S6PR*mext3 + CpcF5F3cS6PL
      *CpF4F5S6PL*mIF4)*mLF5)*Den(Sqr(mext3),Sqr(mIF4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpF2cF4S1PR*(-(b1tmp2*mext3*(CpcF5F3cS6PR*CpF4F5S6PL*mext3 + CpcF5F3cS6PL
      *CpF4F5S6PR*mIF4)) + b0tmp1*(CpcF5F3cS6PL*CpF4F5S6PL*mext3 + CpcF5F3cS6PR
      *CpF4F5S6PR*mIF4)*mLF5)*Den(Sqr(mext3),Sqr(mIF4)));

   return result;
}

/**
 * @brief Evaluates T3G2N10 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpF2cF4S1PL coupling Cp[F[2], -F[4], S[1]][PL]
 * @param[in] CpF2cF4S1PR coupling Cp[F[2], -F[4], S[1]][PR]
 * @param[in] CpF4F5V6PL coupling Cp[F[4], F[5], V[6]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF4F5V6PR coupling Cp[F[4], F[5], V[6]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpcF5F3cV6PL coupling Cp[-F[5], F[3], -V[6]][LorentzProduct[gamma
    [lt3], PL]]
 * @param[in] CpcF5F3cV6PR coupling Cp[-F[5], F[3], -V[6]][LorentzProduct[gamma
    [lt3], PR]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t3g2n10_FFV(
   double mext1, double mext2, double mext3,
   double mIF4, double mLF5, double mLV6,
   const std::complex<double>& CpF2cF4S1PL, const std::complex<double>&
      CpF2cF4S1PR, const std::complex<double>& CpF4F5V6PL, const std::complex<
      double>& CpF4F5V6PR, const std::complex<double>& CpcF5F3cV6PL, const std
      ::complex<double>& CpcF5F3cV6PR,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext3*mext3, mLF5*mLF5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext3*mext3, mLF5*mLF5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpF2cF4S1PL*Den(Sqr(mext3),Sqr(mIF4))*(-2*CpcF5F3cV6PR*CpF4F5V6PL*(-2*
      b0tmp1 + finite)*mext3*mLF5 + CpF4F5V6PR*mIF4*(CpcF5F3cV6PR*(2*b1tmp2 +
      finite)*mext3 - 2*CpcF5F3cV6PL*(-2*b0tmp1 + finite)*mLF5) + CpcF5F3cV6PL*
      CpF4F5V6PL*(2*b1tmp2 + finite)*Sqr(mext3)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpF2cF4S1PR*Den(Sqr(mext3),Sqr(mIF4))*(CpcF5F3cV6PL*CpF4F5V6PL*(2*b1tmp2
      + finite)*mext3*mIF4 - 2*CpcF5F3cV6PL*CpF4F5V6PR*(-2*b0tmp1 + finite)*
      mext3*mLF5 - 2*CpcF5F3cV6PR*CpF4F5V6PL*(-2*b0tmp1 + finite)*mIF4*mLF5 +
      CpcF5F3cV6PR*CpF4F5V6PR*(2*b1tmp2 + finite)*Sqr(mext3)));

   return result;
}

/**
 * @brief Evaluates T4G1N11 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpcF4F3S1PL coupling Cp[-F[4], F[3], S[1]][PL]
 * @param[in] CpcF4F3S1PR coupling Cp[-F[4], F[3], S[1]][PR]
 * @param[in] CpF2cF5cS6PL coupling Cp[F[2], -F[5], -S[6]][PL]
 * @param[in] CpF2cF5cS6PR coupling Cp[F[2], -F[5], -S[6]][PR]
 * @param[in] CpF5F4S6PL coupling Cp[F[5], F[4], S[6]][PL]
 * @param[in] CpF5F4S6PR coupling Cp[F[5], F[4], S[6]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t4g1n11_FFS(
   double mext1, double mext2, double mext3,
   double mIF4, double mLF5, double mLS6,
   const std::complex<double>& CpcF4F3S1PL, const std::complex<double>&
      CpcF4F3S1PR, const std::complex<double>& CpF2cF5cS6PL, const std::complex
      <double>& CpF2cF5cS6PR, const std::complex<double>& CpF5F4S6PL, const std
      ::complex<double>& CpF5F4S6PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLF5*mLF5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLF5*mLF5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpcF4F3S1PL*(-(b1tmp2*mext2*(CpF2cF5cS6PL*CpF5F4S6PR*mext2 + CpF2cF5cS6PR
      *CpF5F4S6PL*mIF4)) + b0tmp1*(CpF2cF5cS6PR*CpF5F4S6PR*mext2 + CpF2cF5cS6PL
      *CpF5F4S6PL*mIF4)*mLF5)*Den(Sqr(mext2),Sqr(mIF4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpcF4F3S1PR*(-(b1tmp2*mext2*(CpF2cF5cS6PR*CpF5F4S6PL*mext2 + CpF2cF5cS6PL
      *CpF5F4S6PR*mIF4)) + b0tmp1*(CpF2cF5cS6PL*CpF5F4S6PL*mext2 + CpF2cF5cS6PR
      *CpF5F4S6PR*mIF4)*mLF5)*Den(Sqr(mext2),Sqr(mIF4)));

   return result;
}

/**
 * @brief Evaluates T4G2N12 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIF4 mass of internal field F[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpcF4F3S1PL coupling Cp[-F[4], F[3], S[1]][PL]
 * @param[in] CpcF4F3S1PR coupling Cp[-F[4], F[3], S[1]][PR]
 * @param[in] CpF2cF5cV6PL coupling Cp[F[2], -F[5], -V[6]][LorentzProduct[gamma
    [lt3], PL]]
 * @param[in] CpF2cF5cV6PR coupling Cp[F[2], -F[5], -V[6]][LorentzProduct[gamma
    [lt3], PR]]
 * @param[in] CpF5F4V6PL coupling Cp[F[5], F[4], V[6]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpF5F4V6PR coupling Cp[F[5], F[4], V[6]][LorentzProduct[gamma[lt3
    ], PR]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t4g2n12_FFV(
   double mext1, double mext2, double mext3,
   double mIF4, double mLF5, double mLV6,
   const std::complex<double>& CpcF4F3S1PL, const std::complex<double>&
      CpcF4F3S1PR, const std::complex<double>& CpF2cF5cV6PL, const std::complex
      <double>& CpF2cF5cV6PR, const std::complex<double>& CpF5F4V6PL, const std
      ::complex<double>& CpF5F4V6PR,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext2*mext2, mLF5*mLF5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext2*mext2, mLF5*mLF5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpcF4F3S1PL*Den(Sqr(mext2),Sqr(mIF4))*(CpF2cF5cV6PL*CpF5F4V6PL*(2*b1tmp2
      + finite)*mext2*mIF4 - 2*CpF2cF5cV6PL*CpF5F4V6PR*(-2*b0tmp1 + finite)*
      mext2*mLF5 - 2*CpF2cF5cV6PR*CpF5F4V6PL*(-2*b0tmp1 + finite)*mIF4*mLF5 +
      CpF2cF5cV6PR*CpF5F4V6PR*(2*b1tmp2 + finite)*Sqr(mext2)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpcF4F3S1PR*Den(Sqr(mext2),Sqr(mIF4))*(-2*CpF2cF5cV6PL*CpF5F4V6PR*(-2*
      b0tmp1 + finite)*mIF4*mLF5 + CpF2cF5cV6PR*mext2*(CpF5F4V6PR*(2*b1tmp2 +
      finite)*mIF4 - 2*CpF5F4V6PL*(-2*b0tmp1 + finite)*mLF5) + CpF2cF5cV6PL*
      CpF5F4V6PL*(2*b1tmp2 + finite)*Sqr(mext2)));

   return result;
}

/**
 * @brief Evaluates T5G1N13 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpF5F6S4PL coupling Cp[F[5], F[6], S[4]][PL]
 * @param[in] CpF5F6S4PR coupling Cp[F[5], F[6], S[4]][PR]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g1n13_SFF(
   double mext1, double mext2, double mext3,
   double mIS4, double mLF5, double mLF6,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpF5F6S4PL, const std::complex<
      double>& CpF5F6S4PR, const std::complex<double>& CpcF6cF5S1PL, const std
      ::complex<double>& CpcF6cF5S1PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLF6*mLF6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpF2F3cS4PL*Den(Sqr(mext1),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S1PR*
      CpF5F6S4PL*mLF5 + CpcF6cF5S1PL*CpF5F6S4PR*mLF5 + CpcF6cF5S1PL*CpF5F6S4PL*
      mLF6 + CpcF6cF5S1PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S1PR*CpF5F6S4PL +
      CpcF6cF5S1PL*CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext1))));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpF2F3cS4PR*Den(Sqr(mext1),Sqr(mIS4))*(b0tmp2*mLF5*(CpcF6cF5S1PR*
      CpF5F6S4PL*mLF5 + CpcF6cF5S1PL*CpF5F6S4PR*mLF5 + CpcF6cF5S1PL*CpF5F6S4PL*
      mLF6 + CpcF6cF5S1PR*CpF5F6S4PR*mLF6) + (CpcF6cF5S1PR*CpF5F6S4PL +
      CpcF6cF5S1PL*CpF5F6S4PR)*(a0tmp1 + b1tmp3*Sqr(mext1))));

   return result;
}

/**
 * @brief Evaluates T5G2N14 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS4S5S6 coupling Cp[S[4], S[5], S[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g2n14_SSS(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLS6,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpS1cS5cS6, const std::complex<
      double>& CpS4S5S6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,0.5)*b0tmp1
      *CpF2F3cS4PL*CpS1cS5cS6*CpS4S5S6*Den(Sqr(mext1),Sqr(mIS4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,0.5)*
      b0tmp1*CpF2F3cS4PR*CpS1cS5cS6*CpS4S5S6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T5G3N15 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpS4U5U6 coupling Cp[S[4], U[5], U[6]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g3n15_SUU(
   double mext1, double mext2, double mext3,
   double mIS4, double mLU5, double mLU6,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpS1cU6cU5, const std::complex<
      double>& CpS4U5U6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-0.5)*
      b0tmp1*CpF2F3cS4PL*CpS1cU6cU5*CpS4U5U6*Den(Sqr(mext1),Sqr(mIS4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-0.5)*
      b0tmp1*CpF2F3cS4PR*CpS1cU6cU5*CpS4U5U6*Den(Sqr(mext1),Sqr(mIS4)));

   return result;
}

/**
 * @brief Evaluates T5G4N16 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS4S5V6 coupling Cp[S[4], S[5], V[6]][Mom[S[4]] - Mom[S[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g4n16_SSV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLS5, double mLV6,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpS1cS5cV6, const std::complex<
      double>& CpS4S5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto a0tmp1 = lib.A0(mLV6*mLV6, scale*scale);
   const auto b0tmp2 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp3 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpF2F3cS4PL*CpS1cS5cV6*CpS4S5V6*Den(Sqr(mext1),Sqr(mIS4))*(a0tmp1 - 2*
      b1tmp3*Sqr(mext1) + b0tmp2*(Sqr(mext1) + Sqr(mLS5))));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*
      CpF2F3cS4PR*CpS1cS5cV6*CpS4S5V6*Den(Sqr(mext1),Sqr(mIS4))*(a0tmp1 - 2*
      b1tmp3*Sqr(mext1) + b0tmp2*(Sqr(mext1) + Sqr(mLS5))));

   return result;
}

/**
 * @brief Evaluates T5G5N17 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIS4 mass of internal field S[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpF2F3cS4PL coupling Cp[F[2], F[3], -S[4]][PL]
 * @param[in] CpF2F3cS4PR coupling Cp[F[2], F[3], -S[4]][PR]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpS4V5V6 coupling Cp[S[4], V[5], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g5n17_SVV(
   double mext1, double mext2, double mext3,
   double mIS4, double mLV5, double mLV6,
   const std::complex<double>& CpF2F3cS4PL, const std::complex<double>&
      CpF2F3cS4PR, const std::complex<double>& CpS1cV5cV6, const std::complex<
      double>& CpS4V5V6,
   double scale, double finite)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpF2F3cS4PL*CpS1cV5cV6*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext1),Sqr(
      mIS4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,-1)*
      CpF2F3cS4PR*CpS1cV5cV6*CpS4V5V6*(-2*b0tmp1 + finite)*Den(Sqr(mext1),Sqr(
      mIS4)));

   return result;
}

/**
 * @brief Evaluates T5G6N18 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLF5 mass of internal field F[5]
 * @param[in] mLF6 mass of internal field F[6]
 * @param[in] CpF2F3cV4PL coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2F3cV4PR coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpF5F6V4PR coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PR]]
 * @param[in] CpF5F6V4PL coupling Cp[F[5], F[6], V[4]][LorentzProduct[gamma[lt3
    ], PL]]
 * @param[in] CpcF6cF5S1PL coupling Cp[-F[6], -F[5], S[1]][PL]
 * @param[in] CpcF6cF5S1PR coupling Cp[-F[6], -F[5], S[1]][PR]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g6n18_VFF(
   double mext1, double mext2, double mext3,
   double mIV4, double mLF5, double mLF6,
   const std::complex<double>& CpF2F3cV4PL, const std::complex<double>&
      CpF2F3cV4PR, const std::complex<double>& CpF5F6V4PR, const std::complex<
      double>& CpF5F6V4PL, const std::complex<double>& CpcF6cF5S1PL, const std
      ::complex<double>& CpcF6cF5S1PR,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLF5*mLF5, mLF6*mLF6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*(
      CpF2F3cV4PL*mext2 - CpF2F3cV4PR*mext3)*(b0tmp1*(CpcF6cF5S1PR*CpF5F6V4PL +
      CpcF6cF5S1PL*CpF5F6V4PR)*mLF5 + b1tmp2*(CpcF6cF5S1PR*CpF5F6V4PL*mLF5 +
      CpcF6cF5S1PL*CpF5F6V4PR*mLF5 + CpcF6cF5S1PL*CpF5F6V4PL*mLF6 +
      CpcF6cF5S1PR*CpF5F6V4PR*mLF6))*Den(Sqr(mext1),Sqr(mIV4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*(
      CpF2F3cV4PR*mext2 - CpF2F3cV4PL*mext3)*(b0tmp1*(CpcF6cF5S1PR*CpF5F6V4PL +
      CpcF6cF5S1PL*CpF5F6V4PR)*mLF5 + b1tmp2*(CpcF6cF5S1PR*CpF5F6V4PL*mLF5 +
      CpcF6cF5S1PL*CpF5F6V4PR*mLF5 + CpcF6cF5S1PL*CpF5F6V4PL*mLF6 +
      CpcF6cF5S1PR*CpF5F6V4PR*mLF6))*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T5G7N19 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLS6 mass of internal field S[6]
 * @param[in] CpF2F3cV4PL coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2F3cV4PR coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpS1cS5cS6 coupling Cp[S[1], -S[5], -S[6]]
 * @param[in] CpS5S6V4 coupling Cp[S[5], S[6], V[4]][Mom[S[5]] - Mom[S[6]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g7n19_VSS(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLS6,
   const std::complex<double>& CpF2F3cV4PL, const std::complex<double>&
      CpF2F3cV4PR, const std::complex<double>& CpS1cS5cS6, const std::complex<
      double>& CpS5S6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLS6*mLS6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,0.5)*(
      b0tmp1 + 2*b1tmp2)*CpS1cS5cS6*CpS5S6V4*(CpF2F3cV4PL*mext2 - CpF2F3cV4PR*
      mext3)*Den(Sqr(mext1),Sqr(mIV4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,0.5)*(
      b0tmp1 + 2*b1tmp2)*CpS1cS5cS6*CpS5S6V4*(CpF2F3cV4PR*mext2 - CpF2F3cV4PL*
      mext3)*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T5G8N20 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLU5 mass of internal field U[5]
 * @param[in] mLU6 mass of internal field U[6]
 * @param[in] CpF2F3cV4PL coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2F3cV4PR coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpS1cU6cU5 coupling Cp[S[1], -U[6], -U[5]]
 * @param[in] CpU5U6V4 coupling Cp[U[5], U[6], V[4]][Mom[U[5]]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g8n20_VUU(
   double mext1, double mext2, double mext3,
   double mIV4, double mLU5, double mLU6,
   const std::complex<double>& CpF2F3cV4PL, const std::complex<double>&
      CpF2F3cV4PR, const std::complex<double>& CpS1cU6cU5, const std::complex<
      double>& CpU5U6V4,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b1tmp1 = lib.B1(mext1*mext1, mLU5*mLU5, mLU6*mLU6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-0.5)*
      b1tmp1*CpS1cU6cU5*CpU5U6V4*(CpF2F3cV4PL*mext2 - CpF2F3cV4PR*mext3)*Den(
      Sqr(mext1),Sqr(mIV4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,0.5)*
      b1tmp1*CpS1cU6cU5*CpU5U6V4*(-(CpF2F3cV4PR*mext2) + CpF2F3cV4PL*mext3)*Den
      (Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T5G9N21 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLS5 mass of internal field S[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpF2F3cV4PL coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2F3cV4PR coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpS1cS5cV6 coupling Cp[S[1], -S[5], -V[6]][Mom[S[1]] - Mom[-S[5]]
    ]
 * @param[in] CpS5V4V6 coupling Cp[S[5], V[4], V[6]][g[lt2, lt3]]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g9n21_VSV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLS5, double mLV6,
   const std::complex<double>& CpF2F3cV4PL, const std::complex<double>&
      CpF2F3cV4PR, const std::complex<double>& CpS1cS5cV6, const std::complex<
      double>& CpS5V4V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLS5*mLS5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1
      - b1tmp2)*CpS1cS5cV6*CpS5V4V6*(CpF2F3cV4PL*mext2 - CpF2F3cV4PR*mext3)*Den
      (Sqr(mext1),Sqr(mIV4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1)*(b0tmp1
       - b1tmp2)*CpS1cS5cV6*CpS5V4V6*(CpF2F3cV4PR*mext2 - CpF2F3cV4PL*mext3)*
      Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}

/**
 * @brief Evaluates T5G10N22 diagram for process S -> FF
 *
 * @param[in] mext1 mass of external field S[1]
 * @param[in] mext2 mass of external field F[2]
 * @param[in] mext3 mass of external field F[3]
 * @param[in] mIV4 mass of internal field V[4]
 * @param[in] mLV5 mass of internal field V[5]
 * @param[in] mLV6 mass of internal field V[6]
 * @param[in] CpF2F3cV4PL coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PL]]
 * @param[in] CpF2F3cV4PR coupling Cp[F[2], F[3], -V[4]][LorentzProduct[gamma[
    lt3], PR]]
 * @param[in] CpS1cV5cV6 coupling Cp[S[1], -V[5], -V[6]][g[lt2, lt3]]
 * @param[in] CpV4V5V6 coupling Cp[V[4], V[5], V[6]][g[lt1, lt2] (-Mom[V[4]] +
    Mom[V[5]]) + g[lt1, lt3] (Mom[V[4]] - Mom[V[6]]) + g[lt2, lt3] (-Mom[V[5]]
    + Mom[V[6]])]
 *
 * @return value of the one-loop diagram
 */
Decay_amplitude_SFF calculate_diagram_SFF_t5g10n22_VVV(
   double mext1, double mext2, double mext3,
   double mIV4, double mLV5, double mLV6,
   const std::complex<double>& CpF2F3cV4PL, const std::complex<double>&
      CpF2F3cV4PR, const std::complex<double>& CpS1cV5cV6, const std::complex<
      double>& CpV4V5V6,
   double scale)
{
   auto& lib = Loop_library::get();

   const auto b0tmp1 = lib.B0(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);
   const auto b1tmp2 = lib.B1(mext1*mext1, mLV5*mLV5, mLV6*mLV6,
      scale*scale);

   Decay_amplitude_SFF result;

   result.m_decay = mext1;
   result.m_fermion_1 = mext2;
   result.m_fermion_2 = mext3;

   result.form_factor_left = oneOver16PiSqr*(std::complex<double>(0,-1.5)*(
      b0tmp1 + 2*b1tmp2)*CpS1cV5cV6*CpV4V5V6*(CpF2F3cV4PL*mext2 - CpF2F3cV4PR*
      mext3)*Den(Sqr(mext1),Sqr(mIV4)));
   result.form_factor_right = oneOver16PiSqr*(std::complex<double>(0,1.5)*(
      b0tmp1 + 2*b1tmp2)*CpS1cV5cV6*CpV4V5V6*(-(CpF2F3cV4PR*mext2) +
      CpF2F3cV4PL*mext3)*Den(Sqr(mext1),Sqr(mIV4)));

   return result;
}


} // namespace flexiblesusy
