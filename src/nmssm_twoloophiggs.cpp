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

#include "nmssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.hpp"
#include "nmssm2loop.h"
#include "config.h"
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_nmssm)
#else
   #define LOCK_MUTEX()
#endif

using namespace flexiblesusy::mssm_twoloophiggs;

namespace flexiblesusy {
namespace nmssm_twoloophiggs {

static std::mutex mtx_nmssm; /// locks MSSM fortran functions

namespace {
template <typename T> T sqr(T a) { return a * a; }
}

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_at_as_nmssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq,
   double amu, double tanb, double vev2, double gs, double svevS)
{
   const double cosb = 1. / std::sqrt(1. + sqr(tanb));

   const Eigen::Matrix<double, 2, 1> t_mssm = tadpole_higgs_2loop_at_as_mssm(
      rmtsq, mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu, tanb, vev2, gs);

   Eigen::Matrix<double, 3, 1> result;
   result.head<2>() = t_mssm;
   // rescale T1 to get TS
   result(2) = t_mssm(0) * std::sqrt(vev2) * cosb / (svevS * std::sqrt(2.));

   return result;
}

Eigen::Matrix<double, 3, 1> tadpole_higgs_2loop_ab_as_nmssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq,
   double amu, double cotbeta, double vev2, double gs, double svevS)
{
   const double tanb = 1./cotbeta;
   const double sinb = tanb / std::sqrt(1. + sqr(tanb));

   const Eigen::Matrix<double, 2, 1> t_mssm = tadpole_higgs_2loop_ab_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu, cotbeta, vev2, gs);

   Eigen::Matrix<double, 3, 1> result;
   result.head<2>() = t_mssm;
   // rescale T1 to get TS
   result(2) = t_mssm(0) * std::sqrt(vev2) * sinb / (svevS * std::sqrt(2.));

   return result;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as, double amu)
{
   const double vev2 = sqr(vevS * std::sqrt(2.));
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
         rmt, mg, mst1sq, mst2sq, sxt, cxt, scalesq, tanb, vevS, lamS, svevS, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_at_as_nmssm(
         sqr(rmt), mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu, tanb, vev2, gs, svevS).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as, double amu)
{
   const double vev2 = sqr(vevS * std::sqrt(2.));
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
         rmb, mg, msb1sq, msb2sq, sxb, cxb, scalesq, cotb, vevS, lamS, svevS, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_ab_as_nmssm(
         sqr(rmb), mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu, cotb, vev2, gs, svevS).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as, double amu)
{
   const double vev2 = sqr(vevS * std::sqrt(2.));
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
         rmt, mg, mst1sq, mst2sq, sxt, cxt, scalesq, tanb, vevS,
         lamS, svevS, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_at_as_nmssm(
         sqr(rmt), mg, mst1sq, mst2sq, sxt, cxt, scalesq, amu, tanb, vev2, gs, svevS).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as, double amu)
{
   const double vev2 = sqr(vevS * std::sqrt(2.));
   const double gs = std::sqrt(as * 4. * M_PI);

   const Eigen::Matrix<double, 3, 3> se =
      self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
         rmb, mg, msb1sq, msb2sq, sxb, cxb, scalesq, cotb, vevS,
         lamS, svevS, as);

   const Eigen::Matrix<double, 3, 3> tadpoles =
      tadpole_higgs_2loop_ab_as_nmssm(
         sqr(rmb), mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu, cotb, vev2, gs, svevS).asDiagonal();

   return se + tadpoles;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as)
{
   int loop = 2;
   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};

   LOCK_MUTEX();

   effpot_(&loop, &rmt, &mg, &mst1sq, &mst2sq, &sxt, &cxt,
           &scalesq, &tanb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);

   Eigen::Matrix<double, 3, 3> result;
   result << DMS[0][0], DMS[0][1], DMS[0][2],
             DMS[1][0], DMS[1][1], DMS[1][2],
             DMS[2][0], DMS[2][1], DMS[2][2];

   return -result;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as)
{
   int loop = 2;
   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};

   LOCK_MUTEX();

   effpot_(&loop, &rmb, &mg, &msb1sq, &msb2sq, &sxb, &cxb,
           &scalesq, &cotb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);

   // Make appropriate substitutions for elements following 0907.4682
   // bottom of page 9
   std::swap(DMS[0][0], DMS[1][1]);
   std::swap(DMS[0][2], DMS[1][2]);

   Eigen::Matrix<double, 3, 3> result;
   result << DMS[0][0], DMS[0][1], DMS[0][2],
             DMS[1][0], DMS[1][1], DMS[1][2],
             DMS[2][0], DMS[2][1], DMS[2][2];

   return -result;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm_with_tadpoles(
   double rmt, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double tanb, double vevS,
   double lamS, double svevS, double as)
{
   int loop = 2;
   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};

   LOCK_MUTEX();

   effpot_(&loop, &rmt, &mg, &mst1sq, &mst2sq, &sxt, &cxt,
           &scalesq, &tanb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);

   Eigen::Matrix<double, 3, 3> result;
   result << DMP[0][0], DMP[0][1], DMP[0][2],
             DMP[1][0], DMP[1][1], DMP[1][2],
             DMP[2][0], DMP[2][1], DMP[2][2];

   return -result;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm_with_tadpoles(
   double rmb, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double cotb, double vevS,
   double lamS, double svevS, double as)
{
   int loop = 2;
   double DMS[3][3] = {{ 0. }}, DMP[3][3] = {{ 0. }};

   LOCK_MUTEX();

   effpot_(&loop, &rmb, &mg, &msb1sq, &msb2sq, &sxb, &cxb,
           &scalesq, &cotb, &vevS, &lamS, &svevS, &as, &DMS, &DMP);

   // Make appropriate substitutions for elements following 0907.4682
   // bottom of page 9
   std::swap(DMP[0][0], DMP[1][1]);
   std::swap(DMP[0][2], DMP[1][2]);

   Eigen::Matrix<double, 3, 3> result;
   result << DMP[0][0], DMP[0][1], DMP[0][2],
             DMP[1][0], DMP[1][1], DMP[1][2],
             DMP[2][0], DMP[2][1], DMP[2][2];

   return -result;
}

} // namespace nmssm_twoloophiggs
} // namespace flexiblesusy
