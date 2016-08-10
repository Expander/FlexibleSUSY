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
#include "nmssm2loop.h"
#include "config.h"
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_nmssm)
#else
   #define LOCK_MUTEX()
#endif

namespace flexiblesusy {

static std::mutex mtx_nmssm; /// locks MSSM fortran functions

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_at_as_nmssm(
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

   return result;
}

Eigen::Matrix<double, 3, 3> self_energy_higgs_2loop_ab_as_nmssm(
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

   return result;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_at_as_nmssm(
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

   return result;
}

Eigen::Matrix<double, 3, 3> self_energy_pseudoscalar_2loop_ab_as_nmssm(
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

   return result;
}

} // namespace flexiblesusy
