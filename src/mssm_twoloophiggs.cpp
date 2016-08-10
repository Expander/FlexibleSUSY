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

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"
#include "config.h"
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_mssm)
#else
   #define LOCK_MUTEX()
#endif

namespace flexiblesusy {

static std::mutex mtx_mssm; /// locks MSSM fortran functions

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq,
   double amu, double tanb, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result;

   LOCK_MUTEX();

   ewsb2loop_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq,
              &amu, &tanb, &vev2, &gs, &result(0), &result(1));

   return result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   double rmtsq, double rmbsq, double mAsq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   LOCK_MUTEX();

   ddstad_(&rmtsq, &rmbsq, &mAsq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
           &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2,
           &result(0), &result(1));

   return result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq,
   double amu, double cotbeta, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result(tadpole_higgs_2loop_at_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq,
      amu, cotbeta, vev2, gs));

   std::swap(result(0), result(1));

   return result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   double rmtausq, double mAsq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   LOCK_MUTEX();

   tausqtad_(&rmtausq, &mAsq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
             &costau, &scalesq, &amu, &tanb, &vev2, &result(0), &result(1));

   return result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double amu,
   double tanb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   LOCK_MUTEX();

   dszhiggs_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
             &tanb, &vev2, &gs, &scheme,
             &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   double rmtsq, double rmbsq, double fmasq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 2> result;

   LOCK_MUTEX();

   ddshiggs_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
             &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2,
             &result(0,0), &result(0,1), &result(1,1));

   result(1,0) = result(0,1);

   return result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double amu,
   double cotbeta, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result(self_energy_higgs_2loop_at_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu,
      cotbeta, vev2, gs, scheme));

   std::swap(result(0,0), result(1,1));

   return result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   double rmtausq, double fmasq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   LOCK_MUTEX();

   tausqhiggs_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
               &costau, &scalesq, &amu, &tanb, &vev2, &scheme,
               &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return result;
}

double self_energy_pseudoscalar_2loop_at_as_mssm(
   double rmtsq, double mg, double mst1sq, double mst2sq,
   double sxt, double cxt, double scalesq, double amu,
   double tanb, double vev2, double gs)
{
   double result;

   LOCK_MUTEX();

   dszodd_(&rmtsq, &mg, &mst1sq, &mst2sq, &sxt, &cxt, &scalesq, &amu,
           &tanb, &vev2, &gs, &result);

   return result;
}

double self_energy_pseudoscalar_2loop_at_at_mssm(
   double rmtsq, double rmbsq, double fmasq, double mst1sq,
   double mst2sq, double msb1sq, double msb2sq,
   double sxt, double cxt, double sxb, double cxb,
   double scalesq, double amu, double tanb, double vev2)
{
   double result;

   LOCK_MUTEX();

   ddsodd_(&rmtsq, &rmbsq, &fmasq, &mst1sq, &mst2sq, &msb1sq, &msb2sq,
           &sxt, &cxt, &sxb, &cxb, &scalesq, &amu, &tanb, &vev2, &result);

   return result;
}

double self_energy_pseudoscalar_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb1sq, double msb2sq,
   double sxb, double cxb, double scalesq, double amu,
   double cotbeta, double vev2, double gs)
{
   return self_energy_pseudoscalar_2loop_at_as_mssm(
      rmbsq, mg, msb1sq, msb2sq, sxb, cxb, scalesq, amu,
      cotbeta, vev2, gs);
}

double self_energy_pseudoscalar_2loop_atau_atau_mssm(
   double rmtausq, double fmasq, double msnusq, double mstau1sq,
   double mstau2sq, double sintau, double costau, double scalesq,
   double amu, double tanb, double vev2)
{
   double result;

   LOCK_MUTEX();

   tausqodd_(&rmtausq, &fmasq, &msnusq, &mstau1sq, &mstau2sq, &sintau,
             &costau, &scalesq, &amu, &tanb, &vev2, &result);

   return result;
}

} // namespace flexiblesusy
