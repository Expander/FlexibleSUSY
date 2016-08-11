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
#include <cmath>
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_mssm)
#else
   #define LOCK_MUTEX()
#endif

namespace flexiblesusy {
namespace mssm_twoloophiggs {

static std::mutex mtx_mssm; /// locks MSSM fortran functions

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result;

   LOCK_MUTEX();

   ewsb2loop_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2,
              &mu, &tanb, &vev2, &gs, &result(0), &result(1));

   return -result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   LOCK_MUTEX();

   ddstad_(&mt2, &rmbsq, &mA2, &mst12, &mst22, &msb12, &msb22,
           &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
           &result(0), &result(1));

   return -result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2,
   double mu, double cotb, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result(tadpole_higgs_2loop_at_as_mssm(
      rmbsq, mg, msb12, msb22, sxb, cxb, scale2,
      mu, cotb, vev2, gs));

   std::swap(result(0), result(1));

   return result;
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msnusq, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   LOCK_MUTEX();

   tausqtad_(&mtau2, &mA2, &msnusq, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result(0), &result(1));

   return -result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   LOCK_MUTEX();

   dszhiggs_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
             &tanb, &vev2, &gs, &scheme,
             &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return -result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 2> result;

   LOCK_MUTEX();

   ddshiggs_(&mt2, &rmbsq, &mA2, &mst12, &mst22, &msb12, &msb22,
             &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
             &result(0,0), &result(0,1), &result(1,1));

   result(1,0) = result(0,1);

   return -result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
   double rmbsq, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result(self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
      rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs, scheme));

   std::swap(result(0,0), result(1,1));

   return result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
   double mtau2, double mA2, double msnusq, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   LOCK_MUTEX();

   tausqhiggs_(&mtau2, &mA2, &msnusq, &mstau12, &mstau22, &sintau,
               &costau, &scale2, &mu, &tanb, &vev2, &scheme,
               &result(0,0), &result(1,1), &result(0,1));

   result(1,0) = result(0,1);

   return -result;
}

double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   double result;

   LOCK_MUTEX();

   dszodd_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
           &tanb, &vev2, &gs, &result);

   return -result;
}

double self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   double result;

   LOCK_MUTEX();

   ddsodd_(&mt2, &rmbsq, &mA2, &mst12, &mst22, &msb12, &msb22,
           &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2, &result);

   return -result;
}

double self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
   double rmbsq, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs)
{
   return self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu,
      cotb, vev2, gs);
}

double self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
   double mtau2, double mA2, double msnusq, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   double result;

   LOCK_MUTEX();

   tausqodd_(&mtau2, &mA2, &msnusq, &mstau12, &mstau22, &sintau,
             &costau, &scale2, &mu, &tanb, &vev2, &result);

   return -result;
}

// self-energies without tadpoles

namespace {
template <typename T> T sqr(T a) { return a * a; }
}

Eigen::Matrix<double, 2, 2> rotate_scalar(
   double self_energy, double tanb)
{
   const double tanb2 = sqr(tanb);
   const double sinb = tanb / std::sqrt(1. + tanb2);
   const double cosb = 1. / std::sqrt(1. + tanb2);

   Eigen::Matrix<double, 2, 2> result;

   result(0,0) = self_energy * sqr(sinb);
   result(0,1) = - self_energy * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = self_energy * sqr(cosb);

   return result;
}

Eigen::Matrix<double, 2, 2> subtract_mssm_tadpoles_scalar(
   double self_energy, const Eigen::Matrix<double, 2, 1>& tadpoles,
   double tanb)
{
   return rotate_scalar(self_energy, tanb) + Eigen::Matrix<double, 2, 2>(tadpoles.asDiagonal());
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu,
         tanb, vev2, gs, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu,
      tanb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_as_mssm(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_at_mssm(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_at_at_mssm_with_tadpoles(
         mt2, rmbsq, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const double dMA = self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
      mt2, rmbsq, mA2, mst12, mst22, msb12, msb22,
      sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_at_mssm(
         mt2, rmbsq, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_ab_as_mssm_with_tadpoles(
         rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu,
         cotb, vev2, gs, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
      rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_ab_as_mssm(
         rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, 1./cotb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msnusq, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
         mtau2, mA2, msnusq, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msnusq, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msnusq, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 2> tM =
      subtract_mssm_tadpoles_scalar(dMA, tadpoles, tanb);

   return result + tM;
}

Eigen::Matrix<double, 2, 2> rotate_pseudoscalar(
   double self_energy, double tanb)
{
   const double tanb2 = sqr(tanb);
   const double sinb = tanb / std::sqrt(1. + tanb2);
   const double cosb = 1. / std::sqrt(1. + tanb2);

   Eigen::Matrix<double, 2, 2> result;

   // see hep-ph/0105096 Eq. (9)
   result(0,0) = self_energy * sqr(sinb);
   result(0,1) = self_energy * sinb * cosb;
   result(1,0) = result(0,1);
   result(1,1) = self_energy * sqr(cosb);

   return result;
}

Eigen::Matrix<double, 2, 2> subtract_mssm_tadpoles_pseudoscalar(
   double self_energy, const Eigen::Matrix<double, 2, 1>& tadpoles,
   double tanb)
{
   return rotate_pseudoscalar(self_energy, tanb) + Eigen::Matrix<double, 2, 2>(tadpoles.asDiagonal());
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   const double se = self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
      mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_as_mssm(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_at_at_mssm(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   const double se = self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
      mt2, rmbsq, mA2, mst12, mst22, msb12, msb22,
      sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_at_at_mssm(
         mt2, rmbsq, mA2, mst12, mst22, msb12, msb22,
         sxt, cxt, sxb, cxb, scale2, mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_ab_as_mssm(
   double rmbsq, double mg, double msb12, double msb22,
   double sxb, double cxb, double scale2, double mu,
   double cotb, double vev2, double gs)
{
   const double se = self_energy_pseudoscalar_2loop_ab_as_mssm_with_tadpoles(
      rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_ab_as_mssm(
         rmbsq, mg, msb12, msb22, sxb, cxb, scale2, mu, cotb, vev2, gs);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, 1./cotb);
}

Eigen::Matrix<double, 2, 2> self_energy_pseudoscalar_2loop_atau_atau_mssm(
   double mtau2, double mA2, double msnusq, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   const double se = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msnusq, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msnusq, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy
