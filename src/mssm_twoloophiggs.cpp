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
#include "dilog.hpp"
#include <cmath>
#include <utility>

#ifdef ENABLE_THREADS
   #include <mutex>
   #define LOCK_MUTEX() std::lock_guard<std::mutex> lg(mtx_mssm)
   namespace flexiblesusy {namespace mssm_twoloophiggs {
      static std::mutex mtx_mssm; /// locks MSSM fortran functions
   }}
#else
   #define LOCK_MUTEX()
#endif

namespace flexiblesusy {
namespace mssm_twoloophiggs {

namespace {

template <typename T> T sqr(T a) { return a * a; }

double phi(double x, double y, double z)
{
   using std::log;
   using gm2calc::dilog;

   const double u = x/z, v = y/z;
   const double lambda = std::sqrt(sqr(1 - u - v) - 4*u*v);
   const double xp = 0.5 * (1 + (u - v) - lambda);
   const double xm = 0.5 * (1 - (u - v) - lambda);

   return 1./lambda * (2*log(xp)*log(xm) - log(u)*log(v) -
                       2*(dilog(xp) + dilog(xm)) + M_PI*M_PI/3.);
}

/// limit st -> 0
Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm_st_0(
   double mt2, double mg, double mst12, double mst22,
   double /* sxt */, double /* cxt */, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   using std::sqrt;
   using std::atan;
   using std::log;
   using std::sin;
   using std::cos;

   constexpr double Pi = M_PI;
   constexpr double Pi4 = M_PI * M_PI * M_PI * M_PI;
   const double g = sqr(mg);
   const double q = scale2;
   const double t = mt2;
   const double T1 = mst12;
   const double T2 = mst22;
   const double v = std::sqrt(vev2);
   const double beta = std::atan(tanb);
   const double v2 = v * std::sin(beta);
   const double v1 = v * std::cos(beta);

   const double t1 = (sqr(gs)*mg*mt2*mu*(T1*T2*(5*(T1 - T2) + (-T1 +
      T2)*log(g/q)*log(t/q) + ((-g + t)*log(t/g) + T1*(-4 +
      log((g*t)/sqr(q))))*log(T1/q) + ((g - t)*log(t/g) - T2*(-4 +
      log((g*t)/sqr(q))))*log(T2/q)) + (sqr(g) + sqr(t - T1) - 2*g*(t
      + T1))*T2*phi(g,t,T1) - T1*(sqr(g) + sqr(t - T2) - 2*g*(t +
      T2))*phi(g,t,T2)))/(16.*Pi4*T1*(T1 - T2)*T2*tanb*sqr(v1));

   const double t2 = (sqr(gs)*mt2*(T1*T2*(-((T1 - T2)*(5*mg*mu + (2*g
      + 5*(-2*t + T1 + T2))*tanb)) + 6*t*(T1 - T2)*tanb*sqr(log(t/q))
      + (4*mg*mu*T1 + 2*(g + t + 2*T1)*(T1 - T2)*tanb + g*(-T1 +
      T2)*tanb*log(g/q) + mg*mu*((g - t)*log(t/g) -
      T1*log((g*t)/sqr(q))))*log(T1/q) + T1*(-T1 +
      T2)*tanb*sqr(log(T1/q)) + log(T2/q)*(-4*mg*mu*T2 + 2*(T1 -
      T2)*(g + t + 2*T2)*tanb + g*(-T1 + T2)*tanb*log(g/q) +
      mg*mu*((-g + t)*log(t/g) + T2*log((g*t)/sqr(q))) + T2*(-T1 +
      T2)*tanb*log(T2/q)) - (T1 - T2)*log(t/q)*(12*t*tanb - (mg*mu +
      2*g*tanb)*log(g/q) + (g + 2*t)*tanb*(log(T1/q) + log(T2/q)))) -
      T2*(mg*mu*(sqr(g) + sqr(t - T1) - 2*g*(t + T1)) - g*(g + t -
      T1)*(T1 - T2)*tanb)*phi(g,t,T1) + T1*(mg*mu*(sqr(g) + sqr(t -
      T2) - 2*g*(t + T2)) + g*(g + t - T2)*(T1 -
      T2)*tanb)*phi(g,t,T2)))/(16.*Pi4*T1*(T1 - T2)*T2*tanb*sqr(v2));

   Eigen::Matrix<double, 2, 1> result;
   result << t1, t2;

   return -result;
}

/// Pietro Slavich implementation
Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm_general(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   Eigen::Matrix<double, 2, 1> result;

   {
      LOCK_MUTEX();

      ewsb2loop_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2,
                 &mu, &tanb, &vev2, &gs, &result(0), &result(1));
   }

   if (!result.allFinite())
      result.setZero();

   return -result;
}

} // anonymous namespace

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_as_mssm(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2,
   double mu, double tanb, double vev2, double gs)
{
   if (std::abs(sxt) < 1e-7)
      return tadpole_higgs_2loop_at_as_mssm_st_0(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
   else
      return tadpole_higgs_2loop_at_as_mssm_general(
         mt2, mg, mst12, mst22, sxt, cxt, scale2, mu, tanb, vev2, gs);
}

Eigen::Matrix<double, 2, 1> tadpole_higgs_2loop_at_at_mssm(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   {
      LOCK_MUTEX();

      ddstad_(&mt2, &rmbsq, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
              &result(0), &result(1));
   }

   if (!result.allFinite())
      result.setZero();

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
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   Eigen::Matrix<double, 2, 1> result;

   {
      LOCK_MUTEX();

      tausqtad_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
                &costau, &scale2, &mu, &tanb, &vev2, &result(0), &result(1));
   }

   if (!result.allFinite())
      result.setZero();

   return -result;
}

Eigen::Matrix<double, 2, 2> self_energy_higgs_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   {
      LOCK_MUTEX();

      dszhiggs_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
                &tanb, &vev2, &gs, &scheme,
                &result(0,0), &result(1,1), &result(0,1));
   }

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

   {
      LOCK_MUTEX();

      ddshiggs_(&mt2, &rmbsq, &mA2, &mst12, &mst22, &msb12, &msb22,
                &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2,
                &result(0,0), &result(0,1), &result(1,1));
   }

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
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   Eigen::Matrix<double, 2, 2> result;

   {
      LOCK_MUTEX();

      tausqhiggs_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
                  &costau, &scale2, &mu, &tanb, &vev2, &scheme,
                  &result(0,0), &result(1,1), &result(0,1));
   }

   result(1,0) = result(0,1);

   return -result;
}

double self_energy_pseudoscalar_2loop_at_as_mssm_with_tadpoles(
   double mt2, double mg, double mst12, double mst22,
   double sxt, double cxt, double scale2, double mu,
   double tanb, double vev2, double gs)
{
   double result;

   {
      LOCK_MUTEX();

      dszodd_(&mt2, &mg, &mst12, &mst22, &sxt, &cxt, &scale2, &mu,
              &tanb, &vev2, &gs, &result);
   }

   return -result;
}

double self_energy_pseudoscalar_2loop_at_at_mssm_with_tadpoles(
   double mt2, double rmbsq, double mA2, double mst12,
   double mst22, double msb12, double msb22,
   double sxt, double cxt, double sxb, double cxb,
   double scale2, double mu, double tanb, double vev2)
{
   double result;

   {
      LOCK_MUTEX();

      ddsodd_(&mt2, &rmbsq, &mA2, &mst12, &mst22, &msb12, &msb22,
              &sxt, &cxt, &sxb, &cxb, &scale2, &mu, &tanb, &vev2, &result);
   }

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
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   double result;

   {
      LOCK_MUTEX();

      tausqodd_(&mtau2, &mA2, &msv2, &mstau12, &mstau22, &sintau,
                &costau, &scale2, &mu, &tanb, &vev2, &result);
   }

   return -result;
}

// self-energies without tadpoles

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
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2, int scheme)
{
   const Eigen::Matrix<double, 2, 2> result =
      self_energy_higgs_2loop_atau_atau_mssm_with_tadpoles(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2, scheme);

   const double dMA = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
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
   double mtau2, double mA2, double msv2, double mstau12,
   double mstau22, double sintau, double costau, double scale2,
   double mu, double tanb, double vev2)
{
   const double se = self_energy_pseudoscalar_2loop_atau_atau_mssm_with_tadpoles(
      mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
      mu, tanb, vev2);

   const Eigen::Matrix<double, 2, 1> tadpoles =
      tadpole_higgs_2loop_atau_atau_mssm(
         mtau2, mA2, msv2, mstau12, mstau22, sintau, costau, scale2,
         mu, tanb, vev2);

   return subtract_mssm_tadpoles_pseudoscalar(se, tadpoles, tanb);
}

} // namespace mssm_twoloophiggs
} // namespace flexiblesusy
