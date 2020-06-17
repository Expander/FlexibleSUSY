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

// This file has been generated at Wed 17 Jun 2020 09:59:53
// with the script "tquark_to_cpp.m".

#include "mssm_twoloop_mt.hpp"
#include "dilog.hpp"
#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace mssm_twoloop_mt {

namespace {
   const double Pi  = 3.1415926535897932384626433832795;
   const double zt2 = 1.6449340668482264364724151666460;
   const double zt3 = 1.2020569031595942853997381615114;
   const double log2 = std::log(2.);

   template <typename T> T pow2(T x)  { return x*x; }
   template <typename T> T pow3(T x)  { return x*x*x; }
   template <typename T> T pow4(T x)  { return pow2(pow2(x)); }
   template <typename T> T pow5(T x)  { return x*pow4(x); }

   const double oneLoop = 1./pow2(4*Pi);
   const double twoLoop = pow2(oneLoop);

   template <typename T>
   bool is_zero(T a, T prec = std::numeric_limits<T>::epsilon())
   {
      return std::fabs(a) < prec;
   }

   template <typename T>
   bool is_equal(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      return is_zero(a - b, prec);
   }

   template <typename T>
   bool is_equal_rel(T a, T b, T prec = std::numeric_limits<T>::epsilon())
   {
      if (is_equal(a, b, std::numeric_limits<T>::epsilon()))
         return true;

      if (std::abs(a) < std::numeric_limits<T>::epsilon() ||
          std::abs(b) < std::numeric_limits<T>::epsilon())
         return false;

      return std::abs((a - b)/a) < prec;
   }

   /**
    * fin[] function from arXiv:hep-ph/0507139 .
    *
    * @param mm1 squared mass \f$m_1^2\f$
    * @param mm2 squared mass \f$m_2^2\f$
    * @param mmu squared renormalization scale
    *
    * @return fin(m12, m22)
    */
   double fin(double mm1, double mm2, double mmu)
   {
      const double PI = 3.14159265358979323846264338327950288;
      const double log1u = std::log(mm1/mmu);
      const double log2u = std::log(mm2/mmu);
      const double log12 = std::log(mm1/mm2);

      return (6*(mm1*log1u + mm2*log2u) +
         (-mm1 - mm2)*(7 + pow2(PI)/6.) +
         (mm1 - mm2)*(2*dilog(1 - mm1/mm2) +
            pow2(log12)/2.) +
         ((mm1 + mm2)*pow2(log12))/2. -
         2*(mm1*pow2(log1u) + mm2*pow2(log2u)))/2.;
   }

   /// shift gluino mass away from mst1 and mst2 if too close
   double shift_mg(double mg, double mst1, double mst2)
   {
      if (is_equal_rel(std::min(mst1, mst2), mg, 0.0003))
         return mg * 0.9995;

      if (is_equal_rel(std::max(mst1, mst2), mg, 0.0003))
         return mg * 1.0005;

      return mg;
   }

} // anonymous namespace

/// 1-loop QCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_qcd(const Parameters& pars)
{
   using std::log;
   const double g32 = pow2(pars.g3);
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);

   const double result = 6.666666666666667 - 4*log(mmt/mmu);

   return result * g32 * oneLoop;
}

/// 1-loop SUSY contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop_susy(const Parameters& pars)
{
   using std::log;
   const double g32    = pow2(pars.g3);
   const double mt     = pars.mt;
   const double mgl    = pars.mg;
   const double mmu    = pow2(pars.Q);
   const double mmgl   = pow2(pars.mg);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double SX     = 2*mt*pars.xt;
   const double s2t    = SX / (mmst1 - mmst2);

   if (is_equal(mmst1, mmst2, 1e-6) && is_equal(mmst1, mmgl, 1e-6)) {
      const double result = (4*log(mmgl/mmu))/3.;
      return result * g32 * oneLoop;
   }

   if (is_equal(mmst1, mmst2, 1e-6)) {
      const double result =
      (2*(4*mmgl*mmst2 - 4*mmgl*mmst2*log(mmst2/mmu) - 3*pow2(mmgl) + 2*log(
         mmgl/mmu)*pow2(mmgl) - pow2(mmst2) + 2*log(mmst2/mmu)*pow2(mmst2)))/(
         3.*pow2(-mmgl + mmst2));

      return result * g32 * oneLoop;
   }

   if (is_equal(mmgl, mmst1, 1e-6)) {
      const double result =
      -1 + (2*mmst2)/(3.*(-mmgl + mmst2)) + (2*log(mmgl/mmu))/3. - (4*mgl*mmgl*
         s2t*log(mmgl/mmu))/(3.*(-mmgl + mmst2)*mt) - (4*mgl*s2t*(1 + log(mmgl/
         mmu)))/(3.*mt) + (4*mgl*mmst2*s2t*log(mmst2/mmu))/(3.*(-mmgl + mmst2)*
         mt) + log(mmst2/mmu)*((4*mmst2)/(3.*(-mmgl + mmst2)) - (2*pow2(mmst2))
         /(3.*pow2(-mmgl + mmst2))) + (2*log(mmgl/mmu)*pow2(mmgl))/(3.*pow2(-
         mmgl + mmst2));

      return result * g32 * oneLoop;
   }

   if (is_equal(mmgl, mmst2, 1e-6)) {
      const double result =
      -1 + (2*mmst1)/(3.*(-mmgl + mmst1)) + (2*log(mmgl/mmu))/3. + (4*mgl*mmgl*
         s2t*log(mmgl/mmu))/(3.*(-mmgl + mmst1)*mt) + (4*mgl*s2t*(1 + log(mmgl/
         mmu)))/(3.*mt) - (4*mgl*mmst1*s2t*log(mmst1/mmu))/(3.*(-mmgl + mmst1)*
         mt) + log(mmst1/mmu)*((4*mmst1)/(3.*(-mmgl + mmst1)) - (2*pow2(mmst1))
         /(3.*pow2(-mmgl + mmst1))) + (2*log(mmgl/mmu)*pow2(mmgl))/(3.*pow2(-
         mmgl + mmst1));
      return result * g32 * oneLoop;

   }

   const double result =
   -2 + (2*mmst1)/(3.*(-mmgl + mmst1)) + (2*mmst2)/(3.*(-mmgl + mmst2)) + (((4*
      mgl*mmst1*s2t)/(3.*(-mmgl + mmst1)) - (4*mgl*mmst2*s2t)/(3.*(-mmgl +
      mmst2)))*log(mmgl/mmu))/mt - (4*mgl*mmst1*s2t*log(mmst1/mmu))/(3.*(-mmgl
      + mmst1)*mt) + (4*mgl*mmst2*s2t*log(mmst2/mmu))/(3.*(-mmgl + mmst2)*mt) +
      log(mmst1/mmu)*((4*mmst1)/(3.*(-mmgl + mmst1)) - (2*pow2(mmst1))/(3.*pow2
      (-mmgl + mmst1))) + log(mmst2/mmu)*((4*mmst2)/(3.*(-mmgl + mmst2)) - (2*
      pow2(mmst2))/(3.*pow2(-mmgl + mmst2))) + log(mmgl/mmu)*(
      1.3333333333333333 - (4*mmst1)/(3.*(-mmgl + mmst1)) - (4*mmst2)/(3.*(-
      mmgl + mmst2)) + (2*pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) + (2*pow2(mmst2
      ))/(3.*pow2(-mmgl + mmst2)));

   return result * g32 * oneLoop;
}

/// 1-loop full SQCD contributions to Delta Mt over mt [hep-ph/0210258]
double dMt_over_mt_1loop(const Parameters& pars)
{
   return dMt_over_mt_1loop_qcd(pars) + dMt_over_mt_1loop_susy(pars);
}

/// 2-loop QCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_qcd(const Parameters& pars)
{
   using std::log;
   const double g34 = pow4(pars.g3);
   const double mmt = pow2(pars.mt);
   const double mmu = pow2(pars.Q);

   const double result =
   111.72222222222223 + (16*zt2)/3. + (32*log2*zt2)/3. - (8*zt3)/3. - 82*log(
      mmt/mmu) + 22*pow2(log(mmt/mmu));

   return result * g34 * twoLoop;
}

/// 2-loop SUSY contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop_susy(const Parameters& pars)
{
   using std::log;
   const double g34    = pow4(pars.g3);
   const double Xt     = pars.xt;
   const double mt     = pars.mt;
   const double mgl    = shift_mg(pars.mg, pars.mst1, pars.mst2);
   const double mmu    = pow2(pars.Q);
   const double mmt    = pow2(pars.mt);
   const double mmgl   = pow2(mgl);
   const double mmst1  = pow2(pars.mst1);
   const double mmst2  = pow2(pars.mst2);
   const double mmsusy = pow2(pars.msusy);
   const double SX     = 2*mt*pars.xt;
   const double s2t    = SX / (mmst1 - mmst2);

   if (is_equal(mmst1, mmst2, mmt) && is_equal(mmst1, mmgl, mmt) &&
       is_equal(mmst1, mmsusy, mmt) && is_equal(std::abs(Xt), 0., 1e-1)) {
      const double result =
      32.31481481481482 + log(mmsusy/mmu)*(50.148148148148145 - (64*log(mmt/mmu
         ))/3.) + (62*pow2(log(mmsusy/mmu)))/9.;

      return result * g34 * twoLoop;
   }

   const double result =
   -10.055555555555555 - (56*mmst1)/(-mmgl + mmst1) - (896*mmst1)/(9.*(mmst1 -
      mmst2)) + (6*mmst2)/(-mmgl + mmst1) + (896*mmst2)/(9.*(mmst1 - mmst2)) -
      (2*mmst1)/(3.*(-mmgl + mmst2)) - (188*mmst2)/(3.*(-mmgl + mmst2)) + (60*
      mmsusy)/(-mmgl + mmst1) + (60*mmsusy)/(-mmgl + mmst2) + (8*zt2)/9. - (91*
      mmst1*zt2)/(9.*(-mmgl + mmst1)) - (128*mmst1*zt2)/(9.*(mmst1 - mmst2)) +
      (mmst2*zt2)/(-mmgl + mmst1) + (128*mmst2*zt2)/(9.*(mmst1 - mmst2)) + (
      mmst1*zt2)/(9.*(-mmgl + mmst2)) - (11*mmst2*zt2)/(-mmgl + mmst2) + (10*
      mmsusy*zt2)/(-mmgl + mmst1) + (10*mmsusy*zt2)/(-mmgl + mmst2) + (((160*
      mgl*mmsusy*s2t)/(3.*(-mmgl + mmst1)) - (160*mgl*mmsusy*s2t)/(3.*(-mmgl +
      mmst2)))*log(mmsusy/mmu))/mt + 82*log(mmt/mmu) + (-74 - (8*mmst1)/(3.*(-
      mmgl + mmst1)) - (8*mmst2)/(3.*(-mmgl + mmst2)))*log(mmt/mmu) + (((-16*
      mgl*mmst1*s2t)/(3.*(-mmgl + mmst1)) + (16*mgl*mmst2*s2t)/(3.*(-mmgl +
      mmst2)))*log(mmgl/mmu)*log(mmt/mmu))/mt + (16*mgl*mmst1*s2t*log(mmst1/mmu
      )*log(mmt/mmu))/(3.*(-mmgl + mmst1)*mt) - (16*mgl*mmst2*s2t*log(mmst2/mmu
      )*log(mmt/mmu))/(3.*(-mmgl + mmst2)*mt) - (20*log(mmsusy/mmu)*log(mmt/mmu
      ))/3. + (896*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (20*pow2
      (mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (128*zt2*pow2(mmst1))/(9.
      *(-mmgl + mmst1)*(mmst1 - mmst2)) + (8*zt2*pow2(mmst1))/(9.*(-mmgl +
      mmst1)*(-mmgl + mmst2)) + (fin(mmst2,mmgl,mmu)*((-16*mgl*s2t)/(9.*(-mmgl
      + mmst1)) + (88*mgl*s2t)/(9.*(-mmgl + mmst2)) - (4*mgl*mmst1*s2t)/(3.*
      pow2(-mmgl + mmst1)) + (4*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst1))))/mt +
      (fin(mmst1,mmsusy,mmu)*((40*mgl*mmst1*s2t)/(3.*pow2(-mmgl + mmst1)) - (40
      *mgl*mmsusy*s2t)/(3.*pow2(-mmgl + mmst1))))/mt + log(mmst1/mmu)*log(mmt/
      mmu)*(-0.6666666666666666 - (16*mmst1)/(3.*(-mmgl + mmst1)) + (8*pow2(
      mmst1))/(3.*pow2(-mmgl + mmst1))) - (2825*mmst1*mmst2)/(324.*pow2(-mmgl +
      mmst1)) - (27275*mmst1*mmsusy)/(324.*pow2(-mmgl + mmst1)) - (4*mmst1*
      mmst2*zt2)/(3.*pow2(-mmgl + mmst1)) - (130*mmst1*mmsusy*zt2)/(9.*pow2(-
      mmgl + mmst1)) - (40*mgl*mmst1*mmsusy*s2t*log(mmst1/mmu)*log(mmsusy/mmu))
      /(3.*mt*pow2(-mmgl + mmst1)) + (886*pow2(mmst1))/(9.*pow2(-mmgl + mmst1))
      + (233*mmst2*pow2(mmst1))/(324.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (
      1355*mmsusy*pow2(mmst1))/(324.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (
      44*zt2*pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) + (10*mmsusy*zt2*pow2(mmst1)
      )/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) - (896*pow2(mmst2))/(9.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) - (128*zt2*pow2(mmst2))/(9.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) - (233*mmst1*pow2(mmst2))/(216.*(mmst1 - mmst2)*
      pow2(-mmgl + mmst1)) + (233*pow2(mmst1)*pow2(mmst2))/(648.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmst2)) + (fin(mmst1,mmst2,mmu)*((4*mgl*s2t)/(9.*(-
      mmgl + mmst1)) - (4*mgl*s2t)/(9.*(-mmgl + mmst2)) + (4*mgl*mmst1*s2t)/(3.
      *pow2(-mmgl + mmst1)) - (4*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) + (4*
      mgl*mmst1*s2t)/(3.*pow2(-mmgl + mmst2)) - (4*mgl*mmst2*s2t)/(3.*pow2(-
      mmgl + mmst2))))/mt + (fin(mmst1,mmgl,mmu)*((-88*mgl*s2t)/(9.*(-mmgl +
      mmst1)) + (16*mgl*s2t)/(9.*(-mmgl + mmst2)) - (4*mgl*mmst1*s2t)/(3.*pow2(
      -mmgl + mmst2)) + (4*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + (fin(
      mmgl,mmsusy,mmu)*((-40*mgl*s2t)/(3.*(-mmgl + mmst1)) + (40*mgl*s2t)/(3.*(
      -mmgl + mmst2)) - (40*mgl*mmst1*s2t)/(3.*pow2(-mmgl + mmst1)) + (40*mgl*
      mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) + (40*mgl*mmst2*s2t)/(3.*pow2(-mmgl
      + mmst2)) - (40*mgl*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + (fin(
      mmst2,mmsusy,mmu)*((-40*mgl*mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) + (40*mgl
      *mmsusy*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + (log(mmgl/mmu)*log(mmsusy/
      mmu)*((-40*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst1)) + (40*mgl*mmsusy*s2t)/(3.
      *(-mmgl + mmst2)) + (40*mgl*mmst1*mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) -
      (40*mgl*mmst2*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2))))/mt + log(mmgl/mmu)*
      log(mmt/mmu)*(-13.333333333333334 + (16*mmst1)/(3.*(-mmgl + mmst1)) + (16
      *mmst2)/(3.*(-mmgl + mmst2)) - (8*pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) -
      (8*pow2(mmst2))/(3.*pow2(-mmgl + mmst2))) + log(mmst2/mmu)*log(mmt/mmu)*(
      -0.6666666666666666 - (16*mmst2)/(3.*(-mmgl + mmst2)) + (8*pow2(mmst2))/(
      3.*pow2(-mmgl + mmst2))) - (8*mmst1*mmst2)/pow2(-mmgl + mmst2) - (27275*
      mmst2*mmsusy)/(324.*pow2(-mmgl + mmst2)) - (4*mmst1*mmst2*zt2)/(3.*pow2(-
      mmgl + mmst2)) - (130*mmst2*mmsusy*zt2)/(9.*pow2(-mmgl + mmst2)) + (40*
      mgl*mmst2*mmsusy*s2t*log(mmst2/mmu)*log(mmsusy/mmu))/(3.*mt*pow2(-mmgl +
      mmst2)) + (886*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)) + (1355*mmsusy*pow2(
      mmst2))/(324.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) + (44*zt2*pow2(mmst2)
      )/(3.*pow2(-mmgl + mmst2)) + (10*mmsusy*zt2*pow2(mmst2))/(9.*(mmst2 -
      mmsusy)*pow2(-mmgl + mmst2)) - (515*mmst1*pow2(mmsusy))/(108.*(mmst1 -
      mmsusy)*pow2(-mmgl + mmst1)) - (20*mmst1*zt2*pow2(mmsusy))/(9.*(mmst1 -
      mmsusy)*pow2(-mmgl + mmst1)) - (515*mmst2*pow2(mmsusy))/(108.*(mmst2 -
      mmsusy)*pow2(-mmgl + mmst2)) - (20*mmst2*zt2*pow2(mmsusy))/(9.*(mmst2 -
      mmsusy)*pow2(-mmgl + mmst2)) + (95*pow2(mmst1)*pow2(mmsusy))/(162.*pow2(-
      mmgl + mmst1)*pow2(mmst1 - mmsusy)) + (10*zt2*pow2(mmst1)*pow2(mmsusy))/(
      9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) + (95*pow2(mmst2)*pow2(
      mmsusy))/(162.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) + (10*zt2*pow2(
      mmst2)*pow2(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (256
      *pow2(s2t))/9. + (63625*mmst1*pow2(s2t))/(1296.*(-mmgl + mmst1)) + (896*
      mmst1*pow2(s2t))/(9.*(mmst1 - mmst2)) - (10535*mmst2*pow2(s2t))/(3888.*(-
      mmgl + mmst1)) - (896*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*
      mmst2*pow2(s2t))/(648.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (56281*mmst1*
      pow2(s2t))/(1296.*(-mmgl + mmst2)) + (10967*mmst2*pow2(s2t))/(1296.*(-
      mmgl + mmst2)) - (85*mmst1*mmst2*pow2(s2t))/(972.*(mmst1 - mmst2)*(-mmgl
      + mmst2)) + (41*mmst1*zt2*pow2(s2t))/(9.*(-mmgl + mmst1)) + (128*mmst1*
      zt2*pow2(s2t))/(9.*(mmst1 - mmst2)) - (128*mmst2*zt2*pow2(s2t))/(9.*(
      mmst1 - mmst2)) + (mmst1*mmst2*zt2*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1
       - mmst2)) - (52*mmst1*zt2*pow2(s2t))/(9.*(-mmgl + mmst2)) - (11*mmst2*
      zt2*pow2(s2t))/(9.*(-mmgl + mmst2)) - (mmst1*mmst2*zt2*pow2(s2t))/(27.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) - (96683*pow2(mmst1)*pow2(s2t))/(972.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) + (52825*pow2(mmst1)*pow2(s2t))/(1296.*(-
      mmgl + mmst1)*(-mmgl + mmst2)) - (85*pow2(mmst1)*pow2(s2t))/(972.*(mmst1
      - mmst2)*(-mmgl + mmst2)) - (383*zt2*pow2(mmst1)*pow2(s2t))/(27.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) + (52*zt2*pow2(mmst1)*pow2(s2t))/(9.*(-mmgl +
      mmst1)*(-mmgl + mmst2)) - (zt2*pow2(mmst1)*pow2(s2t))/(27.*(mmst1 - mmst2
      )*(-mmgl + mmst2)) + (16*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst1))
      - (2*pow2(mmst1)*pow2(s2t))/pow2(-mmgl + mmst1) + (85*pow2(mmst2)*pow2(
      s2t))/(972.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (387233*pow2(mmst2)*pow2(
      s2t))/(3888.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (zt2*pow2(mmst2)*pow2(s2t
      ))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (383*zt2*pow2(mmst2)*pow2(s2t)
      )/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (167*mmst1*pow2(mmst2)*pow2(s2t
      ))/(3888.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (16*mmst1*mmst2*pow2(s2t
      ))/(9.*pow2(-mmgl + mmst2)) - (2*pow2(mmst2)*pow2(s2t))/pow2(-mmgl +
      mmst2) + (((-20*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst1)) + (20*mgl*mmsusy*s2t
      )/(3.*(-mmgl + mmst2)))*pow2(log(mmsusy/mmu)))/mt + (85*pow2(s2t)*pow3(
      mmst1))/(972.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl + mmst2)) + (zt2*
      pow2(s2t)*pow3(mmst1))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + fin(mmst1,mmsusy,mmu)*(20/(3.*(-mmgl + mmst1)) - (220*mmst1)/(
      9.*pow2(-mmgl + mmst1)) + (10*mmsusy)/pow2(-mmgl + mmst1) - (20*mmst1*
      mmsusy)/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (10*pow2(mmst1))/(9.*
      (mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (10*mmsusy*pow2(mmst1))/(9.*pow2(
      -mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (10*mmst1*pow2(mmsusy))/(9.*pow2(-
      mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (40*mmst1*mmsusy)/(3.*pow3(-mmgl +
      mmst1)) + (40*pow2(mmst1))/(3.*pow3(-mmgl + mmst1))) - (112*pow3(mmst1))/
      (3.*pow3(-mmgl + mmst1)) - (16*zt2*pow3(mmst1))/(3.*pow3(-mmgl + mmst1))
      + (271*pow3(mmst2))/(648.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (zt2*
      pow3(mmst2))/(9.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (7*mmst1*pow3(
      mmst2))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (mmst1*zt2*pow3(
      mmst2))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (169*pow2(s2t)*
      pow3(mmst2))/(1296.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (85*pow2(s2t)*
      pow3(mmst2))/(972.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) + (zt2*pow2(s2t)*
      pow3(mmst2))/(27.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (zt2*pow2(s2t)*
      pow3(mmst2))/(27.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) + fin(mmgl,mmsusy,
      mmu)*(10/(3.*(-mmgl + mmst1)) + 10/(3.*(-mmgl + mmst2)) + (10*mmst1)/pow2
      (-mmgl + mmst1) - (10*mmsusy)/pow2(-mmgl + mmst1) + (10*mmst2)/pow2(-mmgl
       + mmst2) - (10*mmsusy)/pow2(-mmgl + mmst2) + (40*mmst1*mmsusy)/(3.*pow3(
      -mmgl + mmst1)) - (40*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (40*mmst2*
      mmsusy)/(3.*pow3(-mmgl + mmst2)) - (40*pow2(mmst2))/(3.*pow3(-mmgl +
      mmst2))) + fin(mmst2,mmgl,mmu)*(7/(9.*(-mmgl + mmst1)) + 128/(9.*(mmst1 -
      mmst2)) - 5/(-mmgl + mmst2) + (8*mmst1)/(9.*(-mmgl + mmst1)*(-mmgl +
      mmst2)) - (128*mmst2)/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) + mmst1/pow2(-
      mmgl + mmst1) - mmst2/pow2(-mmgl + mmst1) + (12*mmst2)/pow2(-mmgl + mmst2
      ) + (11*pow2(s2t))/(9.*(-mmgl + mmst1)) - (128*pow2(s2t))/(9.*(mmst1 -
      mmst2)) - (11*pow2(s2t))/(9.*(-mmgl + mmst2)) + (26*mmst1*pow2(s2t))/(9.*
      (-mmgl + mmst1)*(-mmgl + mmst2)) + (128*mmst2*pow2(s2t))/(9.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) - (26*mmst1*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) +
      (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst1)) - (4*pow2(mmst1))/(3.*pow3(-mmgl
       + mmst1)) - (16*pow2(mmst2))/(3.*pow3(-mmgl + mmst2))) + fin(mmst1,mmgl,
      mmu)*(-37/(9.*(-mmgl + mmst1)) - 128/(9.*(mmst1 - mmst2)) + (128*mmst1)/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - 1/(9.*(-mmgl + mmst2)) + (8*mmst1)/
      (9.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (12*mmst1)/pow2(-mmgl + mmst1) -
      mmst1/pow2(-mmgl + mmst2) + mmst2/pow2(-mmgl + mmst2) + (5*pow2(s2t))/(3.
      *(-mmgl + mmst1)) + (128*pow2(s2t))/(9.*(mmst1 - mmst2)) - (128*mmst1*
      pow2(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (5*pow2(s2t))/(3.*(-
      mmgl + mmst2)) + (26*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)
      ) - (26*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) - (16*pow2(mmst1))/(3.*
      pow3(-mmgl + mmst1)) + (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst2)) - (4*pow2
      (mmst2))/(3.*pow3(-mmgl + mmst2))) + fin(mmst1,mmst2,mmu)*(2/(9.*(-mmgl +
      mmst1)) + 10/(9.*(-mmgl + mmst2)) - (8*mmst1)/(9.*(-mmgl + mmst1)*(-mmgl
      + mmst2)) - (7*mmst1)/(3.*pow2(-mmgl + mmst1)) + mmst2/pow2(-mmgl + mmst1
      ) + pow2(mmst2)/(9.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (mmst1*pow2(
      mmst2))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + mmst1/pow2(-mmgl +
      mmst2) - (7*mmst2)/(3.*pow2(-mmgl + mmst2)) - (67*pow2(s2t))/(54.*(-mmgl
      + mmst1)) - (mmst2*pow2(s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (65
      *pow2(s2t))/(54.*(-mmgl + mmst2)) - (mmst1*pow2(s2t))/(54.*(-mmgl + mmst1
      )*(-mmgl + mmst2)) + (mmst2*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + (26*mmst1*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) + (26*mmst2*pow2(
      s2t))/(9.*pow2(-mmgl + mmst2)) - (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst1))
      + (4*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + pow3(mmst2)/(9.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmst2)) - (4*mmst1*mmst2)/(3.*pow3(-mmgl + mmst2)) +
      (4*pow2(mmst2))/(3.*pow3(-mmgl + mmst2))) + fin(mmst2,mmsusy,mmu)*(20/(3.
      *(-mmgl + mmst2)) - (220*mmst2)/(9.*pow2(-mmgl + mmst2)) + (10*mmsusy)/
      pow2(-mmgl + mmst2) - (20*mmst2*mmsusy)/(9.*(mmst2 - mmsusy)*pow2(-mmgl +
      mmst2)) + (10*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2(-mmgl + mmst2)) + (
      10*mmsusy*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (
      10*mmst2*pow2(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (
      40*mmst2*mmsusy)/(3.*pow3(-mmgl + mmst2)) + (40*pow2(mmst2))/(3.*pow3(-
      mmgl + mmst2))) + log(mmgl/mmu)*log(mmsusy/mmu)*((10*mmsusy)/(-mmgl +
      mmst1) + (10*mmsusy)/(-mmgl + mmst2) - (70*mmst1*mmsusy)/(3.*pow2(-mmgl +
      mmst1)) - (70*mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) + (40*mmsusy*pow2(
      mmst1))/(3.*pow3(-mmgl + mmst1)) + (40*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl
       + mmst2))) + (pow2(log(mmgl/mmu))*((-10*mgl*mmst1*s2t)/(-mmgl + mmst1) -
      (2*mgl*mmst2*s2t)/(3.*(-mmgl + mmst1)) - (2*mgl*mmst1*s2t)/(9.*(-mmgl +
      mmst2)) + (98*mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) - (20*mgl*mmsusy*s2t)/(
      3.*(-mmgl + mmst1)) + (20*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst2)) + (8*mgl*
      s2t*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (2*mgl*mmst1*
      mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) + (20*mgl*mmst1*mmsusy*s2t)/(3.*pow2(
      -mmgl + mmst1)) - (218*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (
      14*mgl*mmst1*mmst2*s2t)/(9.*pow2(-mmgl + mmst2)) - (20*mgl*mmst2*mmsusy*
      s2t)/(3.*pow2(-mmgl + mmst2)) - (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst2)) + (70*mgl*s2t*pow2(mmst2))/(3.*pow2(-mmgl + mmst2)) - (8*mgl*s2t*
      pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) + (8*mgl*s2t*pow3(
      mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (8*mgl*s2t*pow3(mmst1)
      )/(9.*pow3(-mmgl + mmst1)) - (8*mgl*s2t*pow3(mmst2))/(9.*pow3(-mmgl +
      mmst2))))/mt + log(mmgl/mmu)*(64.88888888888889 + (98*mmst1)/(3.*(-mmgl +
      mmst1)) + (640*mmst1)/(9.*(mmst1 - mmst2)) - (2*mmst2)/(-mmgl + mmst1) -
      (640*mmst2)/(9.*(mmst1 - mmst2)) + (14*mmst1)/(3.*(-mmgl + mmst2)) + (118
      *mmst2)/(3.*(-mmgl + mmst2)) - (20*mmsusy)/(-mmgl + mmst1) - (20*mmsusy)/
      (-mmgl + mmst2) - (640*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2))
      - (20*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (14*mmst1*mmst2
      )/(3.*pow2(-mmgl + mmst1)) + (140*mmst1*mmsusy)/(3.*pow2(-mmgl + mmst1))
      - (116*pow2(mmst1))/pow2(-mmgl + mmst1) + (640*pow2(mmst2))/(9.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) + (38*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) + (
      140*mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) - (4*pow2(mmst1))/(9.*pow2(-
      mmgl + mmst2)) - (1048*pow2(mmst2))/(9.*pow2(-mmgl + mmst2)) + (256*pow2(
      s2t))/9. - (376*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)) - (640*mmst1*pow2(
      s2t))/(9.*(mmst1 - mmst2)) + (16*mmst2*pow2(s2t))/(9.*(-mmgl + mmst1)) +
      (640*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (40*mmst1*pow2(s2t))/(3.*(-
      mmgl + mmst2)) - (272*mmst2*pow2(s2t))/(9.*(-mmgl + mmst2)) + (640*pow2(
      mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (104*pow2(mmst1)
      *pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (32*mmst1*mmst2*pow2(
      s2t))/(9.*pow2(-mmgl + mmst1)) + (170*pow2(mmst1)*pow2(s2t))/(9.*pow2(-
      mmgl + mmst1)) - (640*pow2(mmst2)*pow2(s2t))/(9.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) - (4*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) + (28*pow2(
      mmst1)*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) + (22*pow2(mmst2)*pow2(s2t))/
      pow2(-mmgl + mmst2) + (4*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl +
      mmst1)) + (4*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (28*
      pow2(s2t)*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (28*
      pow2(s2t)*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (8*
      mmst2*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) - (80*mmsusy*pow2(mmst1))/(3.
      *pow3(-mmgl + mmst1)) + (16*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl +
      mmst1)) + (508*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + (4*pow2(s2t)*pow3(
      mmst1))/(3.*pow3(-mmgl + mmst1)) - (8*mmst1*pow2(mmst2))/(3.*pow3(-mmgl +
      mmst2)) - (80*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) + (16*mmst1*
      pow2(mmst2)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (508*pow3(mmst2))/(9.*
      pow3(-mmgl + mmst2)) + (4*pow2(s2t)*pow3(mmst2))/(3.*pow3(-mmgl + mmst2))
      ) - (112*pow3(mmst2))/(3.*pow3(-mmgl + mmst2)) - (16*zt2*pow3(mmst2))/(3.
      *pow3(-mmgl + mmst2)) - (95*mmst1*pow3(mmsusy))/(162.*pow2(-mmgl + mmst1)
      *pow2(mmst1 - mmsusy)) - (10*mmst1*zt2*pow3(mmsusy))/(9.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmsusy)) - (95*mmst2*pow3(mmsusy))/(162.*pow2(-mmgl +
      mmst2)*pow2(mmst2 - mmsusy)) - (10*mmst2*zt2*pow3(mmsusy))/(9.*pow2(-mmgl
       + mmst2)*pow2(mmst2 - mmsusy)) + log(mmst1/mmu)*log(mmsusy/mmu)*((80*
      mmst1*mmsusy)/(9.*pow2(-mmgl + mmst1)) + (10*mmsusy*pow2(mmst1))/(9.*(
      mmst1 - mmsusy)*pow2(-mmgl + mmst1)) - (20*mmst1*pow2(mmsusy))/(9.*(mmst1
       - mmsusy)*pow2(-mmgl + mmst1)) + (10*pow2(mmst1)*pow2(mmsusy))/(9.*pow2(
      -mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (40*mmsusy*pow2(mmst1))/(3.*pow3(-
      mmgl + mmst1)) - (10*mmst1*pow3(mmsusy))/(9.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmsusy))) + log(mmst2/mmu)*log(mmsusy/mmu)*((80*mmst2*mmsusy)/(9.
      *pow2(-mmgl + mmst2)) + (10*mmsusy*pow2(mmst2))/(9.*(mmst2 - mmsusy)*pow2
      (-mmgl + mmst2)) - (20*mmst2*pow2(mmsusy))/(9.*(mmst2 - mmsusy)*pow2(-
      mmgl + mmst2)) + (10*pow2(mmst2)*pow2(mmsusy))/(9.*pow2(-mmgl + mmst2)*
      pow2(mmst2 - mmsusy)) - (40*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2))
      - (10*mmst2*pow3(mmsusy))/(9.*pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)))
      + pow2(log(mmsusy/mmu))*(3.3333333333333335 + (5*mmsusy)/(-mmgl + mmst1)
      + (5*mmsusy)/(-mmgl + mmst2) - (65*mmst1*mmsusy)/(9.*pow2(-mmgl + mmst1))
      + (5*mmsusy*pow2(mmst1))/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) - (65*
      mmst2*mmsusy)/(9.*pow2(-mmgl + mmst2)) + (5*mmsusy*pow2(mmst2))/(9.*(
      mmst2 - mmsusy)*pow2(-mmgl + mmst2)) - (10*mmst1*pow2(mmsusy))/(9.*(mmst1
       - mmsusy)*pow2(-mmgl + mmst1)) - (10*mmst2*pow2(mmsusy))/(9.*(mmst2 -
      mmsusy)*pow2(-mmgl + mmst2)) + (5*pow2(mmst1)*pow2(mmsusy))/(9.*pow2(-
      mmgl + mmst1)*pow2(mmst1 - mmsusy)) + (5*pow2(mmst2)*pow2(mmsusy))/(9.*
      pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy)) - (5*mmst1*pow3(mmsusy))/(9.*
      pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (5*mmst2*pow3(mmsusy))/(9.*
      pow2(-mmgl + mmst2)*pow2(mmst2 - mmsusy))) + log(mmsusy/mmu)*(
      21.11111111111111 - (40*mmsusy)/(-mmgl + mmst1) - (40*mmsusy)/(-mmgl +
      mmst2) + (170*mmst1*mmsusy)/(3.*pow2(-mmgl + mmst1)) - (10*mmsusy*pow2(
      mmst1))/(3.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (170*mmst2*mmsusy)/(
      3.*pow2(-mmgl + mmst2)) - (10*mmsusy*pow2(mmst2))/(3.*(mmst2 - mmsusy)*
      pow2(-mmgl + mmst2)) + (20*mmst1*pow2(mmsusy))/(3.*(mmst1 - mmsusy)*pow2(
      -mmgl + mmst1)) + (20*mmst2*pow2(mmsusy))/(3.*(mmst2 - mmsusy)*pow2(-mmgl
       + mmst2)) - (10*pow2(mmst1)*pow2(mmsusy))/(3.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmsusy)) - (10*pow2(mmst2)*pow2(mmsusy))/(3.*pow2(-mmgl + mmst2)*
      pow2(mmst2 - mmsusy)) + (10*mmst1*pow3(mmsusy))/(3.*pow2(-mmgl + mmst1)*
      pow2(mmst1 - mmsusy)) + (10*mmst2*pow3(mmsusy))/(3.*pow2(-mmgl + mmst2)*
      pow2(mmst2 - mmsusy))) + (log(mmst1/mmu)*((64*mgl*s2t)/9. - (32*mgl*mmst1
      *s2t)/(9.*(-mmgl + mmst1)) + (64*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1)) - (
      64*mgl*mmst1*mmst2*s2t)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (40*mgl*
      mmst1*s2t)/(9.*(-mmgl + mmst2)) - (8*mgl*s2t*pow2(mmst1))/(9.*(-mmgl +
      mmst1)*(-mmgl + mmst2)) + (8*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst1)
      ) + (80*mgl*mmst1*mmsusy*s2t)/(3.*pow2(-mmgl + mmst1)) - (40*mgl*s2t*pow2
      (mmst1))/pow2(-mmgl + mmst1) + (64*mgl*s2t*pow2(mmst2))/(9.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) + (32*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) -
      (16*mgl*mmst2*pow3(s2t))/(3.*(-mmgl + mmst1)) + (64*mgl*mmst1*mmst2*pow3(
      s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (16*mgl*mmst1*pow3(s2t))/(9.
      *(-mmgl + mmst2)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst1
      )) + (16*mgl*pow2(mmst1)*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) - (64*mgl*
      pow2(mmst2)*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2))))/mt + (pow2(
      log(mmst1/mmu))*((22*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1)) - (16*mgl*mmst2*
      s2t)/(9.*(-mmgl + mmst1)) + (16*mgl*mmst1*mmst2*s2t)/(3.*(-mmgl + mmst1)*
      (mmst1 - mmst2)) + (2*mgl*mmst1*s2t)/(3.*(-mmgl + mmst2)) - (2*mgl*mmst1*
      mmst2*s2t)/(3.*pow2(-mmgl + mmst1)) - (20*mgl*mmst1*mmsusy*s2t)/(3.*pow2(
      -mmgl + mmst1)) + (62*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (16
      *mgl*s2t*pow2(mmst2))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (8*mgl*s2t*
      pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) - (16*mgl*mmst1*pow3(s2t))/(9.*(-
      mmgl + mmst1)) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (16*mgl*
      mmst1*mmst2*pow3(s2t))/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*
      pow2(mmst1)*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) + (16*mgl*pow2(mmst2)*
      pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2))))/mt + ((-220*mgl*mmst1*
      s2t)/(3.*(-mmgl + mmst1)) - (56*mgl*mmst2*s2t)/(3.*(-mmgl + mmst1)) + (32
      *mgl*mmst1*mmst2*s2t)/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (8*mgl*mmst1
      *s2t)/(-mmgl + mmst2) + (220*mgl*mmst2*s2t)/(3.*(-mmgl + mmst2)) - (80*
      mgl*mmsusy*s2t)/(-mmgl + mmst1) + (80*mgl*mmsusy*s2t)/(-mmgl + mmst2) - (
      88*mgl*mmst1*s2t*zt2)/(9.*(-mmgl + mmst1)) - (44*mgl*mmst2*s2t*zt2)/(9.*(
      -mmgl + mmst1)) + (32*mgl*mmst1*mmst2*s2t*zt2)/(9.*(-mmgl + mmst1)*(mmst1
       - mmst2)) + (4*mgl*mmst1*s2t*zt2)/(3.*(-mmgl + mmst2)) + (88*mgl*mmst2*
      s2t*zt2)/(9.*(-mmgl + mmst2)) - (40*mgl*mmsusy*s2t*zt2)/(3.*(-mmgl +
      mmst1)) + (40*mgl*mmsusy*s2t*zt2)/(3.*(-mmgl + mmst2)) - (32*mgl*s2t*pow2
      (mmst2))/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (32*mgl*s2t*zt2*pow2(
      mmst2))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*mmst1*pow3(s2t))/(
      9.*(-mmgl + mmst1)) + (112*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (
      32*mgl*mmst1*mmst2*pow3(s2t))/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*
      mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.*
      (-mmgl + mmst2)) + (32*mgl*mmst2*zt2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (
      32*mgl*mmst1*mmst2*zt2*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) +
      (32*mgl*pow2(mmst2)*pow3(s2t))/(3.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (32
      *mgl*zt2*pow2(mmst2)*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)))/mt
      + (log(mmst1/mmu)*log(mmst2/mmu)*((8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1))
      - (32*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1)) - (8*mgl*mmst1*s2t)/(9.*(-mmgl
      + mmst2)) - (32*mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) + (8*mgl*s2t*pow2(
      mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (4*mgl*mmst1*mmst2*s2t)/(
      3.*pow2(-mmgl + mmst1)) - (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)
      ) - (32*mgl*s2t*pow2(mmst2))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (32*
      mgl*s2t*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (4*mgl*mmst1*
      mmst2*s2t)/(9.*pow2(-mmgl + mmst2)) - (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-
      mmgl + mmst2)) - (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl +
      mmst1)) + (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)
      ) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) + (16*mgl*mmst1*pow3(
      s2t))/(9.*(-mmgl + mmst2)) + (32*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)
      ) + (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) + (32*mgl*
      pow2(mmst2)*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (32*mgl*
      pow2(mmst2)*pow3(s2t))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (16*mgl*
      mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst2))))/mt + (log(mmgl/mmu)*log
      (mmst1/mmu)*((-64*mgl*s2t)/9. - (12*mgl*mmst1*s2t)/(-mmgl + mmst1) - (32*
      mgl*mmst1*mmst2*s2t)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (20*mgl*mmst1
      *s2t)/(9.*(-mmgl + mmst2)) + (32*mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) - (8
      *mgl*s2t*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (52*mgl*s2t*
      pow2(mmst1))/(3.*pow2(-mmgl + mmst1)) + (32*mgl*s2t*pow2(mmst2))/(9.*(
      mmst1 - mmst2)*(-mmgl + mmst2)) - (4*mgl*mmst1*mmst2*s2t)/(9.*pow2(-mmgl
      + mmst2)) + (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl + mmst2)) + (8*mgl*s2t
      *pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (8*mgl*s2t*pow3(
      mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (16*mgl*s2t*pow3(mmst1
      ))/(9.*pow3(-mmgl + mmst1)) + (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst1
      )) + (32*mgl*mmst1*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2))
      - (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) - (32*mgl*mmst2*pow3(s2t)
      )/(9.*(-mmgl + mmst2)) + (16*mgl*pow2(mmst1)*pow3(s2t))/(9.*pow2(-mmgl +
      mmst1)) - (32*mgl*pow2(mmst2)*pow3(s2t))/(9.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst2))))/mt +
      (log(mmgl/mmu)*log(mmst2/mmu)*((64*mgl*s2t)/9. - (8*mgl*mmst1*s2t)/(9.*(-
      mmgl + mmst1)) - (4*mgl*mmst2*s2t)/(3.*(-mmgl + mmst1)) + (32*mgl*mmst1*
      mmst2*s2t)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (8*mgl*mmst1*s2t)/(9.*(
      -mmgl + mmst2)) + (68*mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) - (8*mgl*s2t*
      pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (4*mgl*mmst1*mmst2*
      s2t)/(3.*pow2(-mmgl + mmst1)) + (8*mgl*s2t*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst1)) - (32*mgl*s2t*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) +
      (8*mgl*mmst1*mmst2*s2t)/(9.*pow2(-mmgl + mmst2)) + (8*mgl*s2t*pow2(mmst1)
      )/(9.*pow2(-mmgl + mmst2)) - (148*mgl*s2t*pow2(mmst2))/(9.*pow2(-mmgl +
      mmst2)) + (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)
      ) - (8*mgl*s2t*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (
      16*mgl*s2t*pow3(mmst2))/(9.*pow3(-mmgl + mmst2)) + (16*mgl*mmst2*pow3(s2t
      ))/(9.*(-mmgl + mmst1)) - (32*mgl*mmst1*mmst2*pow3(s2t))/(9.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) + (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) -
      (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst1)) + (32*mgl*pow2(
      mmst2)*pow3(s2t))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (16*mgl*pow2(
      mmst2)*pow3(s2t))/(9.*pow2(-mmgl + mmst2))))/mt + (log(mmst2/mmu)*((-64*
      mgl*s2t)/9. + (8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1)) + (112*mgl*mmst2*s2t
      )/(9.*(-mmgl + mmst1)) - (64*mgl*mmst1*mmst2*s2t)/(9.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) - (8*mgl*mmst1*s2t)/(9.*(-mmgl + mmst2)) + (8*mgl*mmst2*
      s2t)/(3.*(-mmgl + mmst2)) + (8*mgl*s2t*pow2(mmst1))/(9.*(-mmgl + mmst1)*(
      -mmgl + mmst2)) + (64*mgl*s2t*pow2(mmst2))/(9.*(-mmgl + mmst1)*(mmst1 -
      mmst2)) - (8*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) - (80*mgl*
      mmst2*mmsusy*s2t)/(3.*pow2(-mmgl + mmst2)) + (40*mgl*s2t*pow2(mmst2))/
      pow2(-mmgl + mmst2) - (80*mgl*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) + (64
      *mgl*mmst1*mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*
      mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) - (32*mgl*mmst2*pow3(s2t))/(9.*
      (-mmgl + mmst2)) - (64*mgl*pow2(mmst2)*pow3(s2t))/(9.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) + (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst2))
      - (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*pow2(-mmgl + mmst2))))/mt + (pow2(
      log(mmst2/mmu))*((-22*mgl*mmst2*s2t)/(9.*(-mmgl + mmst1)) + (16*mgl*mmst1
      *mmst2*s2t)/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (10*mgl*mmst2*s2t)/(9.
      *(-mmgl + mmst2)) - (16*mgl*s2t*pow2(mmst2))/(9.*(-mmgl + mmst1)*(mmst1 -
      mmst2)) + (32*mgl*s2t*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) +
      (2*mgl*mmst1*mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) + (20*mgl*mmst2*mmsusy*
      s2t)/(3.*pow2(-mmgl + mmst2)) - (62*mgl*s2t*pow2(mmst2))/(9.*pow2(-mmgl +
      mmst2)) - (8*mgl*s2t*pow3(mmst2))/(9.*pow3(-mmgl + mmst2)) + (16*mgl*
      mmst2*pow3(s2t))/(9.*(-mmgl + mmst1)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(
      9.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (16*mgl*mmst2*pow3(s2t))/(9.*(-mmgl
       + mmst2)) + (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*(-mmgl + mmst1)*(mmst1 -
      mmst2)) - (32*mgl*pow2(mmst2)*pow3(s2t))/(9.*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + (16*mgl*pow2(mmst2)*pow3(s2t))/(9.*pow2(-mmgl + mmst2))))/mt +
      (log(mmgl/mmu)*((544*mgl*mmst1*s2t)/(9.*(-mmgl + mmst1)) + (8*mgl*mmst2*
      s2t)/(3.*(-mmgl + mmst1)) - (8*mgl*mmst1*s2t)/(3.*(-mmgl + mmst2)) - (544
      *mgl*mmst2*s2t)/(9.*(-mmgl + mmst2)) + (80*mgl*mmsusy*s2t)/(3.*(-mmgl +
      mmst1)) - (80*mgl*mmsusy*s2t)/(3.*(-mmgl + mmst2)) - (8*mgl*mmst1*mmst2*
      s2t)/(3.*pow2(-mmgl + mmst1)) - (80*mgl*mmst1*mmsusy*s2t)/(3.*pow2(-mmgl
      + mmst1)) + (40*mgl*s2t*pow2(mmst1))/pow2(-mmgl + mmst1) + (8*mgl*mmst1*
      mmst2*s2t)/(3.*pow2(-mmgl + mmst2)) + (80*mgl*mmst2*mmsusy*s2t)/(3.*pow2(
      -mmgl + mmst2)) - (40*mgl*s2t*pow2(mmst2))/pow2(-mmgl + mmst2) - (16*mgl*
      mmst1*pow3(s2t))/(9.*(-mmgl + mmst1)) - (16*mgl*mmst2*pow3(s2t))/(9.*(-
      mmgl + mmst1)) + (16*mgl*mmst1*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16*mgl*
      mmst2*pow3(s2t))/(9.*(-mmgl + mmst2)) + (16*mgl*mmst1*mmst2*pow3(s2t))/(
      9.*pow2(-mmgl + mmst1)) - (16*mgl*pow2(mmst1)*pow3(s2t))/(9.*pow2(-mmgl +
      mmst1)) - (16*mgl*mmst1*mmst2*pow3(s2t))/(9.*pow2(-mmgl + mmst2)) + (16*
      mgl*pow2(mmst2)*pow3(s2t))/(9.*pow2(-mmgl + mmst2))))/mt + log(mmgl/mmu)*
      log(mmst1/mmu)*(7.111111111111111 - (59*mmst1)/(9.*(-mmgl + mmst1)) + (13
      *mmst1)/(9.*(-mmgl + mmst2)) + (4*pow2(mmst1))/(9.*(-mmgl + mmst1)*(-mmgl
       + mmst2)) - (109*pow2(mmst1))/(9.*pow2(-mmgl + mmst1)) - (29*mmst1*mmst2
      )/(9.*pow2(-mmgl + mmst2)) - (4*pow2(mmst1))/(3.*pow2(-mmgl + mmst2)) + (
      64*pow2(s2t))/9. - (67*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)) - (61*mmst1*
      pow2(s2t))/(9.*(-mmgl + mmst2)) + (13*pow2(mmst1)*pow2(s2t))/((-mmgl +
      mmst1)*(-mmgl + mmst2)) - (97*pow2(mmst1)*pow2(s2t))/(9.*pow2(-mmgl +
      mmst1)) + (56*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) - (4*pow2(
      mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) + (16*pow3(mmst1))/(9.*(-mmgl
      + mmst1)*pow2(-mmgl + mmst2)) - (32*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl +
      mmst2)*pow2(-mmgl + mmst1)) + (16*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl +
      mmst1)*pow2(-mmgl + mmst2)) + (116*pow3(mmst1))/(9.*pow3(-mmgl + mmst1))
      + (160*pow2(s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + (4*mmst1*pow2(
      mmst2))/(3.*pow3(-mmgl + mmst2)) - (32*mmst2*pow2(mmst1)*pow2(s2t))/(9.*
      pow3(-mmgl + mmst2)) - (16*mmst1*pow2(mmst2)*pow2(s2t))/(3.*pow3(-mmgl +
      mmst2)) - (32*pow2(s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst2)) - (4*pow4(
      mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (4*pow2(s2t)*pow4(
      mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (32*pow2(s2t)*pow4
      (mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl + mmst1)) + (32*pow2(s2t)*pow4(
      mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl + mmst2)) - (4*pow4(mmst1))/(3.*
      pow4(-mmgl + mmst1)) - (20*pow2(s2t)*pow4(mmst1))/(3.*pow4(-mmgl + mmst1)
      )) + (271*pow4(mmst2))/(648.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (
      zt2*pow4(mmst2))/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + log(mmst2
      /mmu)*(-5 + (4*mmst1)/(9.*(-mmgl + mmst1)) - (4*mmst2)/(-mmgl + mmst1) -
      (128*mmst2)/(9.*(mmst1 - mmst2)) - (4*mmst1)/(9.*(-mmgl + mmst2)) + (178*
      mmst2)/(9.*(-mmgl + mmst2)) + (16*mmst1*mmst2)/(3.*pow2(-mmgl + mmst1)) +
      (128*pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (14*mmst1*mmst2)
      /(9.*pow2(-mmgl + mmst2)) - (515*mmst2*mmsusy)/(27.*pow2(-mmgl + mmst2))
      + (4*pow2(mmst1))/(9.*pow2(-mmgl + mmst2)) + (248*pow2(mmst2))/(9.*pow2(-
      mmgl + mmst2)) - (25*mmsusy*pow2(mmst2))/(27.*(mmst2 - mmsusy)*pow2(-mmgl
       + mmst2)) - (5*mmst2*pow2(mmsusy))/(9.*(mmst2 - mmsusy)*pow2(-mmgl +
      mmst2)) + (40*pow2(mmst2)*pow2(mmsusy))/(27.*pow2(-mmgl + mmst2)*pow2(
      mmst2 - mmsusy)) - (64*pow2(s2t))/9. - (4385*mmst1*pow2(s2t))/(324.*(-
      mmgl + mmst1)) + (871*mmst2*pow2(s2t))/(324.*(-mmgl + mmst1)) + (128*
      mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) - (11*mmst1*mmst2*pow2(s2t))/(162.*
      (-mmgl + mmst1)*(mmst1 - mmst2)) + (3809*mmst1*pow2(s2t))/(324.*(-mmgl +
      mmst2)) + (8993*mmst2*pow2(s2t))/(324.*(-mmgl + mmst2)) + (11*mmst1*mmst2
      *pow2(s2t))/(162.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (11*pow2(mmst1)*pow2
      (s2t))/(162.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (5681*pow2(mmst1)*pow2(
      s2t))/(324.*(-mmgl + mmst1)*(-mmgl + mmst2)) + (11*pow2(mmst1)*pow2(s2t))
      /(162.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (16*mmst1*mmst2*pow2(s2t))/(9.*
      pow2(-mmgl + mmst1)) + (32*pow2(mmst1)*pow2(s2t))/(9.*pow2(-mmgl + mmst1)
      ) - (5*pow2(mmst2)*pow2(s2t))/(108.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (
      1531*pow2(mmst2)*pow2(s2t))/(108.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (4*
      mmst1*mmst2*pow2(s2t))/pow2(-mmgl + mmst2) + (4*pow2(mmst1)*pow2(s2t))/(
      9.*pow2(-mmgl + mmst2)) - (182*pow2(mmst2)*pow2(s2t))/(9.*pow2(-mmgl +
      mmst2)) - (4*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (11*
      pow2(s2t)*pow3(mmst1))/(162.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + (32*pow2(s2t)*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl +
      mmst1)) - (4*pow2(s2t)*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) - pow3(mmst2)/(3.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (mmst1*
      pow3(mmst2))/(3.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (11*pow2(s2t)
      *pow3(mmst2))/(162.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (11*pow2(s2t)*
      pow3(mmst2))/(162.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) + (8*mmst1*pow2(
      mmst2))/(3.*pow3(-mmgl + mmst2)) + (80*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl
       + mmst2)) - (16*mmst1*pow2(mmst2)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) -
      (220*pow3(mmst2))/(9.*pow3(-mmgl + mmst2)) - (4*pow2(s2t)*pow3(mmst2))/(
      3.*pow3(-mmgl + mmst2)) - (40*mmst2*pow3(mmsusy))/(27.*pow2(-mmgl + mmst2
      )*pow2(mmst2 - mmsusy)) - pow4(mmst2)/(3.*pow2(-mmgl + mmst1)*pow2(mmst1
      - mmst2))) + log(mmst1/mmu)*(-5 + (62*mmst1)/(3.*(-mmgl + mmst1)) + (128*
      mmst1)/(9.*(mmst1 - mmst2)) - (44*mmst1)/(9.*(-mmgl + mmst2)) - (128*pow2
      (mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (4*pow2(mmst1))/(3.*(-
      mmgl + mmst1)*(-mmgl + mmst2)) - (41*mmst1*mmst2)/(27.*pow2(-mmgl + mmst1
      )) - (515*mmst1*mmsusy)/(27.*pow2(-mmgl + mmst1)) + (244*pow2(mmst1))/(9.
      *pow2(-mmgl + mmst1)) - (13*mmst2*pow2(mmst1))/(27.*(mmst1 - mmst2)*pow2(
      -mmgl + mmst1)) - (25*mmsusy*pow2(mmst1))/(27.*(mmst1 - mmsusy)*pow2(-
      mmgl + mmst1)) + (13*mmst1*pow2(mmst2))/(18.*(mmst1 - mmst2)*pow2(-mmgl +
      mmst1)) - (13*pow2(mmst1)*pow2(mmst2))/(54.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmst2)) + (16*mmst1*mmst2)/(3.*pow2(-mmgl + mmst2)) - (5*mmst1*
      pow2(mmsusy))/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) + (40*pow2(mmst1)
      *pow2(mmsusy))/(27.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) - (64*pow2(
      s2t))/9. + (2959*mmst1*pow2(s2t))/(162.*(-mmgl + mmst1)) - (128*mmst1*
      pow2(s2t))/(9.*(mmst1 - mmst2)) - (16*mmst2*pow2(s2t))/(9.*(-mmgl + mmst1
      )) - (2*mmst1*mmst2*pow2(s2t))/(81.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (
      1973*mmst1*pow2(s2t))/(162.*(-mmgl + mmst2)) - (7*mmst2*pow2(s2t))/(162.*
      (-mmgl + mmst2)) + (11*mmst1*mmst2*pow2(s2t))/(162.*(mmst1 - mmst2)*(-
      mmgl + mmst2)) + (2293*pow2(mmst1)*pow2(s2t))/(162.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) - (893*pow2(mmst1)*pow2(s2t))/(162.*(-mmgl + mmst1)*(-
      mmgl + mmst2)) + (11*pow2(mmst1)*pow2(s2t))/(162.*(mmst1 - mmst2)*(-mmgl
      + mmst2)) + (32*mmst1*mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst1)) - (62*
      pow2(mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst1)) - (29*pow2(mmst2)*pow2(
      s2t))/(324.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (2*pow2(mmst2)*pow2(s2t))/
      (81.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (7*mmst1*pow2(mmst2)*pow2(s2t))/(
      324.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (16*mmst1*mmst2*pow2(s2t))/(
      3.*pow2(-mmgl + mmst2)) - (32*pow2(mmst1)*pow2(s2t))/(9.*pow2(-mmgl +
      mmst2)) - (4*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (11*
      pow2(s2t)*pow3(mmst1))/(162.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl +
      mmst2)) - (4*pow2(s2t)*pow3(mmst1))/(9.*(-mmgl + mmst2)*pow2(-mmgl +
      mmst1)) + (32*pow2(s2t)*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) + (8*mmst2*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (80*mmsusy*
      pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) - (16*mmst2*pow2(mmst1)*pow2(s2t))/
      (9.*pow3(-mmgl + mmst1)) - (220*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) - (
      4*pow2(s2t)*pow3(mmst1))/(3.*pow3(-mmgl + mmst1)) - (5*pow3(mmst2))/(54.*
      (mmst1 - mmst2)*pow2(-mmgl + mmst1)) + (mmst1*pow3(mmst2))/(3.*pow2(-mmgl
       + mmst1)*pow2(mmst1 - mmst2)) - (29*pow2(s2t)*pow3(mmst2))/(324.*(-mmgl
      + mmst1)*pow2(mmst1 - mmst2)) + (11*pow2(s2t)*pow3(mmst2))/(162.*(-mmgl +
      mmst2)*pow2(mmst1 - mmst2)) - (40*mmst1*pow3(mmsusy))/(27.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmsusy)) - (5*pow4(mmst2))/(54.*pow2(-mmgl + mmst1)*
      pow2(mmst1 - mmst2))) + pow2(log(mmst1/mmu))*(0.3333333333333333 - (35*
      mmst1)/(18.*(-mmgl + mmst1)) - (64*mmst1)/(9.*(mmst1 - mmst2)) + mmst1/(
      2.*(-mmgl + mmst2)) + (64*pow2(mmst1))/(9.*(-mmgl + mmst1)*(mmst1 - mmst2
      )) + (mmst1*mmst2)/(2.*pow2(-mmgl + mmst1)) + (40*mmst1*mmsusy)/(9.*pow2(
      -mmgl + mmst1)) - (25*pow2(mmst1))/(6.*pow2(-mmgl + mmst1)) + (5*mmsusy*
      pow2(mmst1))/(9.*(mmst1 - mmsusy)*pow2(-mmgl + mmst1)) - (2*mmst1*mmst2)/
      (3.*pow2(-mmgl + mmst2)) - (10*mmst1*pow2(mmsusy))/(9.*(mmst1 - mmsusy)*
      pow2(-mmgl + mmst1)) + (5*pow2(mmst1)*pow2(mmsusy))/(9.*pow2(-mmgl +
      mmst1)*pow2(mmst1 - mmsusy)) - (25*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1))
      + (64*mmst1*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*pow2(s2t))/(
      54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (13*mmst1*pow2(s2t))/(9.*(-mmgl +
      mmst2)) - (mmst1*mmst2*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl + mmst2)) -
      (383*pow2(mmst1)*pow2(s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (13*
      pow2(mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (pow2(mmst1
      )*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (185*pow2(mmst1)*
      pow2(s2t))/(18.*pow2(-mmgl + mmst1)) + (pow2(mmst2)*pow2(s2t))/(54.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) - (pow2(mmst2)*pow2(s2t))/(54.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) + (pow2(s2t)*pow3(mmst1))/(54.*(-mmgl + mmst1)*(
      mmst1 - mmst2)*(-mmgl + mmst2)) - (2*mmst2*pow2(mmst1))/(3.*pow3(-mmgl +
      mmst1)) - (20*mmsusy*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (10*pow3(
      mmst1))/(3.*pow3(-mmgl + mmst1)) - (8*pow2(s2t)*pow3(mmst1))/pow3(-mmgl +
      mmst1) + pow3(mmst2)/(18.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (mmst1*
      pow3(mmst2))/(18.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (pow2(s2t)*
      pow3(mmst2))/(54.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (pow2(s2t)*pow3(
      mmst2))/(54.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) - (5*mmst1*pow3(mmsusy)
      )/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmsusy)) + (2*pow4(mmst1))/(3.*
      pow4(-mmgl + mmst1)) + (10*pow2(s2t)*pow4(mmst1))/(3.*pow4(-mmgl + mmst1)
      ) + pow4(mmst2)/(18.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2))) + log(
      mmst1/mmu)*log(mmst2/mmu)*((4*mmst1)/(9.*(-mmgl + mmst1)) - (4*mmst1)/(9.
      *(-mmgl + mmst2)) - (4*pow2(mmst1))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2))
      + (mmst1*mmst2)/pow2(-mmgl + mmst1) - (4*pow2(mmst1))/(9.*pow2(-mmgl +
      mmst1)) + (17*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) + (4*pow2(mmst1))/(3.
      *pow2(-mmgl + mmst2)) - (19*mmst1*pow2(s2t))/(9.*(-mmgl + mmst1)) + (16*
      mmst2*pow2(s2t))/(9.*(-mmgl + mmst1)) + (mmst1*mmst2*pow2(s2t))/(27.*(-
      mmgl + mmst1)*(mmst1 - mmst2)) + (35*mmst1*pow2(s2t))/(9.*(-mmgl + mmst2)
      ) - (mmst1*mmst2*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (pow2
      (mmst1)*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (91*pow2(mmst1
      )*pow2(s2t))/(9.*(-mmgl + mmst1)*(-mmgl + mmst2)) - (pow2(mmst1)*pow2(s2t
      ))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (32*mmst1*mmst2*pow2(s2t))/(9.
      *pow2(-mmgl + mmst1)) + (20*pow2(mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst1
      )) + (pow2(mmst2)*pow2(s2t))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (
      pow2(mmst2)*pow2(s2t))/(27.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (56*mmst1*
      mmst2*pow2(s2t))/(9.*pow2(-mmgl + mmst2)) + (4*pow2(mmst1)*pow2(s2t))/(3.
      *pow2(-mmgl + mmst2)) - (16*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) + (pow2(s2t)*pow3(mmst1))/(27.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-
      mmgl + mmst2)) + (32*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl + mmst2)*pow2(-
      mmgl + mmst1)) - (16*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl + mmst1)*pow2(-
      mmgl + mmst2)) - (4*mmst2*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (16*
      mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst1)) - (32*pow2(s2t)*
      pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + pow3(mmst2)/(9.*(mmst1 - mmst2)*
      pow2(-mmgl + mmst1)) - (mmst1*pow3(mmst2))/(9.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmst2)) + (pow2(s2t)*pow3(mmst2))/(27.*(-mmgl + mmst1)*pow2(mmst1
       - mmst2)) - (pow2(s2t)*pow3(mmst2))/(27.*(-mmgl + mmst2)*pow2(mmst1 -
      mmst2)) - (4*mmst1*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) + (32*mmst2*pow2
      (mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (16*mmst1*pow2(mmst2)*pow2(
      s2t))/(3.*pow3(-mmgl + mmst2)) + (32*pow2(s2t)*pow3(mmst1))/(9.*pow3(-
      mmgl + mmst2)) + (4*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) + (4*pow2(s2t)*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) - (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl +
      mmst1)) - (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl +
      mmst2)) + pow4(mmst2)/(9.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2))) + log
      (mmgl/mmu)*log(mmst2/mmu)*(7.111111111111111 - (4*mmst1)/(9.*(-mmgl +
      mmst1)) + mmst2/(-mmgl + mmst1) + (4*mmst1)/(9.*(-mmgl + mmst2)) - (55*
      mmst2)/(9.*(-mmgl + mmst2)) + (4*pow2(mmst1))/(9.*(-mmgl + mmst1)*(-mmgl
      + mmst2)) - (7*mmst1*mmst2)/(3.*pow2(-mmgl + mmst1)) + (4*pow2(mmst1))/(
      9.*pow2(-mmgl + mmst1)) - (8*mmst1*mmst2)/(9.*pow2(-mmgl + mmst2)) - (4*
      pow2(mmst1))/(3.*pow2(-mmgl + mmst2)) - (113*pow2(mmst2))/(9.*pow2(-mmgl
      + mmst2)) + (64*pow2(s2t))/9. + (5*mmst1*pow2(s2t))/(-mmgl + mmst1) - (16
      *mmst2*pow2(s2t))/(9.*(-mmgl + mmst1)) - (5*mmst1*pow2(s2t))/(-mmgl +
      mmst2) - (112*mmst2*pow2(s2t))/(9.*(-mmgl + mmst2)) + (13*pow2(mmst1)*
      pow2(s2t))/((-mmgl + mmst1)*(-mmgl + mmst2)) + (32*mmst1*mmst2*pow2(s2t))
      /(9.*pow2(-mmgl + mmst1)) - (20*pow2(mmst1)*pow2(s2t))/(3.*pow2(-mmgl +
      mmst1)) + (8*mmst1*mmst2*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) - (4*pow2(
      mmst1)*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) - (37*pow2(mmst2)*pow2(s2t))/(
      9.*pow2(-mmgl + mmst2)) + (16*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl
       + mmst2)) - (32*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl + mmst2)*pow2(-mmgl +
      mmst1)) + (16*pow2(s2t)*pow3(mmst1))/(3.*(-mmgl + mmst1)*pow2(-mmgl +
      mmst2)) + (4*mmst2*pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) - (16*mmst2*pow2
      (mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst1)) + (32*pow2(s2t)*pow3(mmst1))/
      (9.*pow3(-mmgl + mmst1)) - (32*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-
      mmgl + mmst2)) - (32*mmst1*pow2(mmst2)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)
      ) - (32*pow2(s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst2)) + (116*pow3(mmst2
      ))/(9.*pow3(-mmgl + mmst2)) + (128*pow2(s2t)*pow3(mmst2))/(9.*pow3(-mmgl
      + mmst2)) - (4*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2))
      - (4*pow2(s2t)*pow4(mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2))
      + (32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl + mmst1)) + (
      32*pow2(s2t)*pow4(mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl + mmst2)) - (4*
      pow4(mmst2))/(3.*pow4(-mmgl + mmst2)) - (20*pow2(s2t)*pow4(mmst2))/(3.*
      pow4(-mmgl + mmst2))) + pow2(log(mmgl/mmu))*(-9.555555555555555 - (211*
      mmst1)/(18.*(-mmgl + mmst1)) - (64*mmst1)/(3.*(mmst1 - mmst2)) + mmst2/(
      2.*(-mmgl + mmst1)) + (64*mmst2)/(3.*(mmst1 - mmst2)) - (31*mmst1)/(18.*(
      -mmgl + mmst2)) - (251*mmst2)/(18.*(-mmgl + mmst2)) + (5*mmsusy)/(-mmgl +
      mmst1) + (5*mmsusy)/(-mmgl + mmst2) + (64*pow2(mmst1))/(3.*(-mmgl + mmst1
      )*(mmst1 - mmst2)) + (4*pow2(mmst1))/(3.*(-mmgl + mmst1)*(-mmgl + mmst2))
      - (7*mmst1*mmst2)/(6.*pow2(-mmgl + mmst1)) - (35*mmst1*mmsusy)/(3.*pow2(-
      mmgl + mmst1)) + (821*pow2(mmst1))/(18.*pow2(-mmgl + mmst1)) - (64*pow2(
      mmst2))/(3.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (5*mmst1*mmst2)/(18.*pow2(
      -mmgl + mmst2)) - (35*mmst2*mmsusy)/(3.*pow2(-mmgl + mmst2)) + (4*pow2(
      mmst1))/(3.*pow2(-mmgl + mmst2)) + (829*pow2(mmst2))/(18.*pow2(-mmgl +
      mmst2)) - (128*pow2(s2t))/9. + (15*mmst1*pow2(s2t))/(-mmgl + mmst1) + (64
      *mmst1*pow2(s2t))/(3.*(mmst1 - mmst2)) - (64*mmst2*pow2(s2t))/(3.*(mmst1
      - mmst2)) - (7*mmst1*pow2(s2t))/(9.*(-mmgl + mmst2)) + (128*mmst2*pow2(
      s2t))/(9.*(-mmgl + mmst2)) - (64*pow2(mmst1)*pow2(s2t))/(3.*(-mmgl +
      mmst1)*(mmst1 - mmst2)) - (65*pow2(mmst1)*pow2(s2t))/(9.*(-mmgl + mmst1)*
      (-mmgl + mmst2)) + (pow2(mmst1)*pow2(s2t))/(2.*pow2(-mmgl + mmst1)) + (64
      *pow2(mmst2)*pow2(s2t))/(3.*(mmst1 - mmst2)*(-mmgl + mmst2)) - (8*mmst1*
      mmst2*pow2(s2t))/(3.*pow2(-mmgl + mmst2)) + (4*pow2(mmst1)*pow2(s2t))/(3.
      *pow2(-mmgl + mmst2)) - (37*pow2(mmst2)*pow2(s2t))/(6.*pow2(-mmgl + mmst2
      )) - (16*pow3(mmst1))/(9.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (32*pow2
      (s2t)*pow3(mmst1))/(3.*(-mmgl + mmst2)*pow2(-mmgl + mmst1)) - (16*pow2(
      s2t)*pow3(mmst1))/(3.*(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (2*mmst2*
      pow2(mmst1))/(3.*pow3(-mmgl + mmst1)) + (20*mmsusy*pow2(mmst1))/(3.*pow3(
      -mmgl + mmst1)) - (242*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) - (88*pow2(
      s2t)*pow3(mmst1))/(9.*pow3(-mmgl + mmst1)) + (2*mmst1*pow2(mmst2))/(3.*
      pow3(-mmgl + mmst2)) + (20*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl + mmst2)) +
      (32*mmst2*pow2(mmst1)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (32*mmst1*
      pow2(mmst2)*pow2(s2t))/(9.*pow3(-mmgl + mmst2)) + (32*pow2(s2t)*pow3(
      mmst1))/(9.*pow3(-mmgl + mmst2)) - (242*pow3(mmst2))/(9.*pow3(-mmgl +
      mmst2)) - (56*pow2(s2t)*pow3(mmst2))/(9.*pow3(-mmgl + mmst2)) + (4*pow4(
      mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2)) + (4*pow2(s2t)*pow4(
      mmst1))/(9.*pow2(-mmgl + mmst1)*pow2(-mmgl + mmst2)) - (32*pow2(s2t)*pow4
      (mmst1))/(9.*(-mmgl + mmst2)*pow3(-mmgl + mmst1)) - (32*pow2(s2t)*pow4(
      mmst1))/(9.*(-mmgl + mmst1)*pow3(-mmgl + mmst2)) + (2*pow4(mmst1))/(3.*
      pow4(-mmgl + mmst1)) + (10*pow2(s2t)*pow4(mmst1))/(3.*pow4(-mmgl + mmst1)
      ) + (2*pow4(mmst2))/(3.*pow4(-mmgl + mmst2)) + (10*pow2(s2t)*pow4(mmst2))
      /(3.*pow4(-mmgl + mmst2))) + pow2(log(mmst2/mmu))*(0.3333333333333333 +
      mmst2/(2.*(-mmgl + mmst1)) + (64*mmst2)/(9.*(mmst1 - mmst2)) - (35*mmst2)
      /(18.*(-mmgl + mmst2)) - (2*mmst1*mmst2)/(3.*pow2(-mmgl + mmst1)) - (64*
      pow2(mmst2))/(9.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (mmst1*mmst2)/(2.*
      pow2(-mmgl + mmst2)) + (40*mmst2*mmsusy)/(9.*pow2(-mmgl + mmst2)) - (25*
      pow2(mmst2))/(6.*pow2(-mmgl + mmst2)) + (5*mmsusy*pow2(mmst2))/(9.*(mmst2
       - mmsusy)*pow2(-mmgl + mmst2)) - (10*mmst2*pow2(mmsusy))/(9.*(mmst2 -
      mmsusy)*pow2(-mmgl + mmst2)) + (5*pow2(mmst2)*pow2(mmsusy))/(9.*pow2(-
      mmgl + mmst2)*pow2(mmst2 - mmsusy)) + (13*mmst1*pow2(s2t))/(9.*(-mmgl +
      mmst1)) - (64*mmst2*pow2(s2t))/(9.*(mmst1 - mmst2)) + (mmst1*mmst2*pow2(
      s2t))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) - (13*mmst1*pow2(s2t))/(9.*(-
      mmgl + mmst2)) - (38*mmst2*pow2(s2t))/(9.*(-mmgl + mmst2)) - (mmst1*mmst2
      *pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl + mmst2)) + (pow2(mmst1)*pow2(s2t
      ))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)) + (13*pow2(mmst1)*pow2(s2t))/(9.
      *(-mmgl + mmst1)*(-mmgl + mmst2)) - (pow2(mmst1)*pow2(s2t))/(54.*(mmst1 -
      mmst2)*(-mmgl + mmst2)) + (pow2(mmst2)*pow2(s2t))/(54.*(-mmgl + mmst1)*(
      mmst1 - mmst2)) + (383*pow2(mmst2)*pow2(s2t))/(54.*(mmst1 - mmst2)*(-mmgl
       + mmst2)) + (185*pow2(mmst2)*pow2(s2t))/(18.*pow2(-mmgl + mmst2)) + (
      pow2(s2t)*pow3(mmst1))/(54.*(-mmgl + mmst1)*(mmst1 - mmst2)*(-mmgl +
      mmst2)) + pow3(mmst2)/(18.*(mmst1 - mmst2)*pow2(-mmgl + mmst1)) - (mmst1*
      pow3(mmst2))/(18.*pow2(-mmgl + mmst1)*pow2(mmst1 - mmst2)) + (pow2(s2t)*
      pow3(mmst2))/(54.*(-mmgl + mmst1)*pow2(mmst1 - mmst2)) - (pow2(s2t)*pow3(
      mmst2))/(54.*(-mmgl + mmst2)*pow2(mmst1 - mmst2)) - (2*mmst1*pow2(mmst2))
      /(3.*pow3(-mmgl + mmst2)) - (20*mmsusy*pow2(mmst2))/(3.*pow3(-mmgl +
      mmst2)) + (10*pow3(mmst2))/(3.*pow3(-mmgl + mmst2)) - (8*pow2(s2t)*pow3(
      mmst2))/pow3(-mmgl + mmst2) - (5*mmst2*pow3(mmsusy))/(9.*pow2(-mmgl +
      mmst2)*pow2(mmst2 - mmsusy)) + pow4(mmst2)/(18.*pow2(-mmgl + mmst1)*pow2(
      mmst1 - mmst2)) + (2*pow4(mmst2))/(3.*pow4(-mmgl + mmst2)) + (10*pow2(s2t
      )*pow4(mmst2))/(3.*pow4(-mmgl + mmst2)));

   return result * g34 * twoLoop;
}

/// 2-loop full SQCD contributions to Delta Mt over mt [hep-ph/0507139]
double dMt_over_mt_2loop(const Parameters& pars)
{
   return dMt_over_mt_2loop_qcd(pars) + dMt_over_mt_2loop_susy(pars);
}

} // namespace mssm_twoloop_mt
} // namespace flexiblesusy
