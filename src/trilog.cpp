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

#include "trilog.hpp"
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {
   template <typename T> T pow2(T x) noexcept { return x*x; }

   // converts -0.0 to 0.0
   template <typename T>
   std::complex<T> clog(std::complex<T> z) noexcept {
      if (std::real(z) == T(0)) { z.real(T(0)); }
      if (std::imag(z) == T(0)) { z.imag(T(0)); }
      return std::log(z);
   }

   template <typename T>
   bool is_close(const std::complex<T>& a, T b, T eps)
   {
      return std::abs(std::real(a) - b) < eps && std::abs(std::imag(a)) < eps;
   }

   template <typename T>
   std::complex<T> cadd(T a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(a + std::real(b), std::imag(b));
   }

   template <typename T>
   std::complex<T> cadd(const std::complex<T>& a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(std::real(a) + std::real(b),
                             std::imag(a) + std::imag(b));
   }

   template <typename T>
   std::complex<T> cmul(const std::complex<T>& a, T b) noexcept
   {
      return std::complex<T>(std::real(a) * b, std::imag(a) * b);
   }

   template <typename T>
   std::complex<T> cmul(const std::complex<T>& a, const std::complex<T>& b) noexcept
   {
      return std::complex<T>(
         std::real(a) * std::real(b) - std::imag(a) * std::imag(b),
         std::real(a) * std::imag(b) + std::imag(a) * std::real(b));
   }

} // anonymous namespace

/**
 * @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$
 * @param z complex argument
 * @return \f$\mathrm{Li}_3(z)\f$
 */
std::complex<double> trilog(const std::complex<double>& z) noexcept
{
   const double eps   = 10.0*std::numeric_limits<double>::epsilon();
   const double PI    = 3.141592653589793;
   const double PI2   = PI*PI;
   const double zeta2 = 1.644934066848226;
   const double zeta3 = 1.202056903159594;
   const double bf[18] = {
      1.0                  , -3.0/8.0              ,
      17.0/216.0           , -5.0/576.0            ,
      1.296296296296296e-04,  8.101851851851851e-05,
     -3.419357160853759e-06, -1.328656462585034e-06,
      8.660871756109851e-08,  2.526087595532039e-08,
     -2.144694468364064e-09, -5.140110622012978e-10,
      5.249582114600829e-11,  1.088775440663631e-11,
     -1.277939609449369e-12, -2.369824177308745e-13,
      3.104357887965462e-14,  5.261758629912506e-15
   };

   if (is_close(z, 0.0, eps)) {
      return { 0.0, 0.0 };
   }
   if (is_close(z, 1.0, eps)) {
      return { zeta3, 0.0 };
   }
   if (is_close(z, -1.0, eps)) {
      return { -0.75*zeta3, 0.0 };
   }
   if (is_close(z, 0.5, eps)) {
      const double ln2  = 0.6931471805599453; // ln(2)
      const double ln23 = 0.3330246519889295; // ln(2)^3
      return { (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.0, 0.0 };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.0) { // |log(z)| < 1
      const auto u  = std::complex<double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta3 + u*(zeta2 - u2/12.0);
      const auto c1 = 0.25 * (3.0 - 2.0*clog(-u));

      const double cs[7] = {
         -3.472222222222222e-03, 1.157407407407407e-05,
         -9.841899722852104e-08, 1.148221634332745e-09,
         -1.581572499080917e-11, 2.419500979252515e-13,
         -3.982897776989488e-15
      };

      return
         cadd(c0,
         cmul(u2, cadd(c1,
         cmul(u2, cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cadd(cs[4],
         cmul(u2, cadd(cs[5],
         cmul(u2, cs[6]))))))))))))))));
   }

   std::complex<double> u(0.0, 0.0), rest(0.0, 0.0);

   if (az <= 1.0) {
      u = -clog(1.0 - z);
   } else { // az > 1
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<double>(lnz, arg); // clog(-z)
      u = -clog(1.0 - 1.0/z);
      rest = -lmz*(pow2(lmz)/6.0 + zeta2);
   }

   return
      cadd(rest,
      cmul(u, cadd(bf[0],
      cmul(u, cadd(bf[1],
      cmul(u, cadd(bf[2],
      cmul(u, cadd(bf[3],
      cmul(u, cadd(bf[4],
      cmul(u, cadd(bf[5],
      cmul(u, cadd(bf[6],
      cmul(u, cadd(bf[7],
      cmul(u, cadd(bf[8],
      cmul(u, cadd(bf[9],
      cmul(u, cadd(bf[10],
      cmul(u, cadd(bf[11],
      cmul(u, cadd(bf[12],
      cmul(u, cadd(bf[13],
      cmul(u, cadd(bf[14],
      cmul(u, cadd(bf[15],
      cmul(u, cadd(bf[16],
      cmul(u, bf[17]))))))))))))))))))))))))))))))))))));
}

/**
 * @brief Complex trilogarithm \f$\mathrm{Li}_3(z)\f$ with long double precision
 * @param z complex argument
 * @return \f$\mathrm{Li}_3(z)\f$
 */
std::complex<long double> trilog(const std::complex<long double>& z) noexcept
{
   const long double eps   = 10.0L*std::numeric_limits<long double>::epsilon();
   const long double PI    = 3.14159265358979323846264338327950288L;
   const long double PI2   = PI*PI;
   const long double zeta2 = 1.64493406684822643647241516664602519L;
   const long double zeta3 = 1.20205690315959428539973816151144999L;
   const long double bf[45] = {
      1.0L,
     -3.0L/8.0L,
      17.0L/216.0L,
     -5.0L/576.0L,
      7.0L/54000.0L,
      7.0L/86400.0L,
     -3.41935716085375949321527552820069827e-06L,
     -1.32865646258503401360544217687074830e-06L,
      8.66087175610985134794658604182413706e-08L,
      2.52608759553203997648442092886537331e-08L,
     -2.14469446836406476093388507573649032e-09L,
     -5.14011062201297891533581769272004962e-10L,
      5.24958211460082943639408880855807284e-11L,
      1.08877544066363183753729715704249107e-11L,
     -1.27793960944936953055818317540722120e-12L,
     -2.36982417730874520997977788101244891e-13L,
      3.10435788796546229428475327046556211e-14L,
      5.26175862991250608413183925112250061e-15L,
     -7.53847954994926536599250143226771028e-16L,
     -1.18623225777522852530825009512459322e-16L,
      1.83169799654913833820892731212815349e-17L,
      2.70681710318373501514907347126169436e-18L,
     -4.45543389782963882643263099217632212e-19L,
     -6.23754849225569465036532224739838641e-20L,
      1.08515215348745349131365609968642833e-20L,
      1.44911748660360819307349049665275324e-21L,
     -2.64663397544589903347408911861443741e-22L,
     -3.38976534885101047219258165860814078e-23L,
      6.46404773360331088903253098219534234e-24L,
      7.97583448960241242420922272590502795e-25L,
     -1.58091787902874833559211176293826770e-25L,
     -1.88614997296228681931102253988531956e-26L,
      3.87155366384184733039971271888313319e-27L,
      4.48011750023456073048653898320511684e-28L,
     -9.49303387191183612641753676027699150e-29L,
     -1.06828138090773812240182143033807908e-29L,
      2.33044789361030518600785199019281371e-30L,
      2.55607757265197540805635698286695865e-31L,
     -5.72742160613725968447274458033057100e-32L,
     -6.13471321379642358258549296897773326e-33L,
      1.40908086040689448401268688489421700e-33L,
      1.47642223976665341443182801167106626e-34L,
     -3.47010516489959160555004020312910903e-35L,
     -3.56210662409746357967357370318293608e-36L,
      8.55369656823692105754731289124468101e-37L
   };

   if (is_close(z, 0.0L, eps)) {
      return { 0.0L, 0.0L };
   }
   if (is_close(z, 1.0L, eps)) {
      return { zeta3, 0.0L };
   }
   if (is_close(z, -1.0L, eps)) {
      return { -0.75L*zeta3, 0.0L };
   }
   if (is_close(z, 0.5L, eps)) {
      const long double ln2  = 0.693147180559945309417232121458176568L; // ln(2)
      const long double ln23 = 0.333024651988929479718853582611730544L; // ln(2)^3
      return { (-2*PI2*ln2 + 4*ln23 + 21*zeta3)/24.0L, 0.0L };
   }

   const auto az  = std::abs(z);
   const auto pz  = std::arg(z);
   const auto lnz = std::log(az);

   if (pow2(lnz) + pow2(pz) < 1.0L) { // |log(z)| < 1
      const auto u  = std::complex<long double>(lnz, pz); // clog(z)
      const auto u2 = u*u;
      const auto c0 = zeta3 + u*(zeta2 - u2/12.0L);
      const auto c1 = 0.25L * (3.0L - 2.0L*clog(-u));

      const long double cs[20] = {
        -3.47222222222222222222222222222222222e-03L,
         1.15740740740740740740740740740740741e-05L,
        -9.84189972285210380448475686570924666e-08L,
         1.14822163433274544385655496766607878e-09L,
        -1.58157249908091658933409775160616911e-11L,
         2.41950097925251519452732701564998016e-13L,
        -3.98289777698948774786517290926462002e-15L,
         6.92336661830592905806820954095065870e-17L,
        -1.25527223044997727545846570912655367e-18L,
         2.35375400276846523056441171414060379e-20L,
        -4.53639890345868701844750708901700830e-22L,
         8.94516967039264316712031170773304472e-24L,
        -1.79828400469549627172020247141015426e-25L,
         3.67549976479373844433604733912674099e-27L,
        -7.62080797156479522953948500963765478e-29L,
         1.60004196436948597517376392257325602e-30L,
        -3.39676114756037558792312060520851852e-32L,
         7.28227228675776469531725636144432664e-34L,
        -1.57502264795800348718497893940378261e-35L,
         3.43354009248058933588797212016527038e-37L
      };

      return
         cadd(c0,
         cmul(u2, cadd(c1,
         cmul(u2, cadd(cs[0],
         cmul(u2, cadd(cs[1],
         cmul(u2, cadd(cs[2],
         cmul(u2, cadd(cs[3],
         cmul(u2, cadd(cs[4],
         cmul(u2, cadd(cs[5],
         cmul(u2, cadd(cs[6],
         cmul(u2, cadd(cs[7],
         cmul(u2, cadd(cs[8],
         cmul(u2, cadd(cs[9],
         cmul(u2, cadd(cs[10],
         cmul(u2, cadd(cs[11],
         cmul(u2, cadd(cs[12],
         cmul(u2, cadd(cs[13],
         cmul(u2, cadd(cs[14],
         cmul(u2, cadd(cs[15],
         cmul(u2, cadd(cs[16],
         cmul(u2, cadd(cs[17],
         cmul(u2, cadd(cs[18],
         cmul(u2, cs[19]))))))))))))))))))))))))))))))))))))))))));
   }

   std::complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);

   if (az <= 1.0L) {
      u = -clog(1.0L - z);
   } else { // az > 1.0L
      auto arg = PI + pz;
      if (arg > PI) { arg -= 2*PI; }
      const auto lmz = std::complex<long double>(lnz, arg); // clog(-z)
      u = -clog(1.0L - 1.0L/z);
      rest = -lmz*(pow2(lmz)/6.0L + zeta2);
   }

   return
      cadd(rest,
      cmul(u, cadd(bf[0],
      cmul(u, cadd(bf[1],
      cmul(u, cadd(bf[2],
      cmul(u, cadd(bf[3],
      cmul(u, cadd(bf[4],
      cmul(u, cadd(bf[5],
      cmul(u, cadd(bf[6],
      cmul(u, cadd(bf[7],
      cmul(u, cadd(bf[8],
      cmul(u, cadd(bf[9],
      cmul(u, cadd(bf[10],
      cmul(u, cadd(bf[11],
      cmul(u, cadd(bf[12],
      cmul(u, cadd(bf[13],
      cmul(u, cadd(bf[14],
      cmul(u, cadd(bf[15],
      cmul(u, cadd(bf[16],
      cmul(u, cadd(bf[17],
      cmul(u, cadd(bf[18],
      cmul(u, cadd(bf[19],
      cmul(u, cadd(bf[20],
      cmul(u, cadd(bf[21],
      cmul(u, cadd(bf[22],
      cmul(u, cadd(bf[23],
      cmul(u, cadd(bf[24],
      cmul(u, cadd(bf[25],
      cmul(u, cadd(bf[26],
      cmul(u, cadd(bf[27],
      cmul(u, cadd(bf[28],
      cmul(u, cadd(bf[29],
      cmul(u, cadd(bf[30],
      cmul(u, cadd(bf[31],
      cmul(u, cadd(bf[32],
      cmul(u, cadd(bf[33],
      cmul(u, cadd(bf[34],
      cmul(u, cadd(bf[35],
      cmul(u, cadd(bf[36],
      cmul(u, cadd(bf[37],
      cmul(u, cadd(bf[38],
      cmul(u, cadd(bf[39],
      cmul(u, cadd(bf[40],
      cmul(u, cadd(bf[41],
      cmul(u, cadd(bf[42],
      cmul(u, cadd(bf[43],
      cmul(u, bf[44]))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))));
}

} // namespace flexiblesusy
