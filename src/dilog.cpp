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

#include "dilog.hpp"
#include "complex.hpp"
#include <cfloat>
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {

   template <typename T, int N>
   T horner(T x, const T (&c)[N]) noexcept
   {
      T p = 0;
      for (int i = N - 1; i >= 0; --i) {
         p = p*x + c[i];
      }
      return p;
   }

   template <int Nstart, typename T, int N>
   Complex<T> horner(const Complex<T>& z, const T (&coeffs)[N]) noexcept
   {
      static_assert(0 <= Nstart && Nstart < N && N >= 2, "invalid array bounds");

      const T r = z.re + z.re;
      const T s = z.re * z.re + z.im * z.im;
      T a = coeffs[N - 1], b = coeffs[N - 2];

      for (int i = N - 3; i >= Nstart; --i) {
         const T t = a;
         a = b + r * a;
         b = coeffs[i] - s * t;
      }

      return Complex<T>(z.re*a + b, z.im*a);
   }

} // anonymous namespace

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
double clausen_2(double x) noexcept
{
   const double PI = 3.141592653589793;
   const std::complex<double> i(0.0, 1.0);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0) {
      x += 2*PI;
   }

   if (std::abs(x) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<double>::epsilon()) {
      return 0.0;
   }

   return std::imag(dilog(std::exp(i*x)));
}

/**
 * @brief Clausen function \f$\mathrm{Cl}_2(\theta) = \mathrm{Im}(\mathrm{Li}_2(e^{i\theta}))\f$ with long double precision
 * @param x real angle
 * @return \f$\mathrm{Cl}_2(\theta)\f$
 */
long double clausen_2(long double x) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const std::complex<long double> i(0.0L, 1.0L);

   while (x >= 2*PI) {
      x -= 2*PI;
   }

   while (x < 0.0L) {
      x += 2*PI;
   }

   if (std::abs(x) < std::numeric_limits<long double>::epsilon() ||
       std::abs(x - PI) < std::numeric_limits<long double>::epsilon() ||
       std::abs(x - 2*PI) < std::numeric_limits<long double>::epsilon()) {
      return 0.0L;
   }

   return std::imag(dilog(std::exp(i*x)));
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(x)\f$
 * @param x real argument
 * @return \f$\mathrm{Li}_2(x)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a
 * maximum error of 4.16e-18.
 */
double dilog(double x) noexcept
{
   const double PI = 3.1415926535897932;
   const double P[] = {
      1.0706105563309304277e+0,
     -4.5353562730201404017e+0,
      7.4819657596286408905e+0,
     -6.0516124315132409155e+0,
      2.4733515209909815443e+0,
     -4.6937565143754629578e-1,
      3.1608910440687221695e-2,
     -2.4630612614645039828e-4
   };
   const double Q[] = {
      1.0000000000000000000e+0,
     -4.5355682121856044935e+0,
      8.1790029773247428573e+0,
     -7.4634190853767468810e+0,
      3.6245392503925290187e+0,
     -8.9936784740041174897e-1,
      9.8554565816757007266e-2,
     -3.2116618742475189569e-3
   };

   double y = 0, r = 0, s = 1;

   // transform to [0, 1/2)
   if (x < -1) {
      const double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5*l - std::log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5*l);
      s = 1;
   } else {
      const double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5*l*l;
      s = -1;
   }

   const double z = y - 0.25;
   const double z2 = z*z;
   const double z4 = z2*z2;

   const double p = P[0] + z * P[1] + z2 * (P[2] + z * P[3]) +
                    z4 * (P[4] + z * P[5] + z2 * (P[6] + z * P[7]));
   const double q = Q[0] + z * Q[1] + z2 * (Q[2] + z * Q[3]) +
                    z4 * (Q[4] + z * Q[5] + z2 * (Q[6] + z * Q[7]));

   return r + s*y*p/q;
}

/**
 * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param x real argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @author Alexander Voigt
 *
 * Implemented as an economized Pade approximation with a maximum
 * error of 2.13e-20 (long double precision) and 1.03e-38 (quadruple
 * precision), respectively.
 */
long double dilog(long double x) noexcept
{
   const long double PI  = 3.14159265358979323846264338327950288L;

#if LDBL_DIG <= 18
   const long double P[] = {
      1.07061055633093042767673531395124630e+0L,
     -5.25056559620492749887983310693176896e+0L,
      1.03934845791141763662532570563508185e+1L,
     -1.06275187429164237285280053453630651e+1L,
      5.95754800847361224707276004888482457e+0L,
     -1.78704147549824083632603474038547305e+0L,
      2.56952343145676978700222949739349644e-1L,
     -1.33237248124034497789318026957526440e-2L,
      7.91217309833196694976662068263629735e-5L
   };
   const long double Q[] = {
      1.00000000000000000000000000000000000e+0L,
     -5.20360694854541370154051736496901638e+0L,
      1.10984640257222420881180161591516337e+1L,
     -1.24997590867514516374467903875677930e+1L,
      7.97919868560471967115958363930214958e+0L,
     -2.87732383715218390800075864637472768e+0L,
      5.49210416881086355164851972523370137e-1L,
     -4.73366369162599860878254400521224717e-2L,
      1.23136575793833628711851523557950417e-3L
   };
#else
   const long double P[] = {
      1.0706105563309304276767353139512463033e+0L,
     -1.0976999936495889694316759019749736514e+1L,
      5.0926847456521671435598196295391041399e+1L,
     -1.4137651316583179225403308056234710093e+2L,
      2.6167438966024905838946321639772104746e+2L,
     -3.4058684964478756354428571563597134889e+2L,
      3.2038320769322335817541819891481680235e+2L,
     -2.2042602382067574460998180123431383994e+2L,
      1.1098617775200190748144566197068215051e+2L,
     -4.0511655788875306581280201485439376229e+1L,
      1.0505482055068989911823839054507173824e+1L,
     -1.8711701614474515075977773430379529915e+0L,
      2.1700009060344649725726270289240541652e-1L,
     -1.5030137964245010798524220531488795883e-2L,
      5.3441650657913964794668981233490297624e-4L,
     -7.0813883289232115572593407126124042503e-6L,
      9.9037749983459843207756017329825486520e-9L
   };
   const long double Q[] = {
      1.0000000000000000000000000000000000000e+0L,
     -1.0552362671556327276432317611336360654e+1L,
      5.0559576537049301822461383847757604795e+1L,
     -1.4554341156550484644439941130107249144e+2L,
      2.8070530313005385434963768448877956414e+2L,
     -3.8295927403204227055447294669165437162e+2L,
      3.8034123327297619837621395664129395552e+2L,
     -2.7878320883603471098333035075019703189e+2L,
      1.5127204204944666467295435798726892255e+2L,
     -6.0401294167088266781532200159750431265e+1L,
      1.7480717870136655299932069333598290551e+1L,
     -3.5731199220918363429317367913920234597e+0L,
      4.9537900060585746351429453828993042241e-1L,
     -4.3750415383481013167382471956794084257e-2L,
      2.2234568413141078158637017098109483650e-3L,
     -5.4149527323164210917629185927908003034e-5L,
      4.1629116823949062782848730828350635143e-7L
   };
#endif

   long double y = 0, r = 0, s = 1;

   // transform to [0, 1/2)
   if (x < -1) {
      const long double l = std::log(1 - x);
      y = 1/(1 - x);
      r = -PI*PI/6 + l*(0.5L*l - std::log(-x));
      s = 1;
   } else if (x == -1) {
      return -PI*PI/12;
   } else if (x < 0) {
      const long double l = std::log1p(-x);
      y = x/(x - 1);
      r = -0.5L*l*l;
      s = -1;
   } else if (x == 0) {
      return 0;
   } else if (x < 0.5L) {
      y = x;
      r = 0;
      s = 1;
   } else if (x < 1) {
      y = 1 - x;
      r = PI*PI/6 - std::log(x)*std::log(1 - x);
      s = -1;
   } else if (x == 1) {
      return PI*PI/6;
   } else if (x < 2) {
      const long double l = std::log(x);
      y = 1 - 1/x;
      r = PI*PI/6 - l*(std::log(1 - 1/x) + 0.5L*l);
      s = 1;
   } else {
      const long double l = std::log(x);
      y = 1/x;
      r = PI*PI/3 - 0.5L*l*l;
      s = -1;
   }

   const long double z = y - 0.25L;

   const long double p = horner(z, P);
   const long double q = horner(z, Q);

   return r + s*y*p/q;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ by Alexander Voigt
 */
std::complex<double> dilog(const std::complex<double>& z_) noexcept
{
   const double PI = 3.1415926535897932;
   const Complex<double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 9}]
   const double bf[10] = {
      - 1.0/4.0,
      + 1.0/36.0,
      - 1.0/3600.0,
      + 1.0/211680.0,
      - 1.0/10886400.0,
      + 1.0/526901760.0,
      - 4.0647616451442255e-11,
      + 8.9216910204564526e-13,
      - 1.9939295860721076e-14,
      + 4.5189800296199182e-16
   };

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return dilog(z.re);
      }
      // z.re > 1
      return { dilog(z.re), -PI*std::log(z.re) };
   }

   const double nz = norm_sqr(z);

   if (nz < std::numeric_limits<double>::epsilon()) {
      return z_;
   }

   Complex<double> u(0.0, 0.0), rest(0.0, 0.0);
   double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5) {
      if (nz > 1) {
         const Complex<double> lz = log(-z);
         u = -log(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      } else { // nz <= 1
         u = -log(1.0 - z);
         rest = 0;
         sgn = 1;
      }
   } else { // z.re > 0.5
      if (nz <= 2*z.re) {
         u = -log(z);
         rest = u*log(1.0 - z) + PI*PI/6;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<double> lz = log(-z);
         u = -log(1.0 - 1.0 / z);
         rest = -0.5*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   const Complex<double> u2(u*u);

   return sgn*(u + u2*(bf[0] + u*horner<1>(u2, bf))) + rest;
}

/**
 * @brief Complex dilogarithm \f$\mathrm{Li}_2(z)\f$ with long double precision
 * @param z_ complex argument
 * @return \f$\mathrm{Li}_2(z)\f$
 * @note Implementation translated from SPheno to C++
 * @author Werner Porod
 * @note translated to C++ and extended to long double precision by Alexander Voigt
 */
std::complex<long double> dilog(const std::complex<long double>& z_) noexcept
{
   const long double PI = 3.14159265358979323846264338327950288L;
   const Complex<long double> z = { std::real(z_), std::imag(z_) };

   // bf[1..N-1] are the even Bernoulli numbers / (2 n + 1)!
   // generated by: Table[BernoulliB[2 n]/(2 n + 1)!, {n, 1, 22}]
   const long double bf[] = {
      -1.0L/4.0L                                 ,
       1.0L/36.0L                                ,
      -1.0L/3600.0L                              ,
       1.0L/211680.0L                            ,
      -1.0L/10886400.0L                          ,
       1.0L/526901760.0L                         ,
      -4.06476164514422552680590938629196667e-11L,
       8.92169102045645255521798731675274885e-13L,
      -1.99392958607210756872364434779378971e-14L,
       4.51898002961991819165047655285559323e-16L,
      -1.03565176121812470144834115422186567e-17L,
#if LDBL_DIG > 18
       2.39521862102618674574028374300098038e-19L,
      -5.58178587432500933628307450562541991e-21L,
       1.30915075541832128581230739918659230e-22L,
      -3.08741980242674029324227976486646243e-24L,
       7.31597565270220342035790560925214859e-26L,
      -1.74084565723400074098905514775970255e-27L,
       4.15763564461389971961789962077522667e-29L,
      -9.96214848828462210319400670245583885e-31L,
       2.39403442489616530052116798789374956e-32L,
      -5.76834735536739008429179316187765424e-34L,
       1.39317947964700797782788660391154833e-35L,
      -3.37212196548508947046847363525493096e-37L
#endif
   };

   // special cases
   if (z.im == 0) {
      if (z.re <= 1) {
         return dilog(z.re);
      }
      // z.re > 1
      return { dilog(z.re), -PI*std::log(z.re) };
   }

   const long double nz = norm_sqr(z);

   if (nz < std::numeric_limits<long double>::epsilon()) {
      return z_;
   }

   Complex<long double> u(0.0L, 0.0L), rest(0.0L, 0.0L);
   long double sgn = 1;

   // transformation to |z|<1, Re(z)<=0.5
   if (z.re <= 0.5L) {
      if (nz > 1) {
         const Complex<long double> lz = log(-z);
         u = -log(1.0L - 1.0L/z);
         rest = -0.5L*lz*lz - PI*PI/6;
         sgn = -1;
      } else { // nz <= 1
         u = -log(1.0L - z);
         rest = 0;
         sgn = 1;
      }
   } else { // z.re > 0.5L
      if (nz <= 2*z.re) {
         u = -log(z);
         rest = u*log(1.0L - z) + PI*PI/6;
         sgn = -1;
      } else { // nz > 2*z.re
         const Complex<long double> lz = log(-z);
         u = -log(1.0L - 1.0L/z);
         rest = -0.5L*lz*lz - PI*PI/6;
         sgn = -1;
      }
   }

   const Complex<long double> u2(u*u);

   return sgn*(u + u2*(bf[0] + u*horner<1>(u2, bf))) + rest;
}

} // namespace flexiblesusy
