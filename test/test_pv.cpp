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

#include <limits>
#include <cmath>
#include "pv.hpp"
#include "wrappers.hpp"
#include "stopwatch.hpp"
#include "numerics.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pv

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

using namespace std;
using namespace flexiblesusy;
using namespace flexiblesusy::passarino_veltman;

BOOST_AUTO_TEST_CASE(test_real_parts)
{
    BOOST_CHECK_CLOSE_FRACTION(ReA0 (2,     1),  0.61370563888010943, 1e-14);
    BOOST_CHECK_CLOSE_FRACTION(ReB0 (2,3,4, 1), -1.14772143349266490, 1e-14);
    BOOST_CHECK_CLOSE_FRACTION(ReB1 (2,3,4, 1),  0.59926550299197479, 1e-14);
    BOOST_CHECK_CLOSE_FRACTION(ReB00(2,3,4, 1), -0.24981786818504056, 1e-14);
}

const double scale  = 100;
const double scale2 = Sqr(scale);
const double p  = 91.0;
const double p2 = Sqr(p);

BOOST_AUTO_TEST_CASE( test_ReA0 )
{
   BOOST_CHECK_EQUAL(ReA0(0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReB0 )
{
   BOOST_CHECK_EQUAL(ReB0(0., 0., 0., scale2), 0.);

   BOOST_CHECK_CLOSE(ReB0(p2, 0., 0., scale2), 2. - log(p2/scale2), 1.0e-10);
   BOOST_CHECK_CLOSE(ReB0(0., p2, 0., scale2), 1. - log(p2/scale2), 1.0e-10);
   BOOST_CHECK_EQUAL(ReB0(0., 0., p2, scale2), ReB0(0., p2, 0., scale2));

   BOOST_CHECK_CLOSE(ReB0(p2, p2, 0., scale2), 2. - log(p2/scale2), 0.005);
   BOOST_CHECK_CLOSE(ReB0(0., p2, p2, scale2), 0. - log(p2/scale2), 0.005);
   BOOST_CHECK_EQUAL(ReB0(p2, 0., p2, scale2), ReB0(p2, p2, 0., scale2));
}

BOOST_AUTO_TEST_CASE( test_ReB1 )
{
   BOOST_CHECK_EQUAL(ReB1(0., 0., 0., scale2), 0.);

   BOOST_CHECK_CLOSE(ReB1(p2, 0,0, scale2), -0.5*ReB0(p2, 0,0, scale2), 1e-10);
   BOOST_CHECK_CLOSE(ReB1(0,p2,p2, scale2), -0.5*ReB0(0,p2,p2, scale2), 1e-10);
}

BOOST_AUTO_TEST_CASE( test_ReB00 )
{
   BOOST_CHECK_EQUAL(ReB00(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReB22 )
{
   BOOST_CHECK_EQUAL(ReB22(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReH0 )
{
   BOOST_CHECK_EQUAL(ReH0(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReF0 )
{
   BOOST_CHECK_EQUAL(ReF0(0., 0., 0., scale2), 0.);
}

BOOST_AUTO_TEST_CASE( test_ReG0 )
{
   BOOST_CHECK_EQUAL(ReG0(0., 0., 0., scale2), 0.);
}

struct B0_args {
   double p, m1, m2, q, prec, speedup;
};

double Softsusy_ReB0(double p2, double m2a, double m2b, double scl2)
{
   return b0(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
}

double Fast_ReB0(double p2, double m2a, double m2b, double scl2)
{
   return b0_fast(sqrt(p2), sqrt(m2a), sqrt(m2b), sqrt(scl2));
}

BOOST_AUTO_TEST_CASE(test_b0_b0_fast)
{
   Stopwatch stopwatch;

   B0_args args[] = {
      {91., 0.  , 0.  , 91., 1.0e-12, 1.0},
      {91., 0.  , 91. , 91., 4.0e-06, 5.0},
      {91., 91. , 91. , 91., 1.0e-12, 5.0},
      {91., 91. , 0.  , 91., 4.0e-06, 5.0},
      {91., 0.  , 100., 91., 1.0e-12, 2.0},
      {91., 100., 100., 91., 1.0e-12, 1.0}
   };

   for (unsigned i = 0; i < sizeof(args)/sizeof(args[0]); i++) {
      const double p2 = Sqr(args[i].p);
      const double m12 = Sqr(args[i].m1);
      const double m22 = Sqr(args[i].m2);
      const double q2 = Sqr(args[i].q);
      const double prec = args[i].prec;
      const double speedup = args[i].speedup;
      volatile double ss, fs;

      BOOST_MESSAGE("testing B0(" << p2 << "," << m12 << "," << m22 << "," << q2 << ")");

      stopwatch.start();
      for (unsigned i = 0; i < 100000; i++)
         ss = Softsusy_ReB0(p2, m12, m22, scale2);
      stopwatch.stop();
      const double ss_time = stopwatch.get_time_in_seconds();

      stopwatch.start();
      for (unsigned i = 0; i < 100000; i++)
         fs = Fast_ReB0(p2, m12, m22, scale2);
      stopwatch.stop();
      const double fs_time = stopwatch.get_time_in_seconds();

      BOOST_CHECK_CLOSE_FRACTION(ss, fs, prec);

      if (speedup > 1.0)
         BOOST_CHECK_GT(ss_time, speedup*fs_time);
   }
}

#if defined(ENABLE_LOOPTOOLS) || defined(ENABLE_FFLITE)

const double eps = numeric_limits<double>::min();
const complex<double> neg_m2 = complex<double>(-0.1, -eps);

BOOST_AUTO_TEST_CASE(test_negative_m2)
{
    BOOST_CHECK_SMALL(abs(A0 (neg_m2,1) -
	complex<double>(-0.3302585092994046,-0.3141592653589793)), 1e-14);

    BOOST_CHECK_SMALL(abs(B0 (2,neg_m2,4, 1) -
	complex<double>(0.04191013584834996,0.14336849876951852)), 1e-14);
    BOOST_CHECK_SMALL(abs(B0 (2,3,neg_m2, 1) -
	complex<double>(0.49731133798115157,0.24955609428693504)), 1e-14);

    BOOST_CHECK_SMALL(abs(B1 (2,neg_m2,4, 1) -
	complex<double>(0.32573255511542315,-0.003271354485747599)), 1e-14);
    BOOST_CHECK_SMALL(abs(B1 (2,3,neg_m2, 1) -
	complex<double>(-0.6254665451021995,-0.2396442038760973)), 1e-14);

    BOOST_CHECK_SMALL(abs(B00(2,neg_m2,4, 1) -
	complex<double>(0.16595591599028536,-0.003633975888972291)), 1e-14);
    BOOST_CHECK_SMALL(abs(B00(2,3,neg_m2, 1) -
	complex<double>(0.2828439119832702,-0.006501356567577554)), 1e-14);
}

#endif
