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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_wrappers

#include <boost/test/unit_test.hpp>
#include "threshold_loop_functions.hpp"
#include "numerics.h"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_C0 )
{
   using namespace flexiblesusy::threshold_loop_functions;

   struct Values {
      double m1, m2, m3;
   };

   Values value[] = {
      {0, 0, 0},
      {0, 0, 1},
      {0, 1, 1},
      {0, 1, 2},
      {1, 1, 1},
      {1, 1, 0},
      {1, 2, 0},
      {0, 1, 0},
      {1, 0, 1},
      {1, 0, 2},
      {1, 2, 3}
   };

   for (int i = 0; i < sizeof(value)/sizeof(value[0]); i++) {
      BOOST_TEST_MESSAGE("> m1 = " << value[i].m1
                    << ", m2 = " << value[i].m2
                    << ", m3 = " << value[i].m3);
      BOOST_CHECK_CLOSE(softsusy::c0(value[i].m1, value[i].m2, value[i].m3),
                        -Iabc(value[i].m1, value[i].m2, value[i].m3),
                        1e-10);
   }
}

namespace {
   template <typename T> T sqr(T x) { return x*x; }
   template <typename T> T cube(T x) { return x*x*x; }
   template <typename T> T quad(T x) { return x*x*x*x; }
   template <typename T> T pow5(T x) { return x*x*x*x*x; }
   template <typename T> T pow6(T x) { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x) { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x) { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x) { return x*x*x*x*x*x*x*x*x; }
   template <typename T> T pow10(T x) { return x*x*x*x*x*x*x*x*x*x; }
}

double f1_bare(double r)
{
   const double r2 = sqr(r);

   return (6*(r2+3)*r2)/(7*sqr(r2-1))
      + (6*(r2-5)*quad(r)*std::log(r2))/(7*cube(r2-1));
}

double f2_bare(double r)
{
   const double r2 = sqr(r);

   return (2*(r2+11)*r2)/(9*sqr(r2-1))
      + (2*(5*r2-17)*quad(r)*std::log(r2))/(9*cube(r2-1));
}

double f3_bare(double r)
{
   const double r2 = sqr(r);
   const double r4 = quad(r);

   return (2*(r4+9*r2+2))/(3*sqr(r2-1))
      + (2*(r4-7*r2-6)*r2*std::log(r2))/(3*cube(r2-1));
}

double f4_bare(double r)
{
   const double r2 = sqr(r);
   const double r4 = quad(r);

   return (2*(5*r4+25*r2+6))/(7*sqr(r2-1))
      + (2*(r4-19*r2-18)*r2*std::log(r2))/(7*cube(r2-1));
}

double f5_bare(double r1, double r2)
{
   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1+sqr(r1+r2)-r12*r22)/((r12-1)*(r22-1))
      + (cube(r1)*(r12+1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (cube(r2)*(r22+1)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 0.75 * result;
}

double f6_bare(double r1, double r2)
{
   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r12+r22+r1*r2-r12*r22)/((r12-1)*(r22-1))
      + (pow5(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (pow5(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6./7. * result;
}

double f7_bare(double r1, double r2)
{
   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1+r1*r2)/((r12-1)*(r22-1))
      + (cube(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (cube(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6. * result;
}

double f8_bare(double r1, double r2)
{
   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r1+r2)/((r12-1)*(r22-1))
      + (quad(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (quad(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 1.5 * result;
}

BOOST_AUTO_TEST_CASE(test_f1)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
}

BOOST_AUTO_TEST_CASE(test_f2)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
}

BOOST_AUTO_TEST_CASE(test_f3)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
}

BOOST_AUTO_TEST_CASE(test_f4)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
}

BOOST_AUTO_TEST_CASE(test_f5)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x, y = 2.;

   x = 1.1;     BOOST_CHECK_CLOSE(f5(x,y), f5_bare(x,y), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f5(x,y), f5_bare(x,y), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f5(x,y), f5_bare(x,y), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f5(x,y), f5_bare(x,y), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f5(x,y), f5_bare(x,y), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f5(y,x), f5_bare(y,x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f5(y,x), f5_bare(y,x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f5(y,x), f5_bare(y,x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f5(y,x), f5_bare(y,x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f5(y,x), f5_bare(y,x), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 1.02;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 1.01;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 1.001;   BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 1.0001;  BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
}

BOOST_AUTO_TEST_CASE(test_f6)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x, y = 2.;

   x = 1.1;     BOOST_CHECK_CLOSE(f6(x,y), f6_bare(x,y), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f6(x,y), f6_bare(x,y), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f6(x,y), f6_bare(x,y), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f6(x,y), f6_bare(x,y), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f6(x,y), f6_bare(x,y), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f6(y,x), f6_bare(y,x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f6(y,x), f6_bare(y,x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f6(y,x), f6_bare(y,x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f6(y,x), f6_bare(y,x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f6(y,x), f6_bare(y,x), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 1.02;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 1.01;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 1.001;   BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 1.0001;  BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
}

BOOST_AUTO_TEST_CASE(test_f7)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x, y = 2.;

   x = 1.1;     BOOST_CHECK_CLOSE(f7(x,y), f7_bare(x,y), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f7(x,y), f7_bare(x,y), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f7(x,y), f7_bare(x,y), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f7(x,y), f7_bare(x,y), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f7(x,y), f7_bare(x,y), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f7(y,x), f7_bare(y,x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f7(y,x), f7_bare(y,x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f7(y,x), f7_bare(y,x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f7(y,x), f7_bare(y,x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f7(y,x), f7_bare(y,x), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 1.02;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 1.01;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 1.001;   BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 1.0001;  BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 7e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 4e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
}

BOOST_AUTO_TEST_CASE(test_f8)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x, y = 2.;

   x = 1.1;     BOOST_CHECK_CLOSE(f8(x,y), f8_bare(x,y), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f8(x,y), f8_bare(x,y), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f8(x,y), f8_bare(x,y), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f8(x,y), f8_bare(x,y), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f8(x,y), f8_bare(x,y), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f8(y,x), f8_bare(y,x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(f8(y,x), f8_bare(y,x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(f8(y,x), f8_bare(y,x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(f8(y,x), f8_bare(y,x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(f8(y,x), f8_bare(y,x), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 1.02;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 1.01;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 1.001;   BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 1.0001;  BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
}
