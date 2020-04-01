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
#define BOOST_TEST_MODULE test_threshold_loop_functions

#include <boost/test/unit_test.hpp>

#include "config.h"
#include "threshold_loop_functions.hpp"
#include "numerics.h"
#include "dilog.hpp"
#include "logger.hpp"
#include "benchmark.hpp"

#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

const char PATH_SEPARATOR =
#ifdef _WIN32
   '\\';
#else
   '/';
#endif

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

   const auto pass_all_1 = [] (double, double) -> bool { return true; };
   const auto pass_all_2 = [] (double, double, double) -> bool { return true; };

   /// tests mono-variate `func' against values from file
   template <typename T>
   void test_1(const char* func_name, T func, double eps,
               std::function<bool(double, double)> filter = pass_all_1)
   {
      const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR +
                                 func_name + ".txt");
      BOOST_TEST_MESSAGE("reading file " << filename);

      std::ifstream fstr(filename);
      std::string line;
      std::istringstream iss;
      std::vector<double> v(2, 0.0);

      while (std::getline(fstr, line)) {
         iss.clear();
         iss.str(line);

         v.assign(std::istream_iterator<double>(iss),
                  std::istream_iterator<double>());

         if (v.size() < 2) {
            continue;
         }

         const auto x          = v.at(0);
         const auto f_expected = v.at(1);

         if (!filter(x, f_expected)) {
            continue;
         }

         const auto f_fs = func(x);

         BOOST_TEST_MESSAGE("x = " << x << ", " << func_name
                                   << "(expected) = " << f_expected << ", "
                                   << func_name << "(FS) = " << f_fs);

         BOOST_CHECK_CLOSE_FRACTION(f_expected, f_fs, eps);
      }
   }

   /// tests bi-variate `func' against values from file
   template <typename T>
   void test_2(const char* func_name, T func, double eps,
               std::function<bool(double, double, double)> filter = pass_all_2)
   {
      const std::string filename(std::string(TEST_DATA_DIR) + PATH_SEPARATOR +
                                 func_name + ".txt");
      BOOST_TEST_MESSAGE("reading file " << filename);

      std::ifstream fstr(filename);
      std::string line;
      std::istringstream iss;
      std::vector<double> v(2, 0.0);

      while (std::getline(fstr, line)) {
         iss.clear();
         iss.str(line);

         v.assign(std::istream_iterator<double>(iss),
                  std::istream_iterator<double>());

         if (v.size() < 3) {
            continue;
         }

         const auto x          = v.at(0);
         const auto y          = v.at(1);
         const auto f_expected = v.at(2);

         if (!filter(x, y, f_expected)) {
            continue;
         }

         const auto f_fs = func(x, y);

         BOOST_TEST_MESSAGE("x = " << x << ", y = " << y << ", " << func_name
                                   << "(expected) = " << f_expected << ", "
                                   << func_name << "(FS) = " << f_fs);

         BOOST_CHECK_CLOSE_FRACTION(f_expected, f_fs, eps);
      }
   }
} // anonymous namespace

double F1_bare(double x)
{
   const double x2 = sqr(x);

   return x*std::log(x2)/(x2-1);
}

double F2_bare(double x)
{
   const double x2 = sqr(x);

   return 6*x2*(2-2*x2+(1+x2)*std::log(x2))/cube(x2-1);
}

double F3_bare(double x)
{
   const double x2 = sqr(x);

   return 2*x*(5*(1-x2)+(1+4*x2)*std::log(x2))/(3*sqr(x2-1));
}

double F4_bare(double x)
{
   const double x2 = sqr(x);

   return 2*x*(x2-1-std::log(x2))/sqr(x2-1);
}

double F5_bare(double x)
{
   const double x2 = sqr(x);
   const double x4 = quad(x);

   return 3*x*(1-x4+2*x2*std::log(x2))/cube(1-x2);
}

double F6_bare(double x)
{
   const double x2 = sqr(x);

   return (x2-3)/(4*(1-x2)) + x2*(x2-2)/(2*sqr(1.-x2))*std::log(x2);
}

double F7_bare(double x)
{
   const double x2 = sqr(x);
   const double x4 = quad(x);

   return (-3*(x4-6*x2+1.))/(2*sqr(x2-1))
      + (3*x4*(x2-3.))/(cube(x2-1.))*std::log(x2);
}

double F8_bare(double x1, double x2)
{
   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return -2. + 2./(x12-x22)
      *(quad(x1)/(x12-1.)*std::log(x12)
        -quad(x2)/(x22-1.)*std::log(x22));
}

double F9_bare(double x1, double x2)
{
   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return 2./(x12-x22)*(x12/(x12-1.)*std::log(x12)-x22/(x22-1.)*std::log(x22));
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

BOOST_AUTO_TEST_CASE(test_F1)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F1(x), F1_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F1(0)));
   BOOST_CHECK(!std::isnan(F1(1)));
}

BOOST_AUTO_TEST_CASE(test_F1_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("F1", [] (double x) { return F1(x); }, 1e-14);
}

BOOST_AUTO_TEST_CASE(test_F2)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F2(x), F2_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F2(0)));
   BOOST_CHECK(!std::isnan(F2(1)));
}

BOOST_AUTO_TEST_CASE(test_F2_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("F2", [] (double x) { return F2(x); }, 1e-11);
}

BOOST_AUTO_TEST_CASE(test_F3)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F3(x), F3_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F3(0)));
   BOOST_CHECK(!std::isnan(F3(1)));
}

BOOST_AUTO_TEST_CASE(test_F3_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("F3", [] (double x) { return F3(x); }, 1e-13);
}

BOOST_AUTO_TEST_CASE(test_F4)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F4(x), F4_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F4(0)));
   BOOST_CHECK(!std::isnan(F4(1)));
}

BOOST_AUTO_TEST_CASE(test_F4_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("F4", [] (double x) { return F4(x); }, 1e-14);
}

BOOST_AUTO_TEST_CASE(test_F5)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 5e-4);

   x = -1.1;    BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 5e-4);

   x = 0.1;     BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F5(x), F5_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F5(0)));
   BOOST_CHECK(!std::isnan(F5(1)));
   BOOST_CHECK(!std::isnan(F5(-1)));
}

BOOST_AUTO_TEST_CASE(test_F5_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("F5", [] (double x) { return F5(x); }, 1e-11);
}

BOOST_AUTO_TEST_CASE(test_F6)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);

   x = -1.1;    BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F6(x), F6_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F6(0)));
   BOOST_CHECK(!std::isnan(F6(1)));
   BOOST_CHECK(!std::isnan(F6(-1)));
}

BOOST_AUTO_TEST_CASE(test_F6_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   auto filter_small = [](double x, double f) {
      return std::abs(x - 1.0) > 1e-3 && std::abs(f) > 0.01;
   };

   test_1("F6", [] (double x) { return F6(x); }, 1e-14, filter_small);
}

BOOST_AUTO_TEST_CASE(test_F7)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x;

   x = 1.1;     BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);

   x = -1.1;    BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_CLOSE(F7(x), F7_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(F7(0)));
   BOOST_CHECK(!std::isnan(F7(1)));
   BOOST_CHECK(!std::isnan(F7(-1)));
}

BOOST_AUTO_TEST_CASE(test_F7_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("F7", [] (double x) { return F7(x); }, 1e-11);
}

BOOST_AUTO_TEST_CASE(test_F8)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x, y = 2.;

   x = 1.1;     BOOST_CHECK_CLOSE(F8(x,y), F8_bare(x,y), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F8(x,y), F8_bare(x,y), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F8(x,y), F8_bare(x,y), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F8(x,y), F8_bare(x,y), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F8(x,y), F8_bare(x,y), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(F8(y,x), F8_bare(y,x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F8(y,x), F8_bare(y,x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F8(y,x), F8_bare(y,x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F8(y,x), F8_bare(y,x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F8(y,x), F8_bare(y,x), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 2e-2);
   x = 1.02;    BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 2e-2);
   x = 1.01;    BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 2e-2);
   x = 1.001;   BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 2e-2);
   x = 1.0001;  BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 2e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(F8(x,x), F8_bare(x,x + 0.0001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F8(x,1), F8_bare(x,1.000001), 2e-2);
   x = 0.02;    BOOST_CHECK_SMALL(F8(x,1) - F8_bare(x,1.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(F8(x,1) - F8_bare(x,1.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(F8(x,1) - F8_bare(x,1.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(F8(x,1) - F8_bare(x,1.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F8(1,x), F8_bare(1.000001,x), 2e-2);
   x = 0.02;    BOOST_CHECK_SMALL(F8(1,x) - F8_bare(1.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(F8(1,x) - F8_bare(1.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(F8(1,x) - F8_bare(1.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(F8(1,x) - F8_bare(1.00001,x), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F8(x,0), F8_bare(x,0.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(F8(x,0), F8_bare(x,0.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(F8(x,0), F8_bare(x,0.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(F8(x,0), F8_bare(x,0.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(F8(x,0), F8_bare(x,0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F8(0,x), F8_bare(0.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(F8(0,x), F8_bare(0.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(F8(0,x), F8_bare(0.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(F8(0,x), F8_bare(0.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(F8(0,x), F8_bare(0.00001,x), 1e-2);

   x = 0.0001;  BOOST_CHECK_CLOSE(F8(0,0), F8_bare(x,x+0.00001), 1e-2);

   BOOST_CHECK(!std::isnan(F8(0,0)));
   BOOST_CHECK(!std::isnan(F8(0,1)));
   BOOST_CHECK(!std::isnan(F8(1,0)));
   BOOST_CHECK(!std::isnan(F8(1,1)));
   BOOST_CHECK(!std::isnan(F8(2,2)));

   BOOST_CHECK_CLOSE_FRACTION(F8(-1,0), F8(1,0), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(0,-1), F8(0,1), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(-1,-1), F8(1,1), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(-2,0), F8(2,0), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(0,-2), F8(0,2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(-2,1), F8(2,1), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(1,-2), F8(1,2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(-2,2), F8(2,-2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F8(-2,-2), F8(2,2), 1e-15);
}

BOOST_AUTO_TEST_CASE(test_F8_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_2("F8", [] (double x, double y) { return F8(x, y); }, 1e-10);
}

BOOST_AUTO_TEST_CASE(test_F9)
{
   using namespace flexiblesusy::threshold_loop_functions;

   double x, y = 2.;

   x = 1.1;     BOOST_CHECK_CLOSE(F9(x,y), F9_bare(x,y), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F9(x,y), F9_bare(x,y), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F9(x,y), F9_bare(x,y), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F9(x,y), F9_bare(x,y), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F9(x,y), F9_bare(x,y), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(F9(y,x), F9_bare(y,x), 1e-5);
   x = 1.02;    BOOST_CHECK_CLOSE(F9(y,x), F9_bare(y,x), 1e-5);
   x = 1.01;    BOOST_CHECK_CLOSE(F9(y,x), F9_bare(y,x), 1e-5);
   x = 1.001;   BOOST_CHECK_CLOSE(F9(y,x), F9_bare(y,x), 1e-5);
   x = 1.0001;  BOOST_CHECK_CLOSE(F9(y,x), F9_bare(y,x), 1e-5);

   x = 1.1;     BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 1.02;    BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 1.01;    BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 1.001;   BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 1.0001;  BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);

   x = -1.1;    BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = -1.02;   BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = -1.01;   BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = -1.001;  BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = -1.0001; BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x - 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(F9(x,x), F9_bare(x,x + 0.0001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F9(x,1), F9_bare(x,1.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(F9(x,1), F9_bare(x,1.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(F9(x,1), F9_bare(x,1.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(F9(x,1), F9_bare(x,1.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(F9(x,1), F9_bare(x,1.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F9(1,x), F9_bare(1.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(F9(1,x), F9_bare(1.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(F9(1,x), F9_bare(1.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(F9(1,x), F9_bare(1.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(F9(1,x), F9_bare(1.00001,x), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F9(x,0), F9_bare(x,0.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(F9(x,0), F9_bare(x,0.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(F9(x,0), F9_bare(x,0.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(F9(x,0), F9_bare(x,0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(F9(0,x), F9_bare(0.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(F9(0,x), F9_bare(0.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(F9(0,x), F9_bare(0.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(F9(0,x), F9_bare(0.00001,x), 1e-2);

   BOOST_CHECK(!std::isnan(F9(0,1)));
   BOOST_CHECK(!std::isnan(F9(1,0)));
   BOOST_CHECK(!std::isnan(F9(1,1)));
   BOOST_CHECK(!std::isnan(F9(-1,-1)));
   BOOST_CHECK(!std::isnan(F9(2,2)));

   BOOST_CHECK_CLOSE_FRACTION(F9(-1,0), F9(1,0), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(0,-1), F9(0,1), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(-1,-1), F9(1,1), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(-2,0), F9(2,0), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(0,-2), F9(0,2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(-2,1), F9(2,1), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(1,-2), F9(1,2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(-2,2), F9(2,-2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(F9(-2,-2), F9(2,2), 1e-15);
}

BOOST_AUTO_TEST_CASE(test_F9_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_2("F9", [] (double x, double y) { return F9(x, y); }, 1e-11);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(f1(x), f1_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_SMALL(f1(x) - f1_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_SMALL(f1(x) - f1_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(f1(0.)));
}

BOOST_AUTO_TEST_CASE(test_f1_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   auto filter_small = [](double x, double f) {
      return std::abs(f) > 1e-5;
   };

   test_1("f1", [] (double x) { return f1(x); }, 1e-11, filter_small);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(f2(x), f2_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_SMALL(f2(x) - f2_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_SMALL(f2(x) - f2_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_SMALL(f2(x) - f2_bare(x), 1e-5);
   x = 1e-8;    BOOST_CHECK_SMALL(f2(x) - f2_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(f2(0.)));
}

BOOST_AUTO_TEST_CASE(test_f2_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   auto filter_small = [](double x, double f) {
      return std::abs(f) > 1e-5;
   };

   test_1("f2", [] (double x) { return f2(x); }, 1e-11, filter_small);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(f3(x), f3_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_SMALL(f3(x) - f3_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_SMALL(f3(x) - f3_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_SMALL(f3(x) - f3_bare(x), 1e-5);
   x = 1e-8;    BOOST_CHECK_SMALL(f3(x) - f3_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(f3(0.)));
}

BOOST_AUTO_TEST_CASE(test_f3_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("f3", [] (double x) { return f3(x); }, 1e-11);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = -1.02;   BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = -1.01;   BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = -1.001;  BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = -1.0001; BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 0.02;    BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 0.01;    BOOST_CHECK_CLOSE(f4(x), f4_bare(x), 1e-5);
   x = 0.001;   BOOST_CHECK_SMALL(f4(x) - f4_bare(x), 1e-5);
   x = 0.0001;  BOOST_CHECK_SMALL(f4(x) - f4_bare(x), 1e-5);
   x = 0.00001; BOOST_CHECK_SMALL(f4(x) - f4_bare(x), 1e-5);
   x = 1e-8;    BOOST_CHECK_SMALL(f4(x) - f4_bare(x), 1e-5);

   BOOST_CHECK(!std::isnan(f4(0.)));
}

BOOST_AUTO_TEST_CASE(test_f4_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_1("f4", [] (double x) { return f4(x); }, 1e-11);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = -1.02;   BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = -1.01;   BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = -1.001;  BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = -1.0001; BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x - 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f5(x,x), f5_bare(x,x + 0.0001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f5(x,0), f5_bare(x,0.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f5(x,0), f5_bare(x,0.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f5(x,0), f5_bare(x,0.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f5(x,0), f5_bare(x,0.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f5(x,0), f5_bare(x,0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f5(0,x), f5_bare(0.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f5(0,x), f5_bare(0.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f5(0,x), f5_bare(0.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f5(0,x), f5_bare(0.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f5(0,x), f5_bare(0.00001,x), 1e-2);

   x = 0.0001;  BOOST_CHECK_CLOSE(f5(0,0), f5_bare(x,x+0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f5(x,1), f5_bare(x,1.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f5(x,1), f5_bare(x,1.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f5(x,1), f5_bare(x,1.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f5(x,1), f5_bare(x,1.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f5(x,1), f5_bare(x,1.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f5(1,x), f5_bare(1.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f5(1,x), f5_bare(1.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f5(1,x), f5_bare(1.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f5(1,x), f5_bare(1.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f5(1,x), f5_bare(1.00001,x), 1e-2);

   BOOST_CHECK(!std::isnan(f5(0,0)));
   BOOST_CHECK(!std::isnan(f5(0,1)));
   BOOST_CHECK(!std::isnan(f5(1,0)));
   BOOST_CHECK(!std::isnan(f5(1,1)));
   BOOST_CHECK(!std::isnan(f5(2,2)));
}

BOOST_AUTO_TEST_CASE(test_f5_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_2("f5", [] (double x, double y) { return f5(x, y); }, 3e-10);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = -1.02;   BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = -1.01;   BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = -1.001;  BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = -1.0001; BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x - 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f6(x,x), f6_bare(x,x + 0.0001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f6(x,0), f6_bare(x,0.000001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f6(x,0), f6_bare(x,0.000001), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(f6(x,0) - f6_bare(x,0.000001), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(f6(x,0) - f6_bare(x,0.000001), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(f6(x,0) - f6_bare(x,0.000001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f6(0,x), f6_bare(0.000001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f6(0,x), f6_bare(0.000001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(f6(0,x) - f6_bare(0.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(f6(0,x) - f6_bare(0.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(f6(0,x) - f6_bare(0.00001,x), 1e-2);

   x = 0.0001;  BOOST_CHECK_SMALL(f6(0,0) - f6_bare(x,x+0.00001), 1e-5);

   x = 0.1;     BOOST_CHECK_CLOSE(f6(x,1), f6_bare(x,1.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f6(x,1), f6_bare(x,1.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f6(x,1), f6_bare(x,1.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f6(x,1), f6_bare(x,1.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f6(x,1), f6_bare(x,1.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f6(1,x), f6_bare(1.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f6(1,x), f6_bare(1.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f6(1,x), f6_bare(1.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f6(1,x), f6_bare(1.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f6(1,x), f6_bare(1.00001,x), 1e-2);

   BOOST_CHECK(!std::isnan(f6(0,0)));
   BOOST_CHECK(!std::isnan(f6(0,1)));
   BOOST_CHECK(!std::isnan(f6(1,0)));
   BOOST_CHECK(!std::isnan(f6(1,1)));
   BOOST_CHECK(!std::isnan(f6(2,2)));
}

BOOST_AUTO_TEST_CASE(test_f6_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   auto filter_small = [](double x, double y, double f) {
      return std::abs(f) > 1e-9;
   };

   test_2("f6", [] (double x, double y) { return f6(x, y); }, 1e-11, filter_small);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = -1.02;   BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = -1.01;   BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = -1.001;  BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = -1.0001; BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1.2e-1);

   x = 2.1;     BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 4e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f7(x,x), f7_bare(x,x + 0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f7(x,0), f7_bare(x,0.000001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f7(x,0), f7_bare(x,0.000001), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(f7(x,0) - f7_bare(x,0.000001), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(f7(x,0) - f7_bare(x,0.000001), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(f7(x,0) - f7_bare(x,0.000001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f7(0,x), f7_bare(0.000001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f7(0,x), f7_bare(0.000001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(f7(0,x) - f7_bare(0.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(f7(0,x) - f7_bare(0.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(f7(0,x) - f7_bare(0.00001,x), 1e-2);

   x = 0.0001;  BOOST_CHECK_CLOSE(f7(0,0), f7_bare(x,x+0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f7(x,1), f7_bare(x,1.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f7(x,1), f7_bare(x,1.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f7(x,1), f7_bare(x,1.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f7(x,1), f7_bare(x,1.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f7(x,1), f7_bare(x,1.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f7(1,x), f7_bare(1.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f7(1,x), f7_bare(1.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f7(1,x), f7_bare(1.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f7(1,x), f7_bare(1.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f7(1,x), f7_bare(1.00001,x), 1e-2);

   BOOST_CHECK(!std::isnan(f7(0,0)));
   BOOST_CHECK(!std::isnan(f7(0,1)));
   BOOST_CHECK(!std::isnan(f7(1,0)));
   BOOST_CHECK(!std::isnan(f7(1,1)));
   BOOST_CHECK(!std::isnan(f7(2,2)));
}

BOOST_AUTO_TEST_CASE(test_f7_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   test_2("f7", [] (double x, double y) { return f7(x, y); }, 5e-10);
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

   x = -1.1;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = -1.02;   BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = -1.01;   BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = -1.001;  BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = -1.0001; BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x - 0.0001), 1e-2);

   x = 2.1;     BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.02;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.01;    BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.001;   BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);
   x = 2.0001;  BOOST_CHECK_CLOSE(f8(x,x), f8_bare(x,x + 0.0001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f8(x,0), f8_bare(x,0.000001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f8(x,0), f8_bare(x,0.000001), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(f8(x,0) - f8_bare(x,0.000001), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(f8(x,0) - f8_bare(x,0.000001), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(f8(x,0) - f8_bare(x,0.000001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f8(0,x), f8_bare(0.000001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f8(0,x), f8_bare(0.000001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_SMALL(f8(0,x) - f8_bare(0.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_SMALL(f8(0,x) - f8_bare(0.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_SMALL(f8(0,x) - f8_bare(0.00001,x), 1e-2);

   x = 0.0001;  BOOST_CHECK_SMALL(f8(0,0) - f8_bare(x,x+0.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f8(x,1), f8_bare(x,1.00001), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f8(x,1), f8_bare(x,1.00001), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f8(x,1), f8_bare(x,1.00001), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f8(x,1), f8_bare(x,1.00001), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f8(x,1), f8_bare(x,1.00001), 1e-2);

   x = 0.1;     BOOST_CHECK_CLOSE(f8(1,x), f8_bare(1.00001,x), 1e-2);
   x = 0.02;    BOOST_CHECK_CLOSE(f8(1,x), f8_bare(1.00001,x), 1e-2);
   x = 0.01;    BOOST_CHECK_CLOSE(f8(1,x), f8_bare(1.00001,x), 1e-2);
   x = 0.001;   BOOST_CHECK_CLOSE(f8(1,x), f8_bare(1.00001,x), 1e-2);
   x = 0.0001;  BOOST_CHECK_CLOSE(f8(1,x), f8_bare(1.00001,x), 1e-2);

   BOOST_CHECK(!std::isnan(f8(0,0)));
   BOOST_CHECK(!std::isnan(f8(0,1)));
   BOOST_CHECK(!std::isnan(f8(1,0)));
   BOOST_CHECK(!std::isnan(f8(1,1)));
   BOOST_CHECK(!std::isnan(f8(2,2)));
}

BOOST_AUTO_TEST_CASE(test_f8_data)
{
   using namespace flexiblesusy::threshold_loop_functions;

   auto filter_small = [](double x, double y, double f) {
      return std::abs(f) > 1e-10;
   };

   test_2("f8", [] (double x, double y) { return f8(x, y); }, 3e-10, filter_small);
}

namespace {

   std::complex<double> complex_sqrt(double x) {
      return std::sqrt(std::complex<double>(x, 0));
   }

/// \f$phi_{xyz}(x,y,z)\f$ function (Author: Emanuele Bagnaschi)
double phixyz(double x, double y, double z)
{
   const double u = x/z, v = y/z, m = x/y;
   double fac = 0., my_x = 0., my_y = 0., my_z = 0.;
   const double devu = std::fabs(u-1), devv = std::fabs(v-1), devm = std::fabs(m-1);
   const double eps = 0.0000001;
   const double PI = M_PI;

   if (std::abs(sqr(1 - u - v) - 4*u*v) <= std::numeric_limits<double>::epsilon()) {
      return 0.0;
   }

   // The defintion that we implement is valid when x/z < 1 and y/z < 1.
   // We have to reshuffle the arguments to obtain the other branches
   // u > 1 && v < 1 --> x > z && z > y -> y < z < x -> we want x as
   // the third argument x * phi(x,y,z) = z*phi(z,y,x) -> phi(z,y,x) =
   // x/z*phi(x,y,z)

   if (devu >= eps && devv >= eps && devm >= eps) {
      if (u > 1 && v < 1) {
         fac = z/x;
         my_x = z;
         my_y = y;
         my_z = x;
      }
      // u > 1 && v > 1 --> x > z && y > z -> z < y/x -> we want
      // max(x,y) as the third argument
      else if (u > 1 && v > 1) {
         // x as a third argument
         if (x > y) {
            fac = z/x;
            my_x = z;
            my_y = y;
            my_z = x;
         }
         // y as a third argument
         else {
            fac = z/y;
            my_x = z;
            my_y = x;
            my_z = y;
         }
      }
      // u < 1 && v > 1 --> x < z < y --> we want y as the third
      // argument
      else if (u < 1 && v > 1) {
         fac = z/y;
         my_x = z;
         my_y = x;
         my_z = y;
      }
      else if (u < 1 && v < 1) {
         fac = 1.;
         my_x = x;
         my_y = y;
         my_z = z;
      }
      else {
         ERROR("unhandled case in phixyz function!");
      }

      const auto u = my_x/my_z;
      const auto v = my_y/my_z;
      const auto lambda_c = complex_sqrt(sqr(1-u-v)-4*u*v);
      const auto xminus_c = 0.5*(1-(u-v)-lambda_c);
      const auto xplus_c  = 0.5*(1+(u-v)-lambda_c);

      return std::real(fac * 1. / lambda_c *
                       (2. * std::log(xplus_c) * std::log(xminus_c) -
                        std::log(u) * std::log(v) -
                        2. * (dilog(xplus_c) + dilog(xminus_c)) +
                        sqr(PI) / 3.));
   }

   // u and v is equal to one -> all arguments are the same
   if (devu < eps && devv < eps)
      return 2.343907238689459; // from Mathematica

   // u = 1
   if (devu < eps) {
      // if v > 1 then y> z (= x) we have again to reshuffle the
      // argument: here we have phi(z,y,z) with z < y and we want
      // phi(z,z,y)
      if (v > 1) {
         fac = z/y;
         my_x = z;
         // my_y = x; // should be equal to
         my_z = y;

         return std::real(
             fac * 1. / complex_sqrt(1. - 4. * my_x / my_z) *
             (sqr(PI) / 3. +
              2. * sqr(std::log(0.5 -
                                0.5 * complex_sqrt(1. - 4. * my_x / my_z))) -
              sqr(std::log(my_x / my_z)) -
              4. * dilog(0.5 *
                         (1. - complex_sqrt(sqr(1 - 2. * my_x / my_z) -
                                            4. * sqr(my_x) / (sqr(my_z)))))));
      }
      // phi(z,y,z) with z > y, just need to reshufle to phi(y,z,z),
      // i.e. do nothing since they are the same
      else {
        fac = 1.;
        my_x = y;
        // my_y = x;
        my_z = z; // should be equal to my_y

        return std::real(
            fac * 1. /
            (3. * complex_sqrt(my_x * (my_x - 4 * my_z) / (sqr(my_z)))) *
            (sqr(PI) +
             6. * std::log((my_x -
                            my_z * complex_sqrt(my_x * (my_x - 4. * my_z) /
                                                sqr(my_z))) /
                           (2. * my_z)) *
                 std::log(1. - my_x / (2. * my_z) -
                          0.5 * complex_sqrt((sqr(my_x) - 4. * my_x * my_z) /
                                             sqr(my_z))) -
             6. * dilog(0.5 * (2. - complex_sqrt(sqr(my_x) / sqr(my_z) -
                                                 4. * my_x / my_z) -
                               my_x / my_z)) -
             6. * dilog(0.5 * (-complex_sqrt(sqr(my_x) / sqr(my_z) -
                                             4. * my_x / my_z) +
                               my_x / my_z))));
      }
   }

   // v = 1
   if (devv < eps) {
      // if u > 1 we have again to reshuffle the arguments: here we
      // start with phi(x,z,z) with z < x and we want phi(z,z,x)
      if (u > 1) {
         fac = z/x;
         my_x = z;
         // my_y = y;
         my_z = x;

         return std::real(
             fac * 1. / complex_sqrt(1. - 4. * my_x / my_z) *
             (sqr(PI) / 3. +
              2. * sqr(std::log(0.5 -
                                0.5 * complex_sqrt(1. - 4. * my_x / my_z))) -
              sqr(std::log(my_x / my_z)) -
              4. * dilog(0.5 *
                         (1. - complex_sqrt(sqr(1 - 2. * my_x / my_z) -
                                            4. * sqr(my_x) / (sqr(my_z)))))));
      }
      // phi(x,z,z) with z > x, ok as it is
      else {
        fac = 1.;
        my_x = x;
        // my_y = y;
        my_z = z; // should be equal to y

        return std::real(
            fac * 1. /
            (3. * complex_sqrt(my_x * (my_x - 4 * my_z) / (sqr(my_z)))) *
            (sqr(PI) +
             6. * std::log((my_x -
                            my_z * complex_sqrt(my_x * (my_x - 4. * my_z) /
                                                sqr(my_z))) /
                           (2. * my_z)) *
                 std::log(1. - my_x / (2. * my_z) -
                          0.5 * complex_sqrt((sqr(my_x) - 4. * my_x * my_z) /
                                             sqr(my_z))) -
             6. * dilog(0.5 * (2. - complex_sqrt(sqr(my_x) / sqr(my_z) -
                                                 4. * my_x / my_z) -
                               my_x / my_z)) -
             6. * dilog(0.5 * (-complex_sqrt(sqr(my_x) / sqr(my_z) -
                                             4. * my_x / my_z) +
                               my_x / my_z))));
      }
   }

   if (devm < eps) {
      // if (v = ) u > 1, we are in the case phi(x,x,z) with x > z. We
      // have to reshufle to phi(z,x,x)
      if (u > 1) {
         fac = z/x;
         my_x = z;
         // my_y = y;
         my_z = x;

         return std::real(
             fac * 1. /
             (3. * complex_sqrt(my_x * (my_x - 4 * my_z) / (sqr(my_z)))) *
             (sqr(PI) +
              6. * std::log((my_x -
                             my_z * complex_sqrt(my_x * (my_x - 4. * my_z) /
                                                 sqr(my_z))) /
                            (2. * my_z)) *
                  std::log(1. - my_x / (2. * my_z) -
                           0.5 * complex_sqrt((sqr(my_x) - 4. * my_x * my_z) /
                                              sqr(my_z))) -
              6. * dilog(0.5 * (2. - complex_sqrt(sqr(my_x) / sqr(my_z) -
                                                  4. * my_x / my_z) -
                                my_x / my_z)) -
              6. * dilog(0.5 * (-complex_sqrt(sqr(my_x) / sqr(my_z) -
                                              4. * my_x / my_z) +
                                my_x / my_z))));
      }
      // if (v = u) < 1 we can directly use phi(x,x,z)
      else {
        fac = 1.;
        my_x = x;
        // my_y = y;
        my_z = z;

        return std::real(
            fac * 1. / complex_sqrt(1. - 4. * my_x / my_z) *
            (sqr(PI) / 3. +
             2. * sqr(std::log(0.5 -
                               0.5 * complex_sqrt(1. - 4. * my_x / my_z))) -
             sqr(std::log(my_x / my_z)) -
             4. * dilog(0.5 *
                        (1. - complex_sqrt(sqr(1 - 2. * my_x / my_z) -
                                           4. * sqr(my_x) / (sqr(my_z)))))));
      }
   }

   FATAL("unhandled case in phixyz function!");

   return std::numeric_limits<double>::quiet_NaN();
}
} // anonymous namespace

struct XYZ {
   double x{}, y{}, z{};
};

BOOST_AUTO_TEST_CASE(test_phixyz)
{
   using namespace flexiblesusy::threshold_loop_functions;
   double x = 0., y = 0., z = 0.;

   const XYZ xyz[] = {
      {1, 2, 3},
      {2, 3, 1},
      {3, 1, 2},
      {2, 1, 3},
      {1, 3, 2},
      {3, 2, 1},

      {1, 1, 2},
      {1, 2, 1},
      {2, 1, 1},

      {1, 1, 10},
      {1, 10, 1},
      {10, 1, 1},

      {1, 10, 20},
      {10, 20, 1},
      {20, 1, 10},

      {1, 2, 2},
      {2, 2, 1},
      {2, 1, 2},

      {862.132647151542, 862.132267190459, 684.729637476883},

      // lambda = 0
      {1, 0, 1},
      {4, 1, 1},

      // lambda ~ 0
      // {4, 1 + 1e-5, 1},
      // {200.220790830763, 599.56612604427, 106.834963457636},

      {1, 1, 1}
   };

   for (const auto s: xyz) {
      const auto x = s.x, y = s.y, z = s.z;

      BOOST_TEST_MESSAGE("x = " << x << ", y = " << y << ", z = " << z);

      const double prec = 1e-15;

      // test identities
      BOOST_CHECK_CLOSE_FRACTION(phixyz(x,y,z)  , phixyz(y,x,z) , prec);
      BOOST_CHECK_CLOSE_FRACTION(phi_xyz(x,y,z) , phi_xyz(y,x,z), prec);

      // test identities
      BOOST_CHECK_CLOSE_FRACTION(x*phixyz(x,y,z) , z*phixyz(z,y,x) , prec);
      BOOST_CHECK_CLOSE_FRACTION(x*phi_xyz(x,y,z), z*phi_xyz(z,y,x), prec);

      // compare implementations
      BOOST_CHECK_CLOSE_FRACTION(phixyz(x,y,z), phi_xyz(x,y,z), prec);
   }
}

#ifdef ENABLE_RANDOM

template <class T>
std::vector<XYZ> generate_random_triples(
   unsigned n, T start, T stop)
{
   const auto x = generate_random_data<T>(n, start, stop);
   const auto y = generate_random_data<T>(n, start, stop);
   const auto z = generate_random_data<T>(n, start, stop);

   std::vector<XYZ> v(5*n);

   for (int i = 0; i < n; i++) {
      v[i] = {x[i], y[i], z[i]};
   }

   for (int i = 0; i < n; i++) {
      v[n + i] = {x[i], x[i], y[i]};
   }

   for (int i = 0; i < n; i++) {
      v[2*n + i] = {x[i], y[i], x[i]};
   }

   for (int i = 0; i < n; i++) {
      v[3*n + i] = {y[i], x[i], x[i]};
   }

   for (int i = 0; i < n; i++) {
      v[4*n + i] = {x[i], x[i], x[i]};
   }

   return v;
}

BOOST_AUTO_TEST_CASE(test_phi_random)
{
   const unsigned N = 10000;
   const auto triples = generate_random_triples(N, 1.0, 1000.0);

   auto phi_fs = [](const XYZ& t) {
      return flexiblesusy::threshold_loop_functions::phi_xyz(t.x, t.y, t.z);
   };

   auto phi_eb = [](const XYZ& t) {
      return phixyz(t.x, t.y, t.z);
   };

   // low testing precision since phi_xyz becomes unstable when lambda ~ 0
   const double prec = 1e-7;

   for (const auto t: triples) {
      BOOST_CHECK_CLOSE_FRACTION(phi_fs(t), phi_eb(t), prec);
   }
}

BOOST_AUTO_TEST_CASE(bench_phi)
{
   const unsigned N = 10000;
   const auto triples = generate_random_triples(N, 1.0, 1000.0);

   auto phi_fs = [](const XYZ& t) {
      return flexiblesusy::threshold_loop_functions::phi_xyz(t.x, t.y, t.z);
   };

   auto phi_eb = [](const XYZ& t) {
      return phixyz(t.x, t.y, t.z);
   };

   const auto time_phi_fs_in_s = time_in_seconds(phi_fs, triples)/N;
   const auto time_phi_eb_in_s = time_in_seconds(phi_eb, triples)/N;

   BOOST_TEST_MESSAGE("average run-time for phi_xyz [FS]: " << time_phi_fs_in_s << " s");
   BOOST_TEST_MESSAGE("average run-time for phi_xyz [EB]: " << time_phi_eb_in_s << " s");
}

#endif
