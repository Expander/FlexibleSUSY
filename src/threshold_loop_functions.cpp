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

#include "threshold_loop_functions.hpp"

#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace threshold_loop_functions {

namespace {
   const double EPS = 1.e-8;
   double sqr(double x) { return x*x; }

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
}

double F1(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return x*std::log(x2)/(x2-1);
}

double F2(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return 6*x2*(2-2*x2+(1+x2)*std::log(x2))/std::pow(x2-1,3);
}

double F3(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return 2*x*(5*(1-x2)+(1+4*x2)*std::log(x2))/(3*sqr(x2-1));
}

double F4(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);

   return 2*x*(x2-1-std::log(x2))/sqr(x2-1);
}

double F5(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);
   const double x4 = std::pow(x,4);

   return 3*x*(1-x4+2*x2*std::log(x2))/std::pow(1-x2,3);
}

double F6(double x)
{
   if (is_equal(x, 1., EPS))
      return 0.;

   const double x2 = sqr(x);

   return (x2-3)/(4*(1-x2)) + x2*(x2-2)/(2*sqr(1.-x2))*std::log(x2);
}

double F7(double x)
{
   if (is_equal(x, 1., EPS))
      return 1.;

   const double x2 = sqr(x);
   const double x4 = std::pow(x,4);

   return (-3*(x4-6*x2+1.))/(2*sqr(x2-1))
      + (3*x4*(x2-3.))/(std::pow(x2-1.,3))*std::log(x2);
}

double F8(double x1, double x2)
{
   if (is_equal(x1, 1., EPS) && is_equal(x2, 1., EPS))
      return 1.;

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return -2. + 2./(x12-x22)
      *(std::pow(x1,4)/(x12-1.)*std::log(x12)
        -std::pow(x2,4)/(x22-1.)*std::log(x22));
}

double F9(double x1, double x2)
{
   if (is_equal(x1, 1., EPS) && is_equal(x2, 1., EPS))
      return 1.;

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return 2./(x12-x22)*(x12/(x12-1.)*std::log(x12)-x22/(x22-1.)*std::log(x22));
}

double f(double r)
{
   return F5(r);
}

double g(double r)
{
   return F7(r);
}

double f1(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);

   return (6*(r2+3)*r2)/(7*sqr(r2-1))
      + (6*(r2-5)*std::pow(r,4)*std::log(r2))/(7*std::pow(r2-1,3));
}

double f2(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);

   return (2*(r2+11)*r2)/(9*sqr(r2-1))
      + (2*(5*r2-17)*std::pow(r,4)*std::log(r2))/(9*std::pow(r2-1,3));
}

double f3(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);
   const double r4 = std::pow(r,4);

   return (2*(r4+9*r2+2))/(3*sqr(r2-1))
      + (2*(r4-7*r2-6)*r2*std::log(r2))/(3*std::pow(r2-1,3));
}

double f4(double r)
{
   if (is_equal(r, 1., EPS))
      return 1.;

   const double r2 = sqr(r);
   const double r4 = std::pow(r,4);

   return (2*(5*r4+25*r2+6))/(7*sqr(r2-1))
      + (2*(r4-19*r2-18)*r2*std::log(r2))/(7*std::pow(r2-1,3));
}

/// f5(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f5_1_1(double r1, double r2)
{
   return 0.772943722943723
      - 0.5524891774891774*r2
      + 0.7870670995670994*sqr(r2)
      - 0.3316558441558441*std::pow(r2,3)
      + 0.056277056277056266*std::pow(r2,4)
      + r1*(-0.5524891774891774
            + 1.0700757575757573*r2
            - 0.6625541125541123*sqr(r2)
            + 0.22483766233766228*std::pow(r2,3)
            - 0.03344155844155843*std::pow(r2,4))
      + std::pow(r1,3)*(-0.33165584415584404
                        + 0.22483766233766223*r2
                        - 0.08755411255411245*sqr(r2)
                        + 0.01650432900432896*std::pow(r2,3)
                        - 0.0007034632034631958*std::pow(r2,4))
      + std::pow(r1,4)*(0.05627705627705626
                        - 0.03344155844155841*r2
                        + 0.010281385281385256*sqr(r2)
                        - 0.0007034632034631921*std::pow(r2,3)
                        - 0.0002705627705627725*std::pow(r2,4))
      + sqr(r1)*(0.7870670995670994 - 0.6625541125541123*r2
                 + 0.32061688311688297*sqr(r2)
                 - 0.08755411255411248*std::pow(r2,3)
                 + 0.01028138528138527*std::pow(r2,4));
}

/// f5(r1,r2) in the limit r1 -> 1
static double f5_1_r2(double r1, double r2)
{
   return (-0.025*std::pow(-1. + r1,3)*(
              4. - 17.*r2 + 4.*sqr(r2)
              - 25.*std::pow(r2,3)
              - 20.*std::pow(r2,4)
              + 41.*std::pow(r2,5)
              + 12.*std::pow(r2,6)
              + std::pow(r2,7)
              + (-30.*std::pow(r2,3) - 30.*std::pow(r2,5))
              *std::log(sqr(r2))))/(std::pow(-1. + r2,6)*std::pow(1. + r2,2))
      - (0.125*std::pow(-1. + r1,2)*(
            1. - 4.*r2 + sqr(r2)
            - 4.*std::pow(r2,3)
            - 5.*std::pow(r2,4)
            + 8.*std::pow(r2,5)
            + 3.*std::pow(r2,6)
            + (-6.*std::pow(r2,3) - 6.*std::pow(r2,5))
            *std::log(sqr(r2))))/(std::pow(-1. + r2,5)*std::pow(1. + r2,2))
      + (0.75*(-1 + r2 + 2*sqr(r2)
               - std::pow(r2,4)
               - std::pow(r2,5)
               + (std::pow(r2,3) + std::pow(r2,5))
               *std::log(sqr(r2))))/(std::pow(-1 + r2,3)*std::pow(1 + r2,2))
      + (0.25*(-1. + r1)*(
            1. - 1.*r2 - 2.*sqr(r2) + 8.*std::pow(r2,3)
            + std::pow(r2,4) - 7.*std::pow(r2,5)
            + (3.*std::pow(r2,3) + 3.*std::pow(r2,5))
            *std::log(sqr(r2))))/(std::pow(-1. + r2,4)*std::pow(1. + r2,2))
      + (0.05*std::pow(-1. + r1,4)*(
            -1. + 4.5*r2 + 2.*sqr(r2)
            + 16.5*std::pow(r2,3) - 16.5*std::pow(r2,5) - 2.*std::pow(r2,6)
            - 4.5*std::pow(r2,7) + 1.*std::pow(r2,8)
            + (15.*std::pow(r2,3) + 15.*std::pow(r2,5))
            *std::log(sqr(r2))))/(std::pow(-1. + r2,7)*std::pow(1. + r2,2));
}

double f5(double r1, double r2)
{
   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f5_1_1(r1, r2);

   if (is_equal(r1, 1., 0.01))
      return f5_1_r2(r1, r2);

   if (is_equal(r2, 1., 0.01))
      return f5_1_r2(r2, r1);

   const double r12 = sqr(r1);

   const double result
      = (1+sqr(r1+r2)-r12*sqr(r2))/((r12-1)*(sqr(r2)-1))
      + (std::pow(r1,3)*(r12+1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,3)*(sqr(r2)+1)*std::log(sqr(r2)))/((r1-r2)*sqr(sqr(r2)-1));

   return 0.75 * result;
}

/// f6(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f6_1_1(double r1, double r2)
{
   return std::pow(r1,3)*(-0.22455163883735313
                          + 0.24582560296846012*r2
                          - 0.15819418676561536*sqr(r2)
                          + 0.058750773036487326*std::pow(r2,3)
                          - 0.0095856524427953*std::pow(r2,4))
      + std::pow(r1,4)*(0.04236239950525665
                        - 0.04223871366728509*r2
                        + 0.026283240568954855*sqr(r2)
                        - 0.009585652442795299*std::pow(r2,3)
                        + 0.001546072974644403*std::pow(r2,4))
      + sqr(r1)*(0.41403834260977124
                 - 0.5990105132962276*r2
                 + 0.41280148423005564*sqr(r2)
                 - 0.15819418676561534*std::pow(r2,3)
                 + 0.026283240568954855*std::pow(r2,4))
      + (-0.08756957328385832
         + 0.6100803957946813*r2
         - 1.1505875077303678*sqr(r2)
         + 0.02838589981447649*std::pow(r2,3)
         + 2.29802102659245*std::pow(r2,4)
         - 2.9131106988249806*std::pow(r2,5)
         + 1.566419294990722*std::pow(r2,6)
         - 0.39400123685837934*std::pow(r2,7)
         + 0.042362399505256616*std::pow(r2,8))/std::pow(-1. + r2,4)
      - (0.04223871366728506*r1*(-6.150805270863825
                                 + 7.9121522693997095*r2
                                 + 44.04099560761344*sqr(r2)
                                 - 138.08931185944363*std::pow(r2,3)
                                 + 169.98243045387994*std::pow(r2,4)
                                 - 112.3367496339678*std::pow(r2,5)
                                 + 43.46120058565154*std::pow(r2,6)
                                 - 9.8199121522694*std::pow(r2,7)
                                 + 1.*std::pow(r2,8)))/std::pow(-1. + r2,4);
}

/// f6(r1,r2) in the limit r1 -> 1
static double f6_1_r2(double r1, double r2)
{
   return (0.8571428571428571*(-1. + r1)*(
              0.6666666666666666 - 1.1666666666666667*r2
              - 1.3333333333333333*sqr(r2) + 3.333333333333333*std::pow(r2,3)
              + 0.6666666666666666*std::pow(r2,4)
              - 2.1666666666666665*std::pow(r2,5)
              + 1.*std::pow(r2,5)*std::log(sqr(r2))))
      /(std::pow(-1. + r2,4)*std::pow(1. + r2,2))
      + (0.42857142857142855*(-1 + 2*sqr(r2) + 2*std::pow(r2,3)
                              - std::pow(r2,4) - 2*std::pow(r2,5)
                              + 2*std::pow(r2,5)*std::log(sqr(r2))))
      /(std::pow(-1 + r2,3)*std::pow(1 + r2,2))
      + (0.14285714285714285*std::pow(-1 + r1,2)*(
            r2 - 4*sqr(r2) + 4*std::pow(r2,3) + 8*std::pow(r2,4)
            - 5*std::pow(r2,5) - 4*std::pow(r2,6)
            + 6*std::pow(r2,5)*std::log(sqr(r2))))
      /(std::pow(-1 + r2,5)*std::pow(1 + r2,2))
      + (0.05714285714285714*std::pow(-1. + r1,3)*(
            -1. + 5.5*r2 - 11.*sqr(r2) + 5.*std::pow(r2,3)
            + 25.*std::pow(r2,4) - 11.5*std::pow(r2,5) - 13.*std::pow(r2,6)
            + 1.*std::pow(r2,7) + 15.*std::pow(r2,5)*std::log(sqr(r2))))
      /(std::pow(-1. + r2,6)*std::pow(1. + r2,2))
      + (0.014285714285714285*std::pow(-1. + r1,4)*(
            -3. + 18.*r2 - 40.*sqr(r2) + 24.*std::pow(r2,3)
            + 90.*std::pow(r2,4) - 42.*std::pow(r2,5) - 48.*std::pow(r2,6)
            + std::pow(r2,8) + 60.*std::pow(r2,5)*std::log(sqr(r2))))
      /(std::pow(-1. + r2,7)*std::pow(1. + r2,2));
}

double f6(double r1, double r2)
{
   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f6_1_1(r1, r2);

   if (is_equal(r1, 1., 0.01))
      return f6_1_r2(r1, r2);

   if (is_equal(r2, 1., 0.01))
      return f6_1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r12+r22+r1*r2-r12*r22)/((r12-1)*(r22-1))
      + (std::pow(r1,5)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,5)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6./7. * result;
}

/// f7(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f7_1_1(double r1, double r2)
{
   return (15700 - 14411*r2 + 7850*sqr(r2) - 2498*std::pow(r2,3)
           + 355*std::pow(r2,4)
           + sqr(r1)*(7850 - 2558*r2 - 750*sqr(r2) + 940*std::pow(r2,3)
                      - 235*std::pow(r2,4))
           + std::pow(r1,4)*(355 + 65*r2 - 235*sqr(r2) + 142*std::pow(r2,3)
                             - 30*std::pow(r2,4))
           + r1*(-14411 + 8375*r2 - 2558*sqr(r2) + 180*std::pow(r2,3)
                 + 65*std::pow(r2,4))
           + std::pow(r1,3)*(-2498 + 180*r2 + 940*sqr(r2)
                             - 645*std::pow(r2,3) + 142*std::pow(r2,4)))
      /2310.;
}

/// f7(r1,r2) in the limit r1 -> 1
static double f7_1_r2(double r1, double r2)
{
   return (-10*(-1 + r1)*std::pow(-1 + r2,3)*(
              2 - 5*r2 - 4*sqr(r2) + 4*std::pow(r2,3) + 2*std::pow(r2,4)
              + std::pow(r2,5) - 6*std::pow(r2,3)*std::log(sqr(r2)))
           - 30*std::pow(-1 + r2,4)*(
              1 - 2*r2 - 2*sqr(r2) + 2*std::pow(r2,3) + std::pow(r2,4)
              - 2*std::pow(r2,3)*std::log(sqr(r2)))
           + 10*std::pow(-1 + r1,2)*std::pow(-1 + r2,2)*(
              -1 + 3*r2 + 3*sqr(r2) - 3*std::pow(r2,4) - 3*std::pow(r2,5)
              + std::pow(r2,6) + 6*std::pow(r2,3)*std::log(sqr(r2)))
           + 2*std::pow(-1 + r1,3)*(-1 + r2)*(
              -2 + 6*r2 + 18*sqr(r2) + 15*std::pow(r2,3) - 30*std::pow(r2,4)
              - 18*std::pow(r2,5) + 14*std::pow(r2,6) - 3*std::pow(r2,7)
              + 30*std::pow(r2,3)*std::log(sqr(r2)))
           + std::pow(-1 + r1,4)*(
              -1 + 48*sqr(r2) + 42*std::pow(r2,3) - 90*std::pow(r2,4)
              - 24*std::pow(r2,5) + 40*std::pow(r2,6) - 18*std::pow(r2,7)
              + 3*std::pow(r2,8) + 60*std::pow(r2,3)*std::log(sqr(r2))))
      /(10.*std::pow(-1 + r2,7)*std::pow(1 + r2,2));
}

double f7(double r1, double r2)
{
   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f7_1_1(r1, r2);

   if (is_equal(r1, 1., 0.01))
      return f7_1_r2(r1, r2);

   if (is_equal(r2, 1., 0.01))
      return f7_1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (1+r1*r2)/((r12-1)*(r22-1))
      + (std::pow(r1,3)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,3)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6. * result;
}

/// f8(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f8_1_1(double r1, double r2)
{
   return 1. - 0.1*std::pow(-1. + r1,2)
      + 0.07500000000000022*std::pow(-1. + r1,3)
      - 0.04285714285714286*std::pow(-1. + r1,4)
      + (0.23928571428571443
         - 0.4642857142857147*r1
         + 0.3321428571428576*sqr(r1)
         - 0.12857142857142873*std::pow(r1,3)
         + 0.02142857142857143*std::pow(r1,4))*(-1 + r2)
      - (5.551115123125783e-18*std::pow(-1. + r1,3))/std::pow(-1. + r2,6)
      + (5.551115123125783e-18*std::pow(-1. + r1,3))/std::pow(-1. + r2,5)
      - (4.163336342344338e-18*std::pow(-1. + r1,3))/std::pow(-1. + r2,4)
      + (2.7755575615628915e-18*std::pow(-1. + r1,3))/std::pow(-1. + r2,3)
      + (1.7590096046404825e-16*std::pow(-1. + r1,3))/std::pow(-1. + r2,2)
      - (2.654126918244515e-16*std::pow(-1. + r1,3))/(-1. + r2)
      - 0.009523809523809525*(26.125000000000007 - 27.62500000000002*r1
                              + 17.250000000000025*sqr(r1)
                              - 6.250000000000009*std::pow(r1,3)
                              + 1.*std::pow(r1,4))*std::pow(-1. + r2,2)
      + 0.0035714285714285713*(42.666666666666686 - 36.00000000000004*r1
                               + 20.000000000000046*sqr(r1)
                               - 6.666666666666682*std::pow(r1,3)
                               + 1.*std::pow(r1,4))*std::pow(-1. + r2,3)
      - 0.0008658008658008658*(90.37500000000001 - 63.1250000000001*r1
                               + 29.37500000000011*sqr(r1)
                               - 8.125000000000037*std::pow(r1,3)
                               + 1.*std::pow(r1,4))*std::pow(-1. + r2,4);
}

/// f8(r1,r2) in the limit r1 -> 1
static double f8_1_r2(double r1, double r2)
{
   return (-0.07500000000000001*std::pow(-1. + r1,4)*(
              -1. + 5.333333333333333*r2 - 8.*sqr(r2) - 16.*std::pow(r2,3)
              + 16.*std::pow(r2,5) + 8.*std::pow(r2,6)
              - 5.333333333333333*std::pow(r2,7) + 1.*std::pow(r2,8)
              - 20.*std::pow(r2,4)*std::log(sqr(r2))))
      /(std::pow(-1. + r2,7)*std::pow(1. + r2,2))
      - (1.*(-1. + r1)*(
            -0.25 + 1.*r2 - 1.*sqr(r2) - 2.*std::pow(r2,3)
            + 1.25*std::pow(r2,4) + 1.*std::pow(r2,5)
            - 1.5*std::pow(r2,4)*std::log(sqr(r2))))
      /(std::pow(-1. + r2,4)*std::pow(1. + r2,2))
      + (0.75*(-1 + 4*sqr(r2) - 3*std::pow(r2,4)
               + 2*std::pow(r2,4)*std::log(sqr(r2))))
      /(std::pow(-1 + r2,3)*std::pow(1 + r2,2))
      + (0.25*std::pow(-1 + r1,2)*(
            1 - 4*r2 + 4*sqr(r2) + 8*std::pow(r2,3) - 5*std::pow(r2,4)
            - 4*std::pow(r2,5) + 6*std::pow(r2,4)*std::log(sqr(r2))))
      /(std::pow(-1 + r2,5)*std::pow(1 + r2,2))
      + (0.1*std::pow(-1. + r1,3)*(
            1.5000000000000002 - 7.*r2 + 9.*sqr(r2) + 15.*std::pow(r2,3)
            - 7.5*std::pow(r2,4) - 9.*std::pow(r2,5)
            - 3.0000000000000004*std::pow(r2,6) + 1.*std::pow(r2,7)
            + 15.*std::pow(r2,4)*std::log(sqr(r2))))
      /(std::pow(-1. + r2,6)*std::pow(1. + r2,2));
}

double f8(double r1, double r2)
{
   if (is_equal(r1, 1., 0.01) && is_equal(r2, 1., 0.01))
      return f8_1_1(r1, r2);

   if (is_equal(r1, 1., 0.01))
      return f8_1_r2(r1, r2);

   if (is_equal(r2, 1., 0.01))
      return f8_1_r2(r2, r1);

   const double r12 = sqr(r1);
   const double r22 = sqr(r2);

   const double result
      = (r1+r2)/((r12-1)*(r22-1))
      + (std::pow(r1,4)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (std::pow(r2,4)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 1.5 * result;
}

} // namespace threshold_loop_functions
} // namespace flexiblesusy
