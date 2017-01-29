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
#include "pv.hpp"
#include "numerics.h"

#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace threshold_loop_functions {

namespace {
   double sqr(double x) { return x*x; }
   template <typename T> T cube(T x) { return x*x*x; }
   template <typename T> T quad(T x) { return x*x*x*x; }
   template <typename T> T pow5(T x) { return x*x*x*x*x; }
   template <typename T> T pow6(T x) { return x*x*x*x*x*x; }
   template <typename T> T pow7(T x) { return x*x*x*x*x*x*x; }
   template <typename T> T pow8(T x) { return x*x*x*x*x*x*x*x; }
   template <typename T> T pow9(T x) { return x*x*x*x*x*x*x*x*x; }
   template <typename T> T pow10(T x) { return x*x*x*x*x*x*x*x*x*x; }

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

      if (std::fabs(a) < std::numeric_limits<T>::epsilon())
         return is_equal(a, b, prec);

      return std::fabs((a - b)/a) < prec;
   }
}

double F1(double x)
{
   const double x2 = sqr(x);

   if (is_equal(x, 1., 0.01))
      return (16 + 41*x - 44*x2 + 21*cube(x) - 4*quad(x))/30.;

   return x*std::log(x2)/(x2-1);
}

double F2(double x)
{
   const double x2 = sqr(x);

   if (is_equal(x, 1., 0.01))
      return (-5 + 216*x - 226*x2 + 104*cube(x) - 19*quad(x))/70.;

   return 6*x2*(2-2*x2+(1+x2)*std::log(x2))/cube(x2-1);
}

double F3(double x)
{
   const double x2 = sqr(x);

   if (is_equal(x, 1., 0.01))
      return (-27 + 218*x - 142*x2 + 48*cube(x) - 7*quad(x))/90.;

   return 2*x*(5*(1-x2)+(1+4*x2)*std::log(x2))/(3*sqr(x2-1));
}

double F4(double x)
{
   const double x2 = sqr(x);

   if (is_equal(x, 1., 0.01))
      return (31 + 22*x - 42*x2 + 24*cube(x) - 5*quad(x))/30.;

   return 2*x*(x2-1-std::log(x2))/sqr(x2-1);
}

double F5(double x)
{
   const double x2 = sqr(x);
   const double x4 = quad(x);

   if (is_equal(x, 1., 0.01))
      return (13 + 165*x - 174*x2 + 81*cube(x) - 15*x4)/70.;

   return 3*x*(1-x4+2*x2*std::log(x2))/cube(1-x2);
}

double F6(double x)
{
   const double x2 = sqr(x);

   if (is_equal(x, 1., 0.01))
      return (-103 + 128*x - 26*x2 + quad(x))/120.;

   return (x2-3)/(4*(1-x2)) + x2*(x2-2)/(2*sqr(1.-x2))*std::log(x2);
}

double F7(double x)
{
   const double x2 = sqr(x);
   const double x4 = quad(x);

   if (is_equal(x, 1., 0.01))
      return -1.8642857142857139 + 2.057142857142856*x
         + 1.7142857142857153*x2 - 1.1428571428571432*cube(x)
         + 0.2357142857142858*x4;

   return (-3*(x4-6*x2+1.))/(2*sqr(x2-1))
      + (3*x4*(x2-3.))/(cube(x2-1.))*std::log(x2);
}

/// F8(x1,x2) in the limit x1 -> 1 and x2 -> 1
static double F8_1_1(double x1, double x2)
{
   return (-2.7015873015872867 + 13.815873015872924*x2
           - 28.54761904761879*sqr(x2) + 29.898412698412283*cube(x2)
           - 15.82222222222179*quad(x2)
           + 2.9904761904758974*pow5(x2)
           + 0.7031746031747295*pow6(x2)
           - 0.39365079365082534*pow7(x2)
           + 0.05714285714286069*pow8(x2)
           + sqr(x1)*(-0.2999999999999769 + 5.352380952380799*x2
                      - 23.09047619047574*sqr(x2)
                      + 47.08571428571352*cube(x2)
                      - 54.404761904761074*quad(x2)
                      + 38.03809523809467*pow5(x2)
                      - 16.18571428571403*pow6(x2)
                      + 3.923809523809459*pow7(x2)
                      - 0.41904761904761173*pow8(x2))
           + quad(x1)*(0.057142857142858015 + 0.050793650793644535*x2
                             - 1.193650793650774*sqr(x2)
                             + 3.352380952380916*cube(x2)
                             - 4.534920634920595*quad(x2)
                             + 3.511111111111082*pow5(x2)
                             - 1.6095238095237963*pow6(x2)
                             + 0.41269841269840923*pow7(x2)
                             - 0.04603174603174563*pow8(x2))
           + cube(x1)*(-0.16507936507937204 - 1.0158730158729674*x2
                             + 7.961904761904614*sqr(x2)
                             - 19.555555555555294*cube(x2)
                             + 24.926984126983836*quad(x2)
                             - 18.590476190475986*pow5(x2)
                             + 8.292063492063399*pow6(x2)
                             - 2.0825396825396587*pow7(x2)
                             + 0.22857142857142587*pow8(x2))
           + x1*(3.009523809523779 - 17.269841269841073*x2
                 + 43.13650793650738*sqr(x2)
                 - 61.71428571428479*cube(x2)
                 + 55.83492063491965*quad(x2)
                 - 33.01587301587234*pow5(x2)
                 + 12.533333333333038*pow6(x2)
                 - 2.793650793650719*pow7(x2)
                 + 0.279365079365071*pow8(x2)))/quad(-1. + x2);
}

/// F8(x1,x2) in the limit x1 -> 1
static double F8_1_x2(double x1, double x2)
{
   return -2. + (4.*(-1. + x1)*(-0.5 + 2.*sqr(x2) - 1.5*quad(x2)
                                + 1.*quad(x2)*std::log(sqr(x2))))
      /cube(-1. + sqr(x2))
      + (2. - 2.*sqr(x2) + 2.*quad(x2)*std::log(sqr(x2)))
      /sqr(-1. + sqr(x2))
      - (1.3333333333333333*cube(-1. + x1)*(
            -1.*sqr(x2) - 9.*quad(x2) + 9.*pow6(x2)
            + 1.*pow8(x2)
            + (-6.*quad(x2) - 6.*pow6(x2))*std::log(sqr(x2))))
      /pow5(-1. + sqr(x2))
      + (2.*sqr(-1. + x1)*(
            -0.16666666666666666 + 1.5*sqr(x2) + 1.5*quad(x2)
            - 2.8333333333333335*pow6(x2)
            + (3.*quad(x2) + 1.*pow6(x2))*std::log(sqr(x2))))
      /quad(-1. + sqr(x2))
      + (0.26666666666666666*quad(-1. + x1)*(
            0.25 + 2.5*sqr(x2) + 80.*quad(x2) - 47.5*pow6(x2)
            - 36.25*pow8(x2) + 1.*pow10(x2)
            + (37.5*quad(x2) + 75.*pow6(x2)
               + 7.5*pow8(x2))*std::log(sqr(x2))))
      /pow6(-1. + sqr(x2));
}

// F8(x1,x2) in the limit x1 -> x2
static double F8_x1_x2(double x1, double x2)
{
   return -2. + (2.*(x1 - x2)*(3.*x2 - 4.*cube(x2) + 1.*pow5(x2) +
        2.*x2*std::log(sqr(x2))))/cube(-1. + sqr(x2)) +
   (2.*sqr(x2)*(-1. + sqr(x2) - 2.*std::log(sqr(x2)) +
        sqr(x2)*std::log(sqr(x2))))/sqr(-1. + sqr(x2)) -
   (0.33333333333333326*sqr(x1 - x2)*
      (17.000000000000007 - 9.000000000000007*sqr(x2) -
        9.*quad(x2) + 1.*pow6(x2) +
        6.000000000000002*std::log(sqr(x2)) +
        18.000000000000004*sqr(x2)*std::log(sqr(x2))))/
    quad(-1. + sqr(x2)) -
   (1.333333333333333*cube(x1 - x2)*
      (-1. - 9.000000000000002*sqr(x2) +
        9.000000000000002*quad(x2) + 1.*pow6(x2) -
        6.*sqr(x2)*std::log(sqr(x2)) -
        6.000000000000003*quad(x2)*std::log(sqr(x2))))/
    (x2*pow5(-1. + sqr(x2)));
}

double F8(double x1, double x2)
{
   if (is_equal(x1, 1., 0.01) && is_equal(x2, 1., 0.01))
      return F8_1_1(x1, x2);

   if (is_equal(x1, 1., 0.01))
      return F8_1_x2(x1, x2);

   if (is_equal(x2, 1., 0.01))
      return F8_1_x2(x2, x1);

   if (is_equal(x1, x2, 0.00001))
      return F8_x1_x2(x1, x2);

   const double x12 = sqr(x1);
   const double x22 = sqr(x2);

   return -2. + 2./(x12-x22)
      *(quad(x1)/(x12-1.)*std::log(x12)
        -quad(x2)/(x22-1.)*std::log(x22));
}

/// F9(x1,x2) in the limit x1 -> 1 and x2 -> 1
static double F9_1_1(double x1, double x2)
{
   return 8.223809523809523 - 12.863492063492064*x2
      + 10.580952380952382*sqr(x2) - 4.609523809523809*cube(x2)
      + 0.8349206349206348*quad(x2)
      + x1*(-12.863492063492064 + 26.260317460317456*x2
            - 24.609523809523807*sqr(x2) + 11.530158730158728*cube(x2)
            - 2.184126984126984*quad(x2))
      + cube(x1)*(-4.60952380952381 + 11.53015873015873*x2
                        - 11.961904761904762*sqr(x2)
                        + 5.942857142857143*cube(x2)
                        - 1.1682539682539683*quad(x2))
      + quad(x1)*(0.8349206349206351 - 2.1841269841269844*x2
                        + 2.3190476190476197*sqr(x2)
                        - 1.1682539682539685*cube(x2)
                        + 0.23174603174603178*quad(x2))
      + sqr(x1)*(10.580952380952379 - 24.609523809523804*x2
                 + 24.6047619047619*sqr(x2)
                 - 11.96190476190476*cube(x2)
                 + 2.319047619047619*quad(x2));
}

/// F9(x1,x2) in the limit x1 -> 1
static double F9_1_x2(double x1, double x2)
{
   return (-2.*(-1. + x1)*cube(-1. + sqr(x2))*(
              -1. + 1.*quad(x2) - 2.*sqr(x2)*std::log(sqr(x2)))
           + 2.*quad(-1. + sqr(x2))*(
              1. - 1.*sqr(x2) + 1.*sqr(x2)*std::log(sqr(x2)))
           - 1.3333333333333333*cube(-1. + x1)*(
              -1. + sqr(x2))*(
                 -1. - 9.*sqr(x2) + 9.*quad(x2) + 1.*pow6(x2)
                 + (-6.*sqr(x2) - 6.*quad(x2))*std::log(sqr(x2)))
           + 0.3333333333333333*sqr(-1. + x1)*sqr(-1. + sqr(x2))
           *(5. + 15.*sqr(x2) - 21.*quad(x2) + 1.*pow6(x2)
             + (18.*sqr(x2) + 6.*quad(x2))*std::log(sqr(x2)))
           - 0.06666666666666667*quad(-1. + x1)*(
              -16. - 305.*sqr(x2) + 170.*quad(x2) + 160.*pow6(x2)
              - 10.*pow8(x2) + 1.*pow10(x2)
              + (-150.*sqr(x2) - 300.*quad(x2)
                 - 30.*pow6(x2))*std::log(sqr(x2))))
      /pow6(-1. + sqr(x2));
}

/// F9(x1,x2) in the limit x1 -> x2
static double F9_x1_x2(double x1, double x2)
{
   return (2.*(-1. + 1.*sqr(x2) - 1.*std::log(sqr(x2))))/
      sqr(-1. + sqr(x2)) -
      (2.*(x1 - x2)*(-1. + 1.*quad(x2) -
                     2.*sqr(x2)*std::log(sqr(x2))))/
      (x2*cube(-1. + sqr(x2))) -
      (1.3333333333333333*cube(x1 - x2)*
       (-1.0000000000000002 - 9.*sqr(x2) + 9.*quad(x2) +
        1.*pow6(x2) - 6.*sqr(x2)*std::log(sqr(x2)) -
        6.*quad(x2)*std::log(sqr(x2))))/
      (x2*pow5(-1. + sqr(x2))) +
      (1.6666666666666665*sqr(x1 - x2)*
       (0.2000000000000001 - 4.200000000000001*sqr(x2) +
        3.000000000000001*quad(x2) + 1.*pow6(x2) -
        1.2000000000000002*sqr(x2)*std::log(sqr(x2)) -
        3.6000000000000005*quad(x2)*std::log(sqr(x2))))/
      (sqr(x2)*quad(-1. + sqr(x2)));
}

double F9(double x1, double x2)
{
   if (is_equal(x1, 1., 0.01) && is_equal(x2, 1., 0.01))
      return F9_1_1(x1, x2);

   if (is_equal(x1, 1., 0.01))
      return F9_1_x2(x1, x2);

   if (is_equal(x2, 1., 0.01))
      return F9_1_x2(x2, x1);

   if (is_equal(x1, x2, 0.00001))
      return F9_x1_x2(x1, x2);

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
   if (is_equal(r, 1., 0.01))
      return (-81 + 464*r + 270*sqr(r)
              - 208*cube(r) + 45*quad(r))/490.;

   const double r2 = sqr(r);

   return (6*(r2+3)*r2)/(7*sqr(r2-1))
      + (6*(r2-5)*quad(r)*std::log(r2))/(7*cube(r2-1));
}

double f2(double r)
{
   if (is_equal(r, 1., 0.01))
      return (-285 + 1616*r + 1230*sqr(r)
              - 848*cube(r) + 177*quad(r))/1890.;

   const double r2 = sqr(r);

   return (2*(r2+11)*r2)/(9*sqr(r2-1))
      + (2*(5*r2-17)*quad(r)*std::log(r2))/(9*cube(r2-1));
}

double f3(double r)
{
   if (is_equal(r, 1., 0.01))
      return (849 - 1184*r + 1566*sqr(r)
              - 736*cube(r) + 135*quad(r))/630.;

   const double r2 = sqr(r);
   const double r4 = quad(r);

   return (2*(r4+9*r2+2))/(3*sqr(r2-1))
      + (2*(r4-7*r2-6)*r2*std::log(r2))/(3*cube(r2-1));
}

double f4(double r)
{
   const double r2 = sqr(r);
   const double r4 = quad(r);

   if (is_equal(r, 1., 0.01))
      return (2589 - 3776*r + 4278*r2 - 1984*cube(r) + 363*r4)/1470.;

   return (2*(5*r4+25*r2+6))/(7*sqr(r2-1))
      + (2*(r4-19*r2-18)*r2*std::log(r2))/(7*cube(r2-1));
}

/// f5(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f5_1_1(double r1, double r2)
{
   return 0.772943722943723
      - 0.5524891774891774*r2
      + 0.7870670995670994*sqr(r2)
      - 0.3316558441558441*cube(r2)
      + 0.056277056277056266*quad(r2)
      + r1*(-0.5524891774891774
            + 1.0700757575757573*r2
            - 0.6625541125541123*sqr(r2)
            + 0.22483766233766228*cube(r2)
            - 0.03344155844155843*quad(r2))
      + cube(r1)*(-0.33165584415584404
                        + 0.22483766233766223*r2
                        - 0.08755411255411245*sqr(r2)
                        + 0.01650432900432896*cube(r2)
                        - 0.0007034632034631958*quad(r2))
      + quad(r1)*(0.05627705627705626
                        - 0.03344155844155841*r2
                        + 0.010281385281385256*sqr(r2)
                        - 0.0007034632034631921*cube(r2)
                        - 0.0002705627705627725*quad(r2))
      + sqr(r1)*(0.7870670995670994 - 0.6625541125541123*r2
                 + 0.32061688311688297*sqr(r2)
                 - 0.08755411255411248*cube(r2)
                 + 0.01028138528138527*quad(r2));
}

/// f5(r1,r2) in the limit r1 -> 1
static double f5_1_r2(double r1, double r2)
{
   return (-0.025*cube(-1. + r1)*(
              4. - 17.*r2 + 4.*sqr(r2)
              - 25.*cube(r2)
              - 20.*quad(r2)
              + 41.*pow5(r2)
              + 12.*pow6(r2)
              + pow7(r2)
              + (-30.*cube(r2) - 30.*pow5(r2))
              *std::log(sqr(r2))))/(pow6(-1. + r2)*sqr(1. + r2))
      - (0.125*sqr(-1. + r1)*(
            1. - 4.*r2 + sqr(r2)
            - 4.*cube(r2)
            - 5.*quad(r2)
            + 8.*pow5(r2)
            + 3.*pow6(r2)
            + (-6.*cube(r2) - 6.*pow5(r2))
            *std::log(sqr(r2))))/(pow5(-1. + r2)*sqr(1. + r2))
      + (0.75*(-1 + r2 + 2*sqr(r2)
               - quad(r2)
               - pow5(r2)
               + (cube(r2) + pow5(r2))
               *std::log(sqr(r2))))/(cube(-1 + r2)*sqr(1 + r2))
      + (0.25*(-1. + r1)*(
            1. - 1.*r2 - 2.*sqr(r2) + 8.*cube(r2)
            + quad(r2) - 7.*pow5(r2)
            + (3.*cube(r2) + 3.*pow5(r2))
            *std::log(sqr(r2))))/(quad(-1. + r2)*sqr(1. + r2))
      + (0.05*quad(-1. + r1)*(
            -1. + 4.5*r2 + 2.*sqr(r2)
            + 16.5*cube(r2) - 16.5*pow5(r2) - 2.*pow6(r2)
            - 4.5*pow7(r2) + 1.*pow8(r2)
            + (15.*cube(r2) + 15.*pow5(r2))
            *std::log(sqr(r2))))/(pow7(-1. + r2)*sqr(1. + r2));
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
      + (cube(r1)*(r12+1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (cube(r2)*(sqr(r2)+1)*std::log(sqr(r2)))/((r1-r2)*sqr(sqr(r2)-1));

   return 0.75 * result;
}

/// f6(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f6_1_1(double r1, double r2)
{
   return 1 + (4*(-1 + r1))/7. - (2*sqr(-1 + r1))/35.
      - cube(-1 + r1)/70. + (9*quad(-1 + r1))/490.
      + (0.5714285714285714 - (2*(-1 + r1))/35. - sqr(-1 + r1)/70.
         + (9*cube(-1 + r1))/490.
         - (3*quad(-1 + r1))/245.)*(-1 + r2)
      + (-0.05714285714285714 + (1 - r1)/70. + (9*sqr(-1 + r1))/490.
         - (3*cube(-1 + r1))/245.
         + quad(-1 + r1)/147.)*sqr(-1 + r2)
      + (-0.014285714285714285 + (9*(-1 + r1))/490.
         - (3*sqr(-1 + r1))/245. + cube(-1 + r1)/147.
         - quad(-1 + r1)/294.)*cube(-1 + r2)
      + (0.018367346938775512 - (3*(-1 + r1))/245. + sqr(-1 + r1)/147.
         - cube(-1 + r1)/294.
         + (5*quad(-1 + r1))/3234.)*quad(-1 + r2);
}

/// f6(r1,r2) in the limit r1 -> 1
static double f6_1_r2(double r1, double r2)
{
   return ((-1 + sqr(r2))*(
              -3 + 16*r2 + 33*sqr(r2) - 332*cube(r2)
              + 573*quad(r2) - 584*pow5(r2) + 297*pow6(r2)
              - 60*pow7(r2) - 2*cube(r1)*(
                 8 - 49*r2 + 121*sqr(r2) - 129*cube(r2)
                 - 99*quad(r2) + 28*pow5(r2))
              + quad(r1)*(
                 3 - 18*r2 + 43*sqr(r2) - 42*cube(r2)
                 - 47*quad(r2) + pow6(r2))
              - 2*r1*(-8 + 10*r2 + 106*sqr(r2) - 359*cube(r2)
                      + 211*quad(r2) - 101*pow5(r2)
                      + 21*pow6(r2))
              - 2*sqr(r1)*(-15 + 98*r2 - 264*sqr(r2)
                                  + 331*cube(r2) + 106*quad(r2)
                                  - 99*pow5(r2) + 23*pow6(r2)))
           + 60*pow5(r2)*(5 + quad(r1) + cube(r1)*(-5 + r2)
                                - 10*r2 + 10*sqr(r2) - 5*cube(r2)
                                + quad(r2) + sqr(r1)*(
                                   10 - 5*r2 + sqr(r2))
                                + r1*(-10 + 10*r2 - 5*sqr(r2)
                                      + cube(r2)))*std::log(sqr(r2)))
      /(70.*pow7(-1 + r2)*sqr(1 + r2));
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
      + (pow5(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (pow5(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6./7. * result;
}

/// f7(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f7_1_1(double r1, double r2)
{
   return (15700 - 14411*r2 + 7850*sqr(r2) - 2498*cube(r2)
           + 355*quad(r2)
           + sqr(r1)*(7850 - 2558*r2 - 750*sqr(r2) + 940*cube(r2)
                      - 235*quad(r2))
           + quad(r1)*(355 + 65*r2 - 235*sqr(r2) + 142*cube(r2)
                             - 30*quad(r2))
           + r1*(-14411 + 8375*r2 - 2558*sqr(r2) + 180*cube(r2)
                 + 65*quad(r2))
           + cube(r1)*(-2498 + 180*r2 + 940*sqr(r2)
                             - 645*cube(r2) + 142*quad(r2)))
      /2310.;
}

/// f7(r1,r2) in the limit r1 -> 1
static double f7_1_r2(double r1, double r2)
{
   return (-10*(-1 + r1)*cube(-1 + r2)*(
              2 - 5*r2 - 4*sqr(r2) + 4*cube(r2) + 2*quad(r2)
              + pow5(r2) - 6*cube(r2)*std::log(sqr(r2)))
           - 30*quad(-1 + r2)*(
              1 - 2*r2 - 2*sqr(r2) + 2*cube(r2) + quad(r2)
              - 2*cube(r2)*std::log(sqr(r2)))
           + 10*sqr(-1 + r1)*sqr(-1 + r2)*(
              -1 + 3*r2 + 3*sqr(r2) - 3*quad(r2) - 3*pow5(r2)
              + pow6(r2) + 6*cube(r2)*std::log(sqr(r2)))
           + 2*cube(-1 + r1)*(-1 + r2)*(
              -2 + 6*r2 + 18*sqr(r2) + 15*cube(r2) - 30*quad(r2)
              - 18*pow5(r2) + 14*pow6(r2) - 3*pow7(r2)
              + 30*cube(r2)*std::log(sqr(r2)))
           + quad(-1 + r1)*(
              -1 + 48*sqr(r2) + 42*cube(r2) - 90*quad(r2)
              - 24*pow5(r2) + 40*pow6(r2) - 18*pow7(r2)
              + 3*pow8(r2) + 60*cube(r2)*std::log(sqr(r2))))
      /(10.*pow7(-1 + r2)*sqr(1 + r2));
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
      + (cube(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (cube(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 6. * result;
}

/// f8(r1,r2) in the limit r1 -> 1 and r2 -> 1
static double f8_1_1(double r1, double r2)
{
   return 1 - sqr(-1 + r1)/10. + (3*cube(-1 + r1))/40.
      - (3*quad(-1 + r1))/70.
      + ((1 - r1)/10. + (3*sqr(-1 + r1))/40.
         - (3*cube(-1 + r1))/70.
         + (3*quad(-1 + r1))/140.)*(-1 + r2)
      + (-0.1 + (3*(-1 + r1))/40. - (3*sqr(-1 + r1))/70.
         + (3*cube(-1 + r1))/140.
         - quad(-1 + r1)/105.)*sqr(-1 + r2)
      + (0.075 - (3*(-1 + r1))/70. + (3*sqr(-1 + r1))/140.
         - cube(-1 + r1)/105.
         + quad(-1 + r1)/280.)*cube(-1 + r2)
      + (-0.04285714285714286 + (3*(-1 + r1))/140.
         - sqr(-1 + r1)/105. + cube(-1 + r1)/280.
         - quad(-1 + r1)/1155.)*quad(-1 + r2);
}

/// f8(r1,r2) in the limit r1 -> 1
static double f8_1_r2(double r1, double r2)
{
   return (30*quad(-1 + r2)*(
              -1 + 4*sqr(r2) - 3*quad(r2)
              + 2*quad(r2)*std::log(sqr(r2)))
           + 10*sqr(-1 + r1)*sqr(-1 + r2)*(
              1 - 4*r2 + 4*sqr(r2) + 8*cube(r2)
              - 5*quad(r2) - 4*pow5(r2)
              + 6*quad(r2)*std::log(sqr(r2)))
           + 10*(-1 + r1)*cube(-1 + r2)*(
              1 - 4*r2 + 4*sqr(r2) + 8*cube(r2)
              - 5*quad(r2) - 4*pow5(r2)
              + 6*quad(r2)*std::log(sqr(r2)))
           + 2*cube(-1 + r1)*(-1 + r2)*(
              3 - 14*r2 + 18*sqr(r2) + 30*cube(r2)
              - 15*quad(r2) - 18*pow5(r2) - 6*pow6(r2)
              + 2*pow7(r2) + 30*quad(r2)*std::log(sqr(r2)))
           + quad(-1 + r1)*(
              3 - 16*r2 + 24*sqr(r2) + 48*cube(r2)
              - 48*pow5(r2) - 24*pow6(r2) + 16*pow7(r2)
              - 3*pow8(r2) + 60*quad(r2)*std::log(sqr(r2))))
      /(40.*pow7(-1 + r2)*sqr(1 + r2));
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
      + (quad(r1)*std::log(r12))/(sqr(r12-1)*(r1-r2))
      - (quad(r2)*std::log(r22))/((r1-r2)*sqr(r22-1));

   return 1.5 * result;
}

/// Iabc(a,a,a)
static double Iaaa(double a, double b, double c)
{
   return (151.*quad(a) + 13.*sqr(b)*sqr(c) - 128.*cube(a)*(b + c) - 40.*a*b*c*(b + c)
           + sqr(a)*(37.*sqr(b) + 128.*b*c + 37.*sqr(c))) / (60.*pow6(a));
}

/// Iabc(a,a,c)
static double Iaac(double a, double b, double c)
{
   return ((sqr(a) - sqr(c))
           * (17.*pow6(a) - 16.*pow5(a)*b - 40.*cube(a)*b*sqr(c)
              + 8.*a*b*quad(c) - sqr(b)*quad(c) + quad(a)*(5.*sqr(b) + 8.*sqr(c))
              + sqr(a)*(20.*sqr(b)*sqr(c) - quad(c)))
           - 6.*sqr(a)*sqr(c) * log(sqr(a)/sqr(c))
           * (6.*quad(a) - 8.*cube(a)*b + 3.*sqr(a)*(sqr(b) - sqr(c)) + sqr(c)*(sqr(b) + sqr(c))))
      / (6.*sqr(a)*quad(sqr(a) - sqr(c)));
}

/// Iabc(a,a,0)
double Iaa0(double a, double b)
{
   return (17.*sqr(a) - 16.*a*b + 5.*sqr(b)) / (6.*quad(a));
}

/// Iabc(0,b,c)
double I0bc(double b, double c)
{
   return log(sqr(b/c))/(sqr(b) - sqr(c));
}

double Iabc(double a, double b, double c) {
   if ((is_zero(a) && is_zero(b) && is_zero(c)) ||
       (is_zero(a) && is_zero(b)) ||
       (is_zero(a) && is_zero(c)) ||
       (is_zero(b) && is_zero(c)))
      return 0.;

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01) && is_equal_rel(std::abs(a), std::abs(c), 0.01))
      return Iaaa(std::abs(a),std::abs(b),std::abs(c));

   if (is_equal_rel(std::abs(a), std::abs(b), 0.01)) {
      if (is_zero(c))
         return Iaa0(std::abs(a),std::abs(b));
      return Iaac(std::abs(a),std::abs(b),c);
   }

   if (is_equal_rel(std::abs(b), std::abs(c), 0.01)) {
      if (is_zero(a))
         return Iaa0(std::abs(b),std::abs(c));
      return Iaac(std::abs(b),std::abs(c),a);
   }

   if (is_equal_rel(std::abs(a), std::abs(c), 0.01)) {
      if (is_zero(b))
         return Iaa0(std::abs(a),std::abs(c));
      return Iaac(std::abs(a),std::abs(c),b);
   }

   if (is_zero(a))
      return I0bc(b,c);

   if (is_zero(b))
      return I0bc(c,a);

   if (is_zero(c))
      return I0bc(a,b);

   return ( (sqr(a * b) * log(sqr(a / b))
           + sqr(b * c) * log(sqr(b / c))
           + sqr(c * a) * log(sqr(c / a)))
           / ((sqr(a) - sqr(b)) * (sqr(b) - sqr(c)) * (sqr(a) - sqr(c))) );
}

/**
 * B0 function for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 * @param scale \f$Q\f$ (not squared)
 *
 * @return \f$B_0(p=0,m_1,m_2,Q)\f$
 */
double B0(double m1, double m2, double scale)
{
   return passarino_veltman::ReB0(0, m1*m1, m2*m2, scale*scale);
}

/**
 * B0' function for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 *
 * @return \f$B_0'(p=0,m_1,m_2)\f$
 */
double DB0(double m1, double m2)
{
   const double m12 = sqr(m1);
   const double m14 = sqr(m12);
   const double m22 = sqr(m2);
   const double m24 = sqr(m22);

   if (is_zero(m12) || is_zero(m22))
      return 0.;

   if (is_equal_rel(m12, m22, 1e-3))
      return 1./(6. * m22);

   return (m14 - m24 + 2*m12*m22*std::log(m22/m12))/
      (2*cube(m12 - m22));
}

/**
 * C0 function for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * \f$C_0(0,m_1,m_2,m_3) = -I_{abc}(m_1,m_2,m_3)\f$
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 * @param m3 \f$m_3\f$ (not squared)
 *
 * @return \f$C_0(p=0,m_1,m_2,m_3)\f$
 */
double C0(double m1, double m2, double m3)
{
   return c0(m1, m2, m3);
}

double D0(double m1, double m2, double m3, double m4)
{
   return d0(m1,m2,m3,m4);
}

/**
 * \f$\tilde{D}_2\f$ for zero momentum, arxiv:0901.2065 Eq (131)
 *
 * @param m1  \f$m_1\f$ (not squared)
 * @param m2  \f$m_2\f$ (not squared)
 * @param m3  \f$m_3\f$ (not squared)
 * @param m4  \f$m_4\f$ (not squared)
 *
 * @return \f$\tilde{D}_2(m_1,m_2,m_3,m_4)\f$
 */
double D2t(double m1, double m2, double m3, double m4)
{
   return C0(m2, m3, m4) + m1*m1 * D0(m1, m2, m3, m4);
}

/**
 * \f$\tilde{D}_4\f$ for zero momentum, arxiv:0901.2065 Eq (131)
 *
 * @param m1  \f$m_1\f$ (not squared)
 * @param m2  \f$m_2\f$ (not squared)
 * @param m3  \f$m_3\f$ (not squared)
 * @param m4  \f$m_4\f$ (not squared)
 * @param scale \f$Q\f$ (not squared)
 *
 * @return \f$\tilde{D}_4(m_1,m_2,m_3,m_4,Q)\f$
 */
double D4t(double m1, double m2, double m3, double m4, double scale)
{
   return B0(m3, m4, scale) + (m1*m1 + m2*m2) * C0(m2, m3, m4)
      + quad(m1) * D0(m1, m2, m3, m4);
}

/**
 * \f$W\f$ for zero momentum, arxiv:0901.2065 Eq (130)
 *
 * @param m1 \f$m_1\f$ (not squared)
 * @param m2 \f$m_2\f$ (not squared)
 * @param scale \f$Q\f$ (not squared)
 *
 * @return \f$W(m_1,m_2,Q)\f$
 */
double W(double m1, double m2, double scale)
{
   const double m12 = sqr(m1);
   const double m14 = sqr(m12);
   const double m22 = sqr(m2);
   const double m24 = sqr(m22);
   const double m26 = m24 * m22;
   const double Q2  = sqr(scale);

   if (is_zero(m12) || is_zero(m22) || is_zero(Q2))
      return 0.;

   if (is_equal_rel(m12,m22,1e-3))
      return 2./3. - 2. * std::log(Q2/m22);

   return (- 2*std::log(Q2/m12)
           - std::log(m22/m12)*(2*m26 - 6*m12*m24)/cube(m12 - m22)
           - (m14 - 6*m22*m12 + m24)/sqr(m12 - m22));
}

} // namespace threshold_loop_functions
} // namespace flexiblesusy
