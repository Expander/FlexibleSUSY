// \file rk.hpp
//    - Project:     SOFTSUSY
//    - Author:      Ben Allanach, Alexander Voigt
//    - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
//    - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
//    - Description: Integration of ODEs by Runge Kutta, minimum finding and
//                 derivative calculation

#ifndef RK_H
#define RK_H

#include <algorithm>
#include <cmath>
#include <functional>

#include <Eigen/Dense>

#include "logger.hpp"
#include "error.hpp"

namespace flexiblesusy {

namespace runge_kutta {

namespace {
/// Returns |a| with sign of b in front
inline double sign(double a, double b) noexcept {
   return b >= 0 ? std::fabs(a) : -std::fabs(a);
}
} // anonymous namespace

/// A single step of Runge Kutta (5th order), input:
/// y and dydx (derivative of y), x is independent variable. yout is value
/// after step. derivs is a user-supplied function
template <typename ArrayType, typename Derivs>
void rungeKuttaStep(const ArrayType& y, const ArrayType& dydx, double x,
		    double h, ArrayType& yout, ArrayType& yerr, Derivs derivs)
{
  const double a2 = 0.2,a3 = 0.3,a4 = 0.6,a5 = 1.0,a6 = 0.875,b21 =
    0.2,b31 = 3.0 / 40.0,b32 = 9.0 / 40.0,b41 = 0.3,b42 = -0.9,b43 = 1.2,
    b51 = -11.0 / 54.0, b52 = 2.5,b53 = -70.0 / 27.0,b54 = 35.0 / 27.0,
    b61 = 1631.0 / 55296.0,b62 = 175.0 / 512.0,b63 = 575.0 / 13824.0,
    b64 = 44275.0 / 110592.0,b65 = 253.0 / 4096.0,c1 = 37.0 / 378.0,
    c3 = 250.0 / 621.0,c4 = 125.0 / 594.0,c6 = 512.0 / 1771.0,
    dc5 = -277.00 / 14336.0;
  const double dc1 = c1-2825.0 / 27648.0,dc3 = c3-18575.0 / 48384.0,
    dc4 = c4-13525.0 / 55296.0,dc6 = c6-0.25;

  ArrayType ytemp = b21 * h * dydx + y;
  const ArrayType ak2 = derivs(x + a2 * h, ytemp);

  // Allowing piece-wise calculating of ytemp for speed reasons
  ytemp = y + h * (b31 * dydx + b32 * ak2);
  const ArrayType ak3 = derivs(x + a3 * h, ytemp);

  ytemp = y + h * (b41 * dydx + b42 * ak2 + b43 * ak3);
  const ArrayType ak4 = derivs(x+a4*h,ytemp);

  ytemp = y + h * (b51 * dydx + b52 * ak2 + b53 * ak3 + b54 * ak4);
  const ArrayType ak5 = derivs(x + a5 * h, ytemp);

  ytemp = y + h * (b61 * dydx + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);
  const ArrayType ak6 = derivs(x + a6 * h, ytemp);

  yout = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
  yerr = h * (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
}

/// organises the variable step-size for Runge-Kutta evolution
template <typename ArrayType, typename Derivs>
void odeStepper(ArrayType& y, const ArrayType& dydx, double *x, double htry,
                double eps, const ArrayType& yscal, double *hdid, double *hnext,
                Derivs derivs, int& max_step_dir)
{
  const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
  const int n = y.size();
  double errmax, h = htry, htemp, xnew;
  ArrayType yerr(n), ytemp(n);

  for (;;) {
    rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, derivs);
    errmax = (yerr / yscal).abs().maxCoeff(&max_step_dir);
    errmax  /= eps;
    if (!std::isfinite(errmax)) {
#ifdef ENABLE_VERBOSE
       ERROR("odeStepper: non-perturbative running at Q = "
             << std::exp(*x) << " GeV of parameter y(" << max_step_dir
             << ") = " << y(max_step_dir) << ", dy(" << max_step_dir
             << ")/dx = " << dydx(max_step_dir));
#endif
       throw NonPerturbativeRunningError(std::exp(*x), max_step_dir, y(max_step_dir));
    }
    if (errmax <= 1.0) break;
    htemp = SAFETY * h * std::pow(errmax, PSHRNK);
    h = (h >= 0.0 ? std::max(htemp ,0.1 * h) : std::min(htemp, 0.1 * h));
    xnew = (*x) + h;
    if (xnew == *x) {
#ifdef ENABLE_VERBOSE
       ERROR("At Q = " << std::exp(*x) << " GeV "
             "stepsize underflow in odeStepper in parameter y("
             << max_step_dir << ") = " << y(max_step_dir) << ", dy("
             << max_step_dir << ")/dx = " << dydx(max_step_dir));
#endif
       throw NonPerturbativeRunningError(std::exp(*x), max_step_dir, y(max_step_dir));
    }
  }
  if (errmax > ERRCON) *hnext = SAFETY * h * std::pow(errmax,PGROW);
  else *hnext = 5.0 * h;
  *x += (*hdid = h);
  y = ytemp;
}

/// Organises integration of 1st order system of ODEs
template <typename ArrayType, typename Derivs, typename Stepper>
void integrateOdes(ArrayType& ystart, double from, double to, double eps,
                   double h1, double hmin, Derivs derivs, Stepper rkqs)
{
  const int nvar = ystart.size();
  const int MAXSTP = 400;
  const double TINY = 1.0e-16;
  double x = from, hnext, hdid, h = sign(h1, to - from);
  ArrayType yscal(nvar), y(ystart), dydx;
  int max_step_dir;

  for (int nstp = 0; nstp < MAXSTP; nstp++) {
    dydx = derivs(x, y);
    yscal = y.abs() + (dydx * h).abs() + TINY;
    if ((x + h - to) * (x + h - from) > 0.0) h = to - x;
    rkqs(y, dydx, &x, h, eps, yscal, &hdid, &hnext, derivs, max_step_dir);

    if ((x - to) * (to - from) >= 0.0) {
      ystart = y;
      return;
    }

    h = hnext;

    if (std::fabs(hnext) <= hmin)
      break;
  }

#ifdef ENABLE_VERBOSE
  ERROR("Bailed out of rk.cpp:too many steps in integrateOdes\n"
        "********** Q = " << std::exp(x) << " *********");
  ERROR("max step in direction of " << max_step_dir);
  for (int i = 0; i < nvar; i++)
    ERROR("y(" << i << ") = " << y(i) << " dydx(" << i <<
          ") = " << dydx(i));
#endif

  throw NonPerturbativeRunningError(std::exp(x), max_step_dir, y(max_step_dir));
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif // RK_H
