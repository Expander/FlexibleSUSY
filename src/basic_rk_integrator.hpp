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

/**
 * @file basic_rk_integrator.hpp
 * @brief Integration of ODEs by Runge-Kutta
 * @author Ben Allanach, Alexander Voigt
 *
 * The implementation of the Runge-Kutta routines have been derived
 * from SOFTSUSY [hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305].
 */

#ifndef BASIC_RK_INTEGRATOR_H
#define BASIC_RK_INTEGRATOR_H

#include <algorithm>
#include <cmath>
#include <functional>

#include "error.hpp"
#include "logger.hpp"

namespace flexiblesusy {

namespace runge_kutta {

namespace {
/// Returns |a| with sign of b in front
inline double sign(double a, double b) noexcept {
   return b >= 0 ? std::fabs(a) : -std::fabs(a);
}
} // anonymous namespace

template <typename StateType, typename Derivs>
class Basic_rk_stepper {
public:
   double step(StateType&, const StateType&, double*, double,
               double, const StateType&, Derivs, int&) const;
private:
   void runge_kutta_step(const StateType&, const StateType&, double,
                         double, StateType&, StateType&, Derivs) const;
};

template <typename StateType, typename Derivs>
void Basic_rk_stepper<StateType, Derivs>::runge_kutta_step(
   const StateType& y, const StateType& dydx, double x,
   double h, StateType& yout, StateType& yerr, Derivs derivs) const
{
   const double a2 = 0.2;
   const double a3 = 0.3;
   const double a4 = 0.6;
   const double a5 = 1.0;
   const double a6 = 0.875;
   const double b21 = 0.2;
   const double b31 = 3.0 / 40.0;
   const double b32 = 9.0 / 40.0;
   const double b41 = 0.3;
   const double b42 = -0.9;
   const double b43 = 1.2;
   const double b51 = -11.0 / 54.0;
   const double b52 = 2.5;
   const double b53 = -70.0 / 27.0;
   const double b54 = 35.0 / 27.0;
   const double b61 = 1631.0 / 55296.0;
   const double b62 = 175.0 / 512.0;
   const double b63 = 575.0 / 13824.0;
   const double b64 = 44275.0 / 110592.0;
   const double b65 = 253.0 / 4096.0;
   const double c1 = 37.0 / 378.0;
   const double c3 = 250.0 / 621.0;
   const double c4 = 125.0 / 594.0;
   const double c6 = 512.0 / 1771.0;
   const double dc5 = -277.00 / 14336.0;
   const double dc1 = c1 - 2825.0 / 27648.0;
   const double dc3 = c3 - 18575.0 / 48384.0;
   const double dc4 = c4 - 13525.0 / 55296.0;
   const double dc6 = c6 - 0.25;

   StateType ytemp = b21 * h * dydx + y;
   const StateType ak2 = derivs(x + a2 * h, ytemp);

   // Allowing piece-wise calculating of ytemp for speed reasons
   ytemp = y + h * (b31 * dydx + b32 * ak2);
   const StateType ak3 = derivs(x + a3 * h, ytemp);

   ytemp = y + h * (b41 * dydx + b42 * ak2 + b43 * ak3);
   const StateType ak4 = derivs(x+a4*h,ytemp);

   ytemp = y + h * (b51 * dydx + b52 * ak2 + b53 * ak3 + b54 * ak4);
   const StateType ak5 = derivs(x + a5 * h, ytemp);

   ytemp = y + h * (b61 * dydx + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);
   const StateType ak6 = derivs(x + a6 * h, ytemp);

   yout = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
   yerr = h * (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
}

template <typename StateType, typename Derivs>
double Basic_rk_stepper<StateType,Derivs>::step(
   StateType& y, const StateType& dydx, double *x, double htry,
   double eps, const StateType& yscal, Derivs derivs,
   int& max_step_dir) const
{
   const double SAFETY = 0.9;
   const double PGROW = -0.2;
   const double PSHRNK = -0.25;
   const double ERRCON = 1.89e-4;
   const int n = y.size();

   double errmax;
   double h = htry;
   StateType yerr(n);
   StateType ytemp(n);

   for (;;) {
      runge_kutta_step(y, dydx, *x, h, ytemp, yerr, derivs);
      errmax = (yerr / yscal).abs().maxCoeff(&max_step_dir);
      errmax /= eps;
      if (!std::isfinite(errmax)) {
#ifdef ENABLE_VERBOSE
         ERROR("Basic_rk_stepper::step: non-perturbative running at Q = "
               << std::exp(*x) << " GeV of parameter y(" << max_step_dir
               << ") = " << y(max_step_dir) << ", dy(" << max_step_dir
               << ")/dx = " << dydx(max_step_dir));
#endif
         throw NonPerturbativeRunningError(std::exp(*x), max_step_dir,
                                           y(max_step_dir));
      }
      if (errmax <= 1.0) break;
      const double htemp = SAFETY * h * std::pow(errmax, PSHRNK);
      h = (h >= 0.0 ? std::max(htemp, 0.1 * h) : std::min(htemp, 0.1 * h));
      if (*x + h == *x) {
#ifdef ENABLE_VERBOSE
         ERROR("At Q = " << std::exp(*x) << " GeV "
               "stepsize underflow in Basic_rk_stepper::step in parameter y("
               << max_step_dir << ") = " << y(max_step_dir) << ", dy("
               << max_step_dir << ")/dx = " << dydx(max_step_dir));
#endif
         throw NonPerturbativeRunningError(std::exp(*x), max_step_dir,
                                           y(max_step_dir));
      }
   }
   *x += h;
   y = ytemp;

   return errmax > ERRCON ? SAFETY * h * std::pow(errmax,PGROW) : 5.0 * h;
}

template <typename StateType,
          typename Derivs
          = std::function<StateType(double, const StateType&)>,
          typename Stepper = Basic_rk_stepper<StateType,Derivs> >
class Basic_rk_integrator {
public:
   void operator()(double, double, StateType&, Derivs, double) const;
private:
   Stepper stepper{};
};

template <typename StateType, typename Derivs, typename Stepper>
void Basic_rk_integrator<StateType, Derivs, Stepper>::operator()(
   double start, double end, StateType& ystart, Derivs derivs,
   double tolerance) const
{
   const double guess = (start - end) * 0.1; // first step size
   const double hmin = (start - end) * tolerance * 1.0e-5;
   const int nvar = ystart.size();
   const int MAXSTP = 400;
   const double TINY = 1.0e-16;

   double x = start;
   double h = sign(guess, end - start);

   StateType yscal(nvar);
   StateType y(ystart);
   StateType dydx;
   int max_step_dir;

   for (int nstp = 0; nstp < MAXSTP; ++nstp) {
      dydx = derivs(x, y);
      yscal = y.abs() + (dydx * h).abs() + TINY;

      if ((x + h - end) * (x + h - start) > 0.0) {
         h = end - x;
      }

      const double hnext = stepper.step(y, dydx, &x, h, tolerance,
                                        yscal, derivs, max_step_dir);

      if ((x - end) * (end - start) >= 0.0) {
         ystart = y;
         return;
      }

      h = hnext;

      if (std::fabs(hnext) <= hmin)
         break;
   }

#ifdef ENABLE_VERBOSE
   ERROR("Basic_rkintegrator: too many steps\n"
         "********** Q = " << std::exp(x) << " *********");
   ERROR("max step in direction of " << max_step_dir);
   for (int i = 0; i < nvar; i++)
      ERROR("y(" << i << ") = " << y(i) << " dydx(" << i <<
            ") = " << dydx(i));
#endif

   throw NonPerturbativeRunningError(std::exp(x), max_step_dir,
                                     y(max_step_dir));
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif
