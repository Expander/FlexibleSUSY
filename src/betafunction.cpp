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

#include "betafunction.hpp"
#include "logger.hpp"
#include <cmath>

namespace flexiblesusy {

Beta_function::Beta_function()
   : numPars(0)
   , loops(0)
   , thresholds(0)
   , scale(0.0)
   , tolerance(1.e-4)
   , min_tolerance(1.0e-11)
{
}

int Beta_function::run_to(double x2, double eps)
{
   const double tol = get_tolerance(eps);
   return run(scale, x2, tol);
}

int Beta_function::run(double x1, double x2, double eps)
{
   using namespace std::placeholders;

   const double tol = get_tolerance(eps);

   if (std::fabs(x1) < tol || std::fabs(x2) < tol) {
      ERROR("Beta_function::run: One of the scales is close to zero:"
            " starting scale = " << x1 << ", end scale = " << x2);
      return 1;
   }

   Eigen::ArrayXd y(display());
   runge_kutta::Derivs derivs = std::bind(&Beta_function::derivatives,
                                          this, _1, _2);

   int err = call_rk(x1, x2, y, derivs, tol);

   if (err == 0) {
      set(y);
      set_scale(x2);
   }

   return err;
}

Eigen::ArrayXd Beta_function::derivatives(double x, const Eigen::ArrayXd& y)
{
   set_scale(exp(x));
   set(y);
   return beta();
}

int Beta_function::call_rk(double x1, double x2, Eigen::ArrayXd & v,
                           runge_kutta::Derivs derivs, double eps)
{
   if (fabs(x1 - x2) < min_tolerance)
      return 0;

   const double tol = get_tolerance(eps);
   const double from = log(fabs(x1));
   const double to = log(fabs(x2));
   const double guess = (from - to) * 0.1; //first step size
   const double hmin = (from - to) * tol * 1.0e-5;

   const int err =
      runge_kutta::integrateOdes(v, from, to, tol, guess, hmin, derivs,
                                 runge_kutta::odeStepper);

   set_scale(x2);

   return err;
}

double Beta_function::get_tolerance(double eps)
{
   double tol;
   if (eps < 0.0)
      tol = tolerance;
   else if (eps < min_tolerance)
      tol = min_tolerance;
   else
      tol = eps;

   return tol;
}

}
