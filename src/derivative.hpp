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
 * @file derivative.hpp
 * @brief contains functions to calculate derivatives numerically
 */

#ifndef DERIVATIVE_H
#define DERIVATIVE_H

#include <functional>
#include <cmath>
#include <limits>

namespace flexiblesusy {

/**
 * Calculates the forward derivative of \f$f(x)\f$ as \f[f'(x) =
 * \frac{f(x+h)-f(x)}{h} + O(h)\f].  This function calls \f$f\f$ once
 * at \f$(x+h)\f$.  The value of \f$f(x)\f$ has to be provided by the
 * caller.  The step size \f$h\f$ is calculated to be \f$h =
 * \sqrt{\epsilon} x\f$ for \f$x\neq 0\f$.  For \f$x=0\f$, the step
 * size is set to \f$h = \epsilon\f$.
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param fx value of \f$f(x)\f$
 * @param eps \f$\epsilon\f$
 *
 * @return derivative
 */
template <class F, class A>
auto derivative_forward_fx(F&& f, A x, decltype(f(x)) fx, A eps = std::numeric_limits<A>::epsilon())
   -> decltype(f(x))
{
   const A h = std::fabs(x) < eps ? eps : std::sqrt(eps) * x;
   volatile const A xph = x + h; // avoid away optimization
   const A dx = xph - x;
   return (f(xph) - fx) / dx;
}

/**
 * Calculates the forward derivative of \f$f(x)\f$ as \f[f'(x) =
 * \frac{f(x+h)-f(x)}{h} + O(h)\f].  This function calculates
 * \f$f(x)\f$ once and calls derivative_forward_fx().
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param eps measure for step size \f$h\f$
 *
 * @return derivative
 */
template <class F, class A>
auto derivative_forward(F&& f, A x, A eps = std::numeric_limits<A>::epsilon()) -> decltype(f(x))
{
   return derivative_forward_fx(std::forward<F>(f), x, f(x), eps);
}

/**
 * Calculates the backward derivative of \f$f(x)\f$ as \f[f'(x) =
 * \frac{f(x)-f(x-h)}{h} + O(h)\f].  This function calls \f$f\f$ once
 * at \f$(x+h)\f$.  The value of \f$f(x)\f$ has to be provided by the
 * caller.  The step size \f$h\f$ is calculated to be \f$h =
 * \sqrt{\epsilon} x\f$ for \f$x\neq 0\f$.  For \f$x=0\f$, the step
 * size is set to \f$h = \epsilon\f$.
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param fx value of \f$f(x)\f$
 * @param eps \f$\epsilon\f$
 *
 * @return derivative
 */
template <class F, class A>
auto derivative_backward_fx(F&& f, A x, decltype(f(x)) fx, A eps = std::numeric_limits<A>::epsilon())
   -> decltype(f(x))
{
   const A h = std::fabs(x) < eps ? eps : std::sqrt(eps) * x;
   volatile const A xph = x - h; // avoid away optimization
   const A dx = x - xph;
   return (fx - f(xph)) / dx;
}

/**
 * Calculates the backward derivative of \f$f(x)\f$ as \f[f'(x) =
 * \frac{f(x)-f(x-h)}{h} + O(h)\f].  This function calculates
 * \f$f(x)\f$ once and calls derivative_backward_fx().
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param eps measure for step size \f$h\f$
 *
 * @return derivative
 */
template <class F, class A>
auto derivative_backward(F&& f, A x, A eps = std::numeric_limits<A>::epsilon()) -> decltype(f(x))
{
   return derivative_backward_fx(std::forward<F>(f), x, f(x), eps);
}

/**
 * Calculates the central derivative of \f$f(x)\f$ as \f[f'(x) =
 * \frac{f(x+h)-f(x-h)}{2h} + O(h^2)\f].  This function calls \f$f\f$
 * twice.
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param eps measure for step size \f$h\f$
 *
 * @return derivative
 */
template <class F, class A>
auto derivative_central(F&& f, A x, A eps = std::numeric_limits<A>::epsilon()) -> decltype(f(x))
{
   const A h = std::fabs(x) < eps ? eps : std::sqrt(eps) * x;
   const A f_x_plus_h = f(x + h);
   const A f_x_minus_h = f(x - h);

   return (f_x_plus_h - f_x_minus_h) / (2*h);
}

/**
 * Calculates the derivative of \f$f(x)\f$ using the five-point
 * stencil method as \f[f'(x) = \frac{-f(x+2h) + 8f(x+h) - 8f(x-h) +
 * f(x-2h)}{12h} + O(h^4)\f].
 *
 * This function calls \f$f\f$ four times.  The step size \f$h\f$ is
 * calculated to be \f$h = \sqrt{\epsilon} x\f$ for \f$x\neq 0\f$.
 * For \f$x=0\f$, the step size is set to \f$h = \epsilon\f$.
 *
 * @param f function
 * @param x point at which derivative is to be calculated
 * @param eps measure for step size \f$h\f$
 *
 * @return derivative
 */
template <class F, class A>
auto derivative_five_point_stencil(F&& f, A x, A eps = std::numeric_limits<A>::epsilon())
   -> decltype(f(x))
{
   const A h = std::fabs(x) < eps ? eps : std::sqrt(eps) * x;
   const A f_x_plus_2h = f(x + 2*h);
   const A f_x_plus_h = f(x + h);
   const A f_x_minus_h = f(x - h);
   const A f_x_minus_2h = f(x - 2*h);

   return (-f_x_plus_2h + 8*f_x_plus_h - 8*f_x_minus_h + f_x_minus_2h) / (12*h);
}

} // namespace flexiblesusy

#endif
