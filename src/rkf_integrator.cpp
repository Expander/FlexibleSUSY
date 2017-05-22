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
 * @file rkf_integrator.cpp
 * @brief Implementation of the RKF_integrator class
 */

#include "rkf_integrator.hpp"
#include "config.h"

#ifdef ENABLE_ODEINT

#include <boost/numeric/odeint/algebra/algebra_dispatcher.hpp>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>

namespace boost {
namespace numeric {
namespace odeint {

template<>
struct vector_space_norm_inf<Eigen::ArrayXd> {
   using result_type = double;
   double operator()(const Eigen::ArrayXd& x) const
      {
         return x.cwiseAbs().maxCoeff();
      }
};

} // namespace odeint
} // namespace numeric
} // namespace boost

namespace flexiblesusy {

namespace runge_kutta {

void RKF_integrator::operator()(double start, double end,
                                Eigen::ArrayXd& pars, Derivs derivs,
                                double) const
{
   using state_type = Eigen::ArrayXd;
   using stepper_type = boost::numeric::odeint::runge_kutta_fehlberg78<
      state_type, double, state_type, double,
      boost::numeric::odeint::vector_space_algebra
      >;

   const double guess = (end - start) * 0.1; // first step size
   const auto derivatives = [derivs] (const state_type& y, state_type& dydt, double t) -> void {
      dydt = derivs(t, y);
   };

   stepper_type stepper;
   boost::numeric::odeint::integrate_adaptive(
      stepper, derivatives, pars, start, end, guess);
}

} // namespace runge_kutta

} // namespace flexiblesusy

#else

namespace flexiblesusy {

namespace runge_kutta {

void RKF_integrator::operator()(double, double, Eigen::ArrayXd&, Derivs,
                                double) const
{
   throw DisabledOdeintError("Cannot call operator(), because odeint support is disabled.");
}

} // namespace runge_kutta

} // namespace flexiblesusy

#endif
