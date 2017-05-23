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
 * @file rkf_integrator.hpp
 * @brief Integration of ODEs using the Runge-Kutta-Fehlberg method
 */

#ifndef RKF_INTEGRATOR_H
#define RKF_INTEGRATOR_H

#include "error.hpp"
#include <string>

#include <Eigen/Core>

namespace flexiblesusy {

namespace runge_kutta {

class RKF_integrator {
public:
   using Derivs = std::function<Eigen::ArrayXd(double, const Eigen::ArrayXd&)>;
   void operator()(double, double, Eigen::ArrayXd&, Derivs, double) const;
private:
   class DisabledOdeintError : Error {
   public:
      explicit DisabledOdeintError(const std::string& msg_) : msg(msg_) {}
      virtual ~DisabledOdeintError() = default;
      virtual std::string what() const override { return msg; }
   private:
      std::string msg;
   };

   struct RKF_observer {
      void operator()(const Eigen::ArrayXd&, double) const;
   };
};

} // namespace runge_kutta

} // namespace flexiblesusy

#endif
