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

#ifndef BETAFUNCTION_H
#define BETAFUNCTION_H

#include "rk.hpp"

class Beta_function {
public:
   Beta_function();
   virtual ~Beta_function() {}

   void set_scale(double s) { scale = s; }
   void set_parameters(unsigned pars) { numPars = pars; }
   void set_loops(unsigned l) { loops = l; }

   double get_scale() { return scale; }
   unsigned get_parameters() { return numPars; }
   unsigned get_loops() { return loops; }

   virtual const Eigen::ArrayXd display() const = 0;
   virtual void set(const Eigen::ArrayXd&) = 0;
   virtual Eigen::ArrayXd beta() const = 0;

   virtual int run(double, double, double eps = -1.0);
   virtual int run_to(double, double eps = -1.0);

private:
   unsigned numPars;     ///< Number of parameters
   unsigned loops;       ///< To what order does the RG evolution run
   double scale;         ///< Renormalisation scale
   double tolerance;     ///< running tolerance
   double min_tolerance; ///< minimum tolerance allowed

   int call_rk(double, double, Eigen::ArrayXd&,
               runge_kutta::Derivs, double eps = -1.0);
   Eigen::ArrayXd derivatives(double, const Eigen::ArrayXd&);
   double get_tolerance(double eps);
};

#endif
