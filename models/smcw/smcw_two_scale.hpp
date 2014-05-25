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

#ifndef SMCW_TWO_SCALE_H
#define SMCW_TWO_SCALE_H

#include "smcw.hpp"
#include "sm_two_scale.hpp"

#include <iostream>

namespace flexiblesusy {

class Two_scale;

/// Formatted output
std::ostream& operator <<(std::ostream&, const StandardModelCW<Two_scale>&);

template<>
class StandardModelCW<Two_scale>: public StandardModel<Two_scale> {
public:
   const static int numStandardModelCWPars =
      StandardModel<Two_scale>::numStandardModelPars + 3;

   StandardModelCW();
   StandardModelCW(const StandardModelCW<Two_scale>&);
   StandardModelCW(const StandardModel<Two_scale>&, double, double, double);
   StandardModelCW(const DoubleMatrix& yu, const DoubleMatrix& yd,
                   const DoubleMatrix& ye, const DoubleVector& g,
                   double lambda, double vs);

   virtual ~StandardModelCW();

   virtual void calculate_spectrum() {}
   virtual std::string name() const { return "SMCW"; }
   virtual void print(std::ostream& s) const { s << *this; }

   /// sets object to be equal to another
   const StandardModelCW<Two_scale> & operator=(const StandardModelCW<Two_scale>& s);

   /// Set a single gauge coupling
   virtual void setGaugeCoupling(int, double);
   /// Set all gauge couplings
   virtual void setAllGauge(const DoubleVector&);
   /// Set lambda parameter
   void setLambda(double l) { lambda = l; }
   /// Set s VEV
   void setVs(double v) { vs = v; }

   /// Returns a single gauge coupling
   virtual double displayGaugeCoupling(int) const;
   /// Returns all gauge couplings
   virtual DoubleVector displayGauge() const;
   /// Returns lambda parameter
   double displayLambda() const { return lambda; }
   /// Returns s VEV
   double displayVs() const { return vs; }
   /// Calculate beta functions
   StandardModelCW calc_beta() const;

   /// Calculate Z' mass
   double calcZprimeMass() const;

protected:
   /// Sets all RGE parameters to elements of vector
   virtual void set(const DoubleVector&);
   /// Returns all parameters as elements of a vector
   virtual const DoubleVector display() const;
   /// Calculate beta functions
   virtual DoubleVector beta() const;

private:
   double g4;     ///< g4 gauge coupling
   double lambda; ///< scalar potential parameter lambda
   double vs;     ///< VEV of singlet field
};

}

#endif
