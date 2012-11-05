
#ifndef SM_TWO_SCALE_H
#define SM_TWO_SCALE_H

#include "rge.h"
#include "sm.hpp"

#include <iostream>

class Two_scale;

template<>
class StandardModel<Two_scale>: public RGE {
private:
   DoubleMatrix yu, yd, ye; ///< Yukawa matrices for ups, downs and leptons
   DoubleVector g;          ///< Gauge couplings (g1 = sqrt(5/3) g_Y)

public:
   const static int numStandardModelPars = 3 * 3 * 3 + 3;
   typedef enum {YU = 1, YD, YE} yukawa;

   StandardModel();
   StandardModel(const StandardModel&);
   StandardModel(const DoubleMatrix& yu, const DoubleMatrix& yd,
                 const DoubleMatrix& ye, const DoubleVector& g);

   virtual ~StandardModel();

   /// sets object to be equal to another
   const StandardModel & operator=(const StandardModel& s);

   /// Sets Yukawa matrix element
   void setYukawaElement(yukawa, int, int, double);
   /// Sets whole Yukawa matrix
   void setYukawaMatrix(yukawa, const DoubleMatrix&);
   /// Set a single gauge coupling
   void setGaugeCoupling(int, double);
   /// Set all gauge couplings
   void setAllGauge(const DoubleVector&);
   /// Sets all RGE parameters to elements of vector
   void set(const DoubleVector&);

   /// Returns a single Yukawa matrix element
   double displayYukawaElement(yukawa, int, int) const;
   /// Returns a whole Yukawa matrix
   DoubleMatrix displayYukawaMatrix(yukawa) const;
   /// Returns a single gauge coupling
   double displayGaugeCoupling(int) const;
   /// Returns all gauge couplings
   DoubleVector displayGauge() const;
   /// Returns all parameters as elements of a vector
   const DoubleVector display() const;
   /// Calculate beta functions
   DoubleVector beta() const;
   /// Calculate beta functions
   StandardModel calcBeta() const;
};

/// Formatted output
std::ostream& operator <<(std::ostream&, const StandardModel<Two_scale>&);

#endif
