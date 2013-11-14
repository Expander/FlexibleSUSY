
#include "wrappers.hpp"

#include <Eigen/SVD>

namespace flexiblesusy {

double MaxRelDiff(double a, double b)
{
   const double sTin = fabs(a), sTout = fabs(b);
   const double maxx = std::max(sTin, sTout);
   const double underflow = 1.0e-20;

   if (maxx < underflow)
      return 0.0;

   return fabs(1.0 - std::min(sTin, sTout) / maxx);
}

}
