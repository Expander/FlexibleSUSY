
#include "wrappers.hpp"

#include <Eigen/SVD>

namespace flexiblesusy {

Eigen::Matrix3d Diag(const Eigen::Matrix3d& m)
{
   Eigen::Matrix3d diag(m);
   for (int i = 0; i < 3; ++i) {
      for (int k = 0; k < 3; ++k) {
         if (i != k)
            diag(i,k) = 0.0;
      }
   }
   return diag;
}

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
