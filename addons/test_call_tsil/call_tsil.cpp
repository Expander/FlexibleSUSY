#include "call_tsil.hpp"
#include "tsil_cpp.h"

double call_A(double m, double q)
{
   TSIL_REAL m2 = m*m;
   TSIL_REAL q2 = q*q;

   TSIL_COMPLEXCPP result = TSIL_A_(m2, q2);

   return std::real(result);
}
