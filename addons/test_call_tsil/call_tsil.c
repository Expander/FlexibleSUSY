#include "call_tsil.h"
#include "tsil.h"

double call_A(double m, double q)
{
   TSIL_REAL m2 = m*m;
   TSIL_REAL q2 = q*q;

   TSIL_COMPLEX result = TSIL_A(m2, q2);

   return TSIL_CREAL(result);
}
