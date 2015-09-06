#include <complex>
#include "wrappers.hpp"
#include "rge.h"

using namespace flexiblesusy;

int main()
{
   std::complex<double> z(1.0, 2.0);

   // Test for unambiguous operator*
   std::complex<double> result = 2 * z;

   return 0;
}
