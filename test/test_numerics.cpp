
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_numerics

#include <boost/test/unit_test.hpp>
#include "numerics.h"

BOOST_AUTO_TEST_CASE( test_bIntegral )
{
   struct Parameters {
      int n;
      double p, m1, m2, mt;
   };

   Parameters parameters[] = {
      {1, 100., 0.  , 0.  , 100.},
      {1, 100., 0.  , 0.  , 100.},
      {1, 100., 0.  , 100., 100.},
      {1, 100., 100., 0.  , 100.},
      {1, 100., 0.  , 100., 100.},
      {1, 100., 100., 100., 100.},
      {1, 100., 10. , 20. , 100.},
      {1, 100., 10. , 200., 100.},
      {1, 100., 10. , 200., 200.}
   };

   for (unsigned i = 0; i < sizeof(parameters)/sizeof(parameters[0]); i++) {
      const double n  = parameters[i].n;
      const double p  = parameters[i].p;
      const double m1 = parameters[i].m1;
      const double m2 = parameters[i].m2;
      const double mt = parameters[i].mt;

      const double b_ss = bIntegral(n, p, m1, m2, mt);
      const double b_fs = bIntegral_threadsave(n, p, m1, m2, mt);

      BOOST_CHECK_EQUAL(b_ss, b_fs);
   }
}
