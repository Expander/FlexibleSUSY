
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_numerics

#include <boost/test/unit_test.hpp>
#include "numerics.h"
#include "numerics2.hpp"

using namespace flexiblesusy;

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
      const double b_fs = b1_general(p, m1, m2, mt);

      BOOST_CHECK_CLOSE(b_ss, b_fs, 3e-3);
   }
}

BOOST_AUTO_TEST_CASE(test_is_equal)
{
   short s = (short)0;
   int i = 0;
   long l = 0L;
   long long ll = 0LL;
   unsigned short us = 0;
   unsigned int ui = 0;
   unsigned long ul = 0UL;
   unsigned long long ull = 0ULL;

   BOOST_CHECK(is_zero(s));
   BOOST_CHECK(is_zero(i));
   BOOST_CHECK(is_zero(l));
   BOOST_CHECK(is_zero(ll));
   BOOST_CHECK(is_zero(us));
   BOOST_CHECK(is_zero(ui));
   BOOST_CHECK(is_zero(ul));
   BOOST_CHECK(is_zero(ull));

   // BOOST_CHECK(is_equal(s, (short)0));
   BOOST_CHECK(is_equal(i, 0));
   BOOST_CHECK(is_equal(l, 0L));
   BOOST_CHECK(is_equal(ll, 0LL));
   // BOOST_CHECK(is_equal(us, (unsigned short)0));
   BOOST_CHECK(is_equal(ui, 0U));
   BOOST_CHECK(is_equal(ul, 0UL));
   BOOST_CHECK(is_equal(ull, 0ULL));
}
