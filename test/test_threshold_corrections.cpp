#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_threshold_corrections

#include <boost/test/unit_test.hpp>
#include "threshold_corrections.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_get )
{
   Threshold_corrections tc;

   tc.alpha_em    = 0;
   tc.sin_theta_w = 1;
   tc.alpha_s     = 2;
   tc.mw          = 3;
   tc.mz          = 4;
   tc.mt          = 5;
   tc.mb          = 6;
   tc.me          = 7;
   tc.mm          = 8;
   tc.mtau        = 9;

   BOOST_CHECK_EQUAL(tc.get(), 9876543210);
}

BOOST_AUTO_TEST_CASE( test_set )
{
   Threshold_corrections tc;
   tc.set(9876543210);

   BOOST_CHECK_EQUAL(tc.alpha_em    , 0);
   BOOST_CHECK_EQUAL(tc.sin_theta_w , 1);
   BOOST_CHECK_EQUAL(tc.alpha_s     , 2);
   BOOST_CHECK_EQUAL(tc.mw          , 3);
   BOOST_CHECK_EQUAL(tc.mz          , 4);
   BOOST_CHECK_EQUAL(tc.mt          , 5);
   BOOST_CHECK_EQUAL(tc.mb          , 6);
   BOOST_CHECK_EQUAL(tc.me          , 7);
   BOOST_CHECK_EQUAL(tc.mm          , 8);
   BOOST_CHECK_EQUAL(tc.mtau        , 9);
}
