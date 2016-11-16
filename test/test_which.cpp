#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_which

#include <boost/test/unit_test.hpp>
#include "which.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE(test_WHICH_1)
{
   BOOST_CHECK_EQUAL(WHICH(true, 0), 0);
}

BOOST_AUTO_TEST_CASE(test_WHICH_2)
{
   const int i = 1;

   BOOST_CHECK_EQUAL(WHICH(i == 1, 1), 1);
   BOOST_CHECK_EQUAL(WHICH(i == 0, 1), 0);
}

BOOST_AUTO_TEST_CASE(test_WHICH_4)
{
   const int i = 0;

   BOOST_CHECK_EQUAL(WHICH(i > 0, 0, true, 1), 1);
}

BOOST_AUTO_TEST_CASE(test_WHICH_10)
{
   const int i = 0;

   BOOST_CHECK_EQUAL(WHICH(i > 0, 0,
                           i > 1, 1,
                           i > 2, 2,
                           i > 3, 3,
                           i > 4, 4,
                           i > 5, 5,
                           i > 6, 6,
                           i > 7, 7,
                           i > 8, 8,
                           i > 9, 9,
                           true, -1), -1);
}
