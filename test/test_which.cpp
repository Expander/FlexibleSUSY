#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_which

#include <boost/test/unit_test.hpp>
#include "which.hpp"

using namespace flexiblesusy;

namespace {

int eval_count = 0;

}

int fail(bool really)
{
   eval_count++;
   if (really)
      BOOST_ERROR(eval_count << ": eager evaluation occurred!");
   else
      BOOST_TEST_MESSAGE(eval_count << ": test evaluation");
   return 0;
}

BOOST_AUTO_TEST_CASE(test_calls_to_fail)
{
   eval_count = 0;
   fail(false);
   fail(false);
   BOOST_CHECK(eval_count == 2);
}

BOOST_AUTO_TEST_CASE(test_WHICH_1)
{
   BOOST_CHECK_EQUAL(WHICH(true, 0), 0);
}

BOOST_AUTO_TEST_CASE(test_WHICH_2)
{
   const int i = 1;

   BOOST_CHECK_EQUAL(WHICH(i == 1, 1), 1);
   BOOST_CHECK_EQUAL(WHICH(i == 0, (fail(true), 1)), 0);
}

BOOST_AUTO_TEST_CASE(test_WHICH_4)
{
   const int i = 0;

   BOOST_CHECK_EQUAL(WHICH(i > 0, (fail(true), 0), true, 1), 1);
   BOOST_CHECK_EQUAL(WHICH(i < 1, 0, true, (fail(true), 1)), 0);
}

BOOST_AUTO_TEST_CASE(test_WHICH_type)
{
   BOOST_CHECK_EQUAL(WHICH(false, 0, true, 0.5), 0.5);
   BOOST_CHECK_EQUAL(WHICH(false, 0, false, 0, true, 0.5), 0.5);
}

BOOST_AUTO_TEST_CASE(test_WHICH_10)
{
   const int i = 0;

   BOOST_CHECK_EQUAL(WHICH(i > 0, (fail(true), 0),
                           i > 1, (fail(true), 1),
                           i > 2, (fail(true), 2),
                           i > 3, (fail(true), 3),
                           i > 4, (fail(true), 4),
                           i > 5, (fail(true), 5),
                           i > 6, (fail(true), 6),
                           i > 7, (fail(true), 7),
                           i > 8, (fail(true), 8),
                           i > 9, (fail(true), 9),
                           true, -1), -1);
}
