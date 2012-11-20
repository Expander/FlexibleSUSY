#include "sm_two_scale.hpp"
#include "sm_two_scale_experimental_constraint.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_two_scale

#include <boost/test/unit_test.hpp>

#define YU StandardModel<Two_scale>::YU
#define YD StandardModel<Two_scale>::YD
#define YE StandardModel<Two_scale>::YE

BOOST_AUTO_TEST_CASE( test_electroweak_constaint )
{
   StandardModel<Two_scale> sm;
   sm.setScale(ewConstants::MZ);
   StandardModelExpConstraint sm_ew_constraint(&sm);
   const DoubleMatrix zero3x3(3,3);

   // check that sm is initialized to zero
   for (int i = 1; i <= 3; ++i)
      BOOST_CHECK_EQUAL(sm.displayGaugeCoupling(i), 0.0);

   BOOST_CHECK_EQUAL(sm.displayYukawaMatrix(YU), zero3x3);
   BOOST_CHECK_EQUAL(sm.displayYukawaMatrix(YD), zero3x3);
   BOOST_CHECK_EQUAL(sm.displayYukawaMatrix(YE), zero3x3);

   sm_ew_constraint.apply();

   // check that ew constraints have been applied
   BOOST_CHECK_EQUAL(sm.displayGaugeCoupling(1), ewConstants::g1);
   BOOST_CHECK_EQUAL(sm.displayGaugeCoupling(2), ewConstants::g2);
   BOOST_CHECK_EQUAL(sm.displayGaugeCoupling(3), ewConstants::g3);

   BOOST_CHECK_EQUAL(sm.displayYukawaElement(YU, 3, 3), ewConstants::yt);
   BOOST_CHECK_EQUAL(sm.displayYukawaElement(YD, 3, 3), ewConstants::yb);
   BOOST_CHECK_EQUAL(sm.displayYukawaElement(YE, 3, 3), ewConstants::ytau);
}
