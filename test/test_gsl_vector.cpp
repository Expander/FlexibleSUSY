#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_gsl_vector

#include <boost/test/unit_test.hpp>

#include "gsl_vector.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_init_default )
{
   GSL_vector v;
   BOOST_CHECK_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_init_0 )
{
   GSL_vector v(0);
   BOOST_CHECK_EQUAL(v.size(), 0);
}

BOOST_AUTO_TEST_CASE( test_init )
{
   GSL_vector v(3);
   BOOST_CHECK_EQUAL(v[0], 0.);
   BOOST_CHECK_EQUAL(v[1], 0.);
   BOOST_CHECK_EQUAL(v[2], 0.);
}

BOOST_AUTO_TEST_CASE( test_operator )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   BOOST_CHECK_EQUAL(v[0], 1.);
   BOOST_CHECK_EQUAL(v[1], 2.);
   BOOST_CHECK_EQUAL(v[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_copy )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2(v);

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_assign )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2;
   v2 = v;

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_assign_chain )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2, v3;
   v3 = v2 = v;

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
   BOOST_CHECK_EQUAL(v3[0], 1.);
   BOOST_CHECK_EQUAL(v3[1], 2.);
   BOOST_CHECK_EQUAL(v3[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_move )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   GSL_vector v2(std::move(v));

   BOOST_CHECK_EQUAL(v2[0], 1.);
   BOOST_CHECK_EQUAL(v2[1], 2.);
   BOOST_CHECK_EQUAL(v2[2], 3.);
}

BOOST_AUTO_TEST_CASE( test_output )
{
   GSL_vector v(3);
   v[0] = 1.;
   v[1] = 2.;
   v[2] = 3.;

   BOOST_MESSAGE(v);
}
