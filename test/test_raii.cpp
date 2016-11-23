#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_raii

#include <boost/test/unit_test.hpp>

#include "raii.hpp"

using namespace flexiblesusy;

struct Real_model {
   double p{};
   Real_model(double p_) : p(p_) {}
};

BOOST_AUTO_TEST_CASE( test_raii_save )
{
   Real_model m(1.);

   BOOST_CHECK_EQUAL(m.p, 1.);

   {
      RAII_save<double> save(m.p);
      m.p = 2.;
      BOOST_CHECK_EQUAL(m.p, 2.);
   }

   BOOST_CHECK_EQUAL(m.p, 1.);
}

BOOST_AUTO_TEST_CASE( test_make_raii_save )
{
   Real_model m(1.);

   BOOST_CHECK_EQUAL(m.p, 1.);

   {
      const auto save = make_raii_save(m.p);
      m.p = 2.;
      BOOST_CHECK_EQUAL(m.p, 2.);
   }

   BOOST_CHECK_EQUAL(m.p, 1.);
}
