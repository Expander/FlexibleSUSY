
#include "two_scale_running_precision.hpp"

#include <cmath>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_running_precision

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_constant_running_precision )
{
   for (int i = 0; i < 5; ++i) {
      const double default_precison = static_cast<double>(i);
      Two_scale_constant_precision const_prec(default_precison);

      // test that all iterations yield the same precison
      for (int iteration = 0; iteration < 10; ++iteration) {
         double precision = const_prec.get_precision(iteration);
         BOOST_CHECK_EQUAL(precision, default_precison);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_increasing_running_precision )
{
   for (int i = 1; i < 10; ++i) {
      const double increasing_factor = static_cast<double>(i);
      const double min_precision = 1.0e-12;
      Two_scale_increasing_precision incr_prec(increasing_factor, min_precision);

      // test that all iterations yield expected precision
      for (int iteration = 0; iteration < 20; ++iteration) {
         double precision = incr_prec.get_precision(iteration);
         double expected_precision = std::pow(1./increasing_factor, iteration + 1);
         if (expected_precision < min_precision)
            expected_precision = min_precision;
         BOOST_CHECK_CLOSE(precision, expected_precision, 1.0e-5);
      }
   }
}
