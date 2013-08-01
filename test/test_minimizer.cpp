
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_minimizer

#include <boost/test/unit_test.hpp>

#include "minimizer.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   struct Parabola {
      static double func(const gsl_vector* x, void*) {
         const double y = gsl_vector_get(x, 0);
         const double z = gsl_vector_get(x, 1);
         return Sqr(y - 5.0) + Sqr(z - 1.0);
      }
   };

   Minimizer<2> minimizer(Parabola::func, NULL, 100, 1.0e-5);
   const double start[2] = { 10, 10 };
   const int status = minimizer.minimize(start);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_CHECK_SMALL(minimizer.get_minimum_value(), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(minimizer.get_minimum_point(0), 5.0, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(minimizer.get_minimum_point(1), 1.0, 1.0e-5);
}
