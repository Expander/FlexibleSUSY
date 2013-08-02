
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_root_finder

#include <boost/test/unit_test.hpp>

#include "root_finder.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   struct Parabola {
      static int func(const gsl_vector* x, void*, gsl_vector* f) {
         const double y = gsl_vector_get(x, 0);
         const double z = gsl_vector_get(x, 1);
         gsl_vector_set(f, 0, y*(y - 5.0));
         gsl_vector_set(f, 1, z*(z - 1.0));
         return GSL_SUCCESS;
      }
   };

   Root_finder<2> root_finder(Parabola::func, NULL, 100, 1.0e-5);
   const double start[2] = { 10, 10 };
   const int status = root_finder.find_root(start);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_CHECK_CLOSE_FRACTION(root_finder.get_root(0), 5.0, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(root_finder.get_root(1), 1.0, 1.0e-5);
}
