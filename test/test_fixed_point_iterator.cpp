
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_fixed_point_iterator

#include <boost/test/unit_test.hpp>

#define ENABLE_VERBOSE 1
#define ENABLE_DEBUG 1

#include "fixed_point_iterator.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

class Parabola {
public:
   static void reset() { number_of_calls = 0; }
   static unsigned get_number_of_calls() { return number_of_calls; }

   /**
    * Finding root of f(x,y) = ((x-5)^2, (y-1)^2) ,
    *
    * => Update steps
    *
    * (x,y) = (-25/(x-10), -1/(y-2))
    *
    * @param x touple (x,y)
    *
    * @return fixed point iteration update steps
    */
   static int func(const gsl_vector* x, void*, gsl_vector* f) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      gsl_vector_set(f, 0, -25./(y - 10.));
      gsl_vector_set(f, 1, -1./(z - 2.));
      number_of_calls++;
      return GSL_SUCCESS;
   }

private:
   static unsigned number_of_calls;
};

unsigned Parabola::number_of_calls = 0;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   const double precision = 1.0e-4;
   const double start[2] = { 9, 9 };
   Fixed_point_iterator<2> fpi(Parabola::func, NULL, 1000, precision);
   int status = GSL_SUCCESS;

   Parabola::reset();

   status = fpi.find_fixed_point(start);

   const double residual_1 = MaxRelDiff(5.0, fpi.get_fixed_point(0));
   const double residual_2 = MaxRelDiff(1.0, fpi.get_fixed_point(1));

   // Note: The convergence criterion
   // MaxRelDiff(x_{n+1}, x_{n}) < precision
   // is not very good: The method converges slowly.  This means
   // subsequent steps are very close to each other, but x_n might not
   // be close to the true fixed point.

   BOOST_REQUIRE(status == GSL_SUCCESS);
   BOOST_CHECK_LT(residual_1, 100*precision);
   BOOST_CHECK_LT(residual_2, 100*precision);
   BOOST_MESSAGE("fixed point iterator used " << Parabola::get_number_of_calls() << " calls");
}
