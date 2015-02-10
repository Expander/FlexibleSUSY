
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

struct Parabola {
   static unsigned number_of_calls;
   static double func(const gsl_vector* x, void*) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      number_of_calls++;
      return Sqr(y - 5.0) + Sqr(z - 1.0);
   }
};

unsigned Parabola::number_of_calls = 0;

BOOST_AUTO_TEST_CASE( test_number_of_calls )
{
   const double precision = 1.0e-5;
   Minimizer<2> minimizer(Parabola::func, NULL, 100, precision);
   const double start[2] = { 10, 10 };
   int status = GSL_SUCCESS;

   const gsl_multimin_fminimizer_type* solvers[] =
      {gsl_multimin_fminimizer_nmsimplex, gsl_multimin_fminimizer_nmsimplex2,
       gsl_multimin_fminimizer_nmsimplex2rand};

   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i) {
      Parabola::number_of_calls = 0;
      minimizer.set_solver_type(solvers[i]);
      status = minimizer.minimize(start);

      BOOST_REQUIRE(status == GSL_SUCCESS);
      BOOST_CHECK_CLOSE_FRACTION(minimizer.get_minimum_point(0), 5.0, precision);
      BOOST_CHECK_CLOSE_FRACTION(minimizer.get_minimum_point(1), 1.0, precision);
      BOOST_MESSAGE("solver type " << i << " used " << Parabola::number_of_calls << " calls");
   }
}
