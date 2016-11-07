
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_minimizer

#include <boost/test/unit_test.hpp>

#include "minimizer.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

typedef Eigen::Matrix<double,2,1> EV2_t;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   auto parabola = [](const EV2_t& x) -> double {
      const double y = x(0);
      const double z = x(1);
      return Sqr(y - 5.0) + Sqr(z - 1.0);
   };

   Minimizer<2> minimizer(parabola, 100, 1.0e-5);
   Eigen::Matrix<double,2,1> start;
   start << 10, 10;
   const int status = minimizer.minimize(start);
   const auto minimum_point = minimizer.get_solution();

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_CHECK_SMALL(minimizer.get_minimum_value(), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(minimum_point(0), 5.0, 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(minimum_point(1), 1.0, 1.0e-5);
}

static unsigned number_of_calls = 0;

auto parabola = [](const EV2_t& x) -> double {
   number_of_calls++;
   const double y = x(0);
   const double z = x(1);
   return Sqr(y - 5.0) + Sqr(z - 1.0);
};

BOOST_AUTO_TEST_CASE( test_number_of_calls )
{
   const double precision = 1.0e-5;
   Minimizer<2> minimizer(parabola, 100, precision);
   Eigen::Matrix<double,2,1> start;
   start << 10, 10;
   int status = GSL_SUCCESS;

   Minimizer<2>::Solver_type solvers[] =
      { Minimizer<2>::GSLSimplex,
        Minimizer<2>::GSLSimplex2,
        Minimizer<2>::GSLSimplex2Rand };

   for (std::size_t i = 0; i < sizeof(solvers)/sizeof(*solvers); ++i) {
      number_of_calls = 0;
      minimizer.set_solver_type(solvers[i]);
      status = minimizer.minimize(start);
      const auto minimum_point = minimizer.get_solution();

      BOOST_REQUIRE(status == GSL_SUCCESS);
      BOOST_CHECK_CLOSE_FRACTION(minimum_point(0), 5.0, precision);
      BOOST_CHECK_CLOSE_FRACTION(minimum_point(1), 1.0, precision);
      BOOST_TEST_MESSAGE("solver type " << i << " used " << number_of_calls << " calls");
   }
}
