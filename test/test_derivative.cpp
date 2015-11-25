#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_derivative

#include <boost/test/unit_test.hpp>
#include <iomanip>

#include "derivative.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_derivative_x_sqr )
{
   auto f = [](double x) { return x*x; };
   auto df = [](double x) { return 2*x; };

   std::vector<double> x_values = {0., 1., 2., 3., 1e1, 1e2, 1e8};

   for (auto x: x_values) {
      if (x == 0.) {
         BOOST_CHECK_SMALL(derivative_forward(f,x), 1e-8);
         BOOST_CHECK_SMALL(derivative_backward(f,x), 1e-8);
         BOOST_CHECK_SMALL(derivative_central<0>(f,x), 2e-9);
         BOOST_CHECK_SMALL(derivative_central<1>(f,x), 2e-9);
         BOOST_CHECK_SMALL(derivative_central<2>(f,x), 2e-9);
         BOOST_CHECK_SMALL(derivative_central<3>(f,x), 2e-9);
      } else {
         BOOST_CHECK_CLOSE_FRACTION(derivative_forward(f,x), df(x), 1e-8);
         BOOST_CHECK_CLOSE_FRACTION(derivative_backward(f,x), df(x), 1e-8);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<0>(f,x), df(x), 2e-9);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<1>(f,x), df(x), 2e-9);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<2>(f,x), df(x), 2e-9);
         BOOST_CHECK_CLOSE_FRACTION(derivative_central<3>(f,x), df(x), 2e-9);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_stability )
{
   auto f = [](double x) { return x*x; };
   auto df = [](double x) { return 2*x; };

   const double x = 1.;

   std::vector<double> eps_values =
      {1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2};

   for (auto eps: eps_values) {
      const double df_true = df(x);
      const double df_forw = derivative_forward(f,x,eps);
      const double df_back = derivative_backward(f,x,eps);
      const double df_cent = derivative_central<0>(f,x,eps);
      const double df_five = derivative_central<1>(f,x,eps);

      BOOST_MESSAGE(std::setprecision(16) << std::scientific
                    << "eps = " << eps << ": true df = " << df_true
                    << ", forward df = " << df_forw
                    << ", backward df = " << df_back
                    << ", central df = " << df_cent
                    << ", stencil df = " << df_five);

      BOOST_CHECK_CLOSE_FRACTION(df_forw, df_true, 6e-2);
      BOOST_CHECK_CLOSE_FRACTION(df_back, df_true, 6e-2);
      BOOST_CHECK_CLOSE_FRACTION(df_cent, df_true, 2e-8);
      BOOST_CHECK_CLOSE_FRACTION(df_five, df_true, 2e-8);
   }
}
