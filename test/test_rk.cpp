
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_rk

#include <boost/test/unit_test.hpp>

#include "rk.hpp"
#include "rk_legacy.hpp"
#include "error.hpp"
#include "wrappers.hpp"
#include <Eigen/Dense>

using namespace flexiblesusy;
using namespace softsusy;

typedef DoubleVector (*Derivative_t)(double, const DoubleVector&);

Eigen::ArrayXd ToEigenArray(const DoubleVector& v)
{
   Eigen::ArrayXd a(v.size());
   for (int i = v.displayStart(); i <= v.displayEnd(); i++)
      a(i - 1) = v(i);
   return a;
}

DoubleVector ToDoubleVector(const Eigen::ArrayXd& a)
{
   DoubleVector v(a.rows());
   for (int i = 0; i < a.rows(); i++)
      v(i + 1) = a(i);
   return v;
}

template <typename Derivs>
void check_integrateOdes(DoubleVector& parameters,
                         double start, double end,
                         Derivative_t beta_legacy,
                         Derivs beta_eigen)
{
   Eigen::ArrayXd parameter_eigen(ToEigenArray(parameters));

   const double from = log(start), to = log(end);
   const double tol = 1.0e-7;
   const double guess = (from - to) * 0.1; // first step size
   const double hmin = (from - to) * tol * 1.0e-5; // min step size

   const int err_legacy
      = integrateOdes(parameters, from, to, tol, guess,
                      hmin, beta_legacy, odeStepper);

   int err_eigen = 0;
   try {
      runge_kutta::integrateOdes(parameter_eigen, from, to, tol, guess,
                                 hmin, beta_eigen);
   } catch (Error&) {
      err_eigen = 1;
   }

   BOOST_CHECK_EQUAL(err_legacy, err_eigen);

   for (std::size_t i = 0; i < parameters.size(); i++)
      BOOST_CHECK_CLOSE(parameters(i+1), parameter_eigen(i), 1.0e-10);
}

#define DEFINE_BETA_EIGEN(beta) \
   Eigen::ArrayXd beta ## _eigen(double x, const Eigen::ArrayXd& parameter) \
   {                                                                    \
      const DoubleVector _parameter_legacy(ToDoubleVector(parameter));  \
      const DoubleVector _beta_legacy(beta ## _legacy(x, _parameter_legacy)); \
      const Eigen::ArrayXd _beta(ToEigenArray(_beta_legacy));           \
      return _beta;                                                     \
   }

// ============================== one dimension ==============================

DoubleVector beta_one_dim_legacy(double /* x */, const DoubleVector& parameter)
{
   DoubleVector beta(parameter.size());
   beta(1) = parameter(1)*parameter(1) * 0.1;
   return beta;
}

DEFINE_BETA_EIGEN(beta_one_dim)

BOOST_AUTO_TEST_CASE( test_one_dim )
{
   DoubleVector parameter_legacy(1);
   parameter_legacy(1) = 0.5;

   check_integrateOdes(parameter_legacy, 100., 1.0e10,
                       beta_one_dim_legacy, beta_one_dim_eigen);
}

// ============================== 10 dimensions ==============================

DoubleVector beta_ten_dim_legacy(double /* x */, const DoubleVector& parameter)
{
   DoubleVector beta(parameter.size());
   for (std::size_t i = 1; i <= parameter.size(); i++)
      beta(i) = parameter(i)*parameter(1) * 0.1 - 0.5 * i;
   return beta;
}

DEFINE_BETA_EIGEN(beta_ten_dim)

BOOST_AUTO_TEST_CASE( test_ten_dim )
{
   DoubleVector parameter_legacy(10);
   for (std::size_t i = 1; i <= parameter_legacy.size(); i++)
      parameter_legacy(i) = 0.5 + 0.1 * i*i;

   check_integrateOdes(parameter_legacy, 100., 1.0e10,
                       beta_ten_dim_legacy, beta_ten_dim_eigen);
}

// ============================== non-perturbative ==============================

DoubleVector beta_non_pert_legacy(double /* x */, const DoubleVector& parameter)
{
   DoubleVector beta(parameter.size());
   for (std::size_t i = 1; i <= parameter.size(); i++)
      beta(i) = parameter(i)*parameter(i) + 0.5 * i;
   return beta;
}

DEFINE_BETA_EIGEN(beta_non_pert)

BOOST_AUTO_TEST_CASE( test_non_perturbative )
{
   DoubleVector parameter_legacy(100);
   for (std::size_t i = 1; i <= parameter_legacy.size(); i++)
      parameter_legacy(i) = 3.0 + 0.1 * i*i;

   check_integrateOdes(parameter_legacy, 100., 1.0e10,
                       beta_non_pert_legacy, beta_non_pert_eigen);
}

Eigen::ArrayXd beta_gauge_one_loop(double x, const Eigen::ArrayXd& parameters)
{
   const int num_pars = parameters.size();
   Eigen::ArrayXd beta(num_pars);

   for (int i = 0; i < num_pars; ++i) {
      beta(i) = (9.0 / (3. * i - 2.)) * Cube(parameters(i));
   }

   return oneOver16PiSqr * beta;
}

BOOST_AUTO_TEST_CASE( test_numerical_solution_gauge_one_loop )
{
   const int num_pars = 3;
   Eigen::ArrayXd init_parameters(num_pars);
   for (int i = 0; i < num_pars; ++i) {
      init_parameters(i) = 1.2 / Power(2., i);
   }

   Eigen::ArrayXd parameters(init_parameters);
   const double start = 100.;
   const double end = 1.0e10;
   const double from = log(start);
   const double to = log(end);
   const double tol = 1.0e-7;
   const double guess = (from - to) * 0.1; // first step size
   const double hmin = (from - to) * tol * 1.0e-5; // min step size

   int error = 0;
   try {
      runge_kutta::integrateOdes(parameters, from, to, tol, guess,
                                 hmin, beta_gauge_one_loop);
   } catch (Error&) {
      error = 1;
   }

   BOOST_CHECK_EQUAL(error, 0);

   for (int i = 0; i < num_pars; ++i) {
      const double g0 = init_parameters(i);
      const double beta = oneOver16PiSqr * 9.0 / (3. * i - 2.);
      const double exact = g0 / Sqrt(1.0 - 2. * beta * Sqr(g0) * (to - from));

      BOOST_CHECK_CLOSE_FRACTION(parameters(i), exact, tol);
   }
}
