
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_rkf_integrator

#include <boost/test/unit_test.hpp>

#include "rkf_integrator.hpp"
#include "basic_rk_integrator.hpp"
#include "error.hpp"

#include <Eigen/Core>

using namespace flexiblesusy;

template <typename Derivs>
void check_odeint_integrator(Eigen::ArrayXd& parameters,
                             double start, double end,
                             Derivs betas)
{
   Eigen::ArrayXd parameters_rkf(parameters);

   const double from = std::log(start);
   const double to = std::log(end);
   const double tol = 1.e-7;

   int err_rk = 0;
   try {
      runge_kutta::Basic_rk_integrator<Eigen::ArrayXd> rk_integrator;
      rk_integrator(from, to, parameters, betas, tol);
   } catch (Error&) {
      err_rk = 1;
   }

   int err_rkf = 0;
   try {
      runge_kutta::RKF_integrator rkf_integrator;
      rkf_integrator(from, to, parameters_rkf, betas, tol);
   } catch (Error&) {
      err_rkf = 1;
   }

   BOOST_CHECK_EQUAL(err_rk, err_rkf);

   if (err_rk == 0 && err_rkf == 0) {
      for (int i = 0; i < parameters.size(); ++i) {
         BOOST_CHECK_CLOSE_FRACTION(parameters(i),
                                    parameters_rkf(i), 1.0e-5);
      }
   }
}

Eigen::ArrayXd beta_one_dim(double /* x */, const Eigen::ArrayXd& parameter)
{
   Eigen::ArrayXd beta(parameter.size());
   beta(0) = parameter(0) * parameter(0) * 0.1;
   return beta;
}

BOOST_AUTO_TEST_CASE( test_one_dim )
{
   Eigen::ArrayXd parameter(1);
   parameter(0) = 0.5;

   check_odeint_integrator(parameter, 100., 1.e10, beta_one_dim);
}

Eigen::ArrayXd beta_ten_dim(double /* x */, const Eigen::ArrayXd& parameter)
{
   Eigen::ArrayXd beta(parameter.size());
   for (int i = 0; i < parameter.size(); ++i) {
      beta(i) = parameter(i) * parameter(0) * 0.1 - 0.5 * i;
   }
   return beta;
}

BOOST_AUTO_TEST_CASE( test_ten_dim )
{
   Eigen::ArrayXd parameter(10);
   for (int i = 0; i < parameter.size(); ++i) {
      parameter(i) = 0.5 + 0.1 * i * i;
   }

   check_odeint_integrator(parameter, 100., 1.e10, beta_ten_dim);
}

Eigen::ArrayXd beta_non_pert(double /* x */, const Eigen::ArrayXd& parameter)
{
   Eigen::ArrayXd beta(parameter.size());
   for (int i = 0; i < parameter.size(); ++i) {
      beta(i) = parameter(i) * parameter(i) + 0.5 * i;
   }
   return beta;
}

BOOST_AUTO_TEST_CASE( test_non_perturbative )
{
   Eigen::ArrayXd parameter(100);
   for (int i = 0; i < parameter.size(); ++i) {
      parameter(i) = 3. + 0.1 * i * i;
   }

   check_odeint_integrator(parameter, 100., 1.e10, beta_non_pert);
}
