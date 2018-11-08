
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pmns

#include <boost/test/unit_test.hpp>

#include "test.hpp"
#include "pmns.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"

using namespace flexiblesusy;

bool is_orthogonal(const Eigen::Matrix<double,3,3>& m,
                   double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<double,3,3> unity(m * m.transpose());

   bool result = is_equal(unity, Eigen::Matrix<double,3,3>::Identity(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not orthogonal!"
                    "\n   Matrix = " << m <<
                    "\n   Matrix * Matrix^T = " << unity <<
                    "\n   Maximum allowed max_deviation from unity = " << max_deviation);

   return result;
}

bool is_unitary(const Eigen::Matrix<double,3,3>& m,
                double max_deviation = std::numeric_limits<double>::epsilon())
{
   return is_orthogonal(m, max_deviation);
}

bool is_unitary(const Eigen::Matrix<std::complex<double>,3,3>& m,
                double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<std::complex<double>,3,3> unity(m * m.adjoint());

   bool result = is_equal(unity, Eigen::Matrix<std::complex<double>,3,3>::Identity(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not unitary!"
                    "\n   Matrix = " << m <<
                    "\n   Matrix * Matrix^+ = " << unity <<
                    "\n   Maximum allowed max_deviation from unity = " << max_deviation);

   return result;
}

double random_angle()
{
   return 2. * M_PI * rand() / RAND_MAX;
}

BOOST_AUTO_TEST_CASE( test_PMNS_unitarity_from_angles )
{
   PMNS_parameters pmns_pars;

   pmns_pars.theta_12 = 0.5764;
   pmns_pars.theta_13 = 0.1472;
   pmns_pars.theta_23 = 0.7101;
   pmns_pars.delta = 4.3354;

   pmns_pars.alpha_1 = 0.;
   pmns_pars.alpha_2 = 0.;

   const Eigen::Matrix<double,3,3> pmns_real_no_alpha(pmns_pars.get_real_pmns());
   BOOST_CHECK(is_unitary(pmns_real_no_alpha));

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex_no_alpha(pmns_pars.get_complex_pmns());
   BOOST_CHECK(is_unitary(pmns_complex_no_alpha));

   pmns_pars.alpha_1 = random_angle();
   pmns_pars.alpha_2 = random_angle();

   const Eigen::Matrix<double,3,3> pmns_real(pmns_pars.get_real_pmns());
   BOOST_CHECK(is_unitary(pmns_real));

   const Eigen::Matrix<std::complex<double>,3,3> pmns_complex(pmns_pars.get_complex_pmns());
   BOOST_CHECK(is_unitary(pmns_complex));
}
