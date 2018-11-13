
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

bool is_symmetric(const Eigen::Matrix<double,3,3>& m,
                  double max_deviation = std::numeric_limits<double>::epsilon())
{
   bool result = is_equal(m, m.transpose(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not symmetric!"
                         "\n   Matrix = " << m <<
                         "\n   Matrix^T = " << m.transpose() <<
                         "\n   Maximum allowed max_deviation from equality = " << max_deviation);

   return result;
}

bool is_symmetric(const Eigen::Matrix<std::complex<double>,3,3>& m,
                  double max_deviation = std::numeric_limits<double>::epsilon())
{
   bool result = is_equal(m, m.transpose(), max_deviation);

   if (!result)
      BOOST_TEST_MESSAGE("matrix is not symmetric!"
                         "\n   Matrix = " << m <<
                         "\n   Matrix^T = " << m.transpose() <<
                         "\n   Maximum allowed max_deviation from equality = " << max_deviation);

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

BOOST_AUTO_TEST_CASE( test_real_PMNS_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<double,3,3> me(Eigen::Matrix<double,3,3>::Random());
   Eigen::Matrix<double,3,3> mv(Eigen::Matrix<double,3,3>::Random());

   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv));

   // mass eigenvalues
   Eigen::Array<double,3,1> se;
   Eigen::Array<double,3,1> sv;

   // mixing matrices
   Eigen::Matrix<double,3,3> ve;
   Eigen::Matrix<double,3,3> ue;
   Eigen::Matrix<double,3,3> vv;

   fs_svd(me, se, ue, ve);
   fs_diagonalize_hermitian(mv, sv, vv);

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve, 1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv, 1.e-10));

   Eigen::Matrix<double,3,3> pmns(ve*vv.adjoint());

   // transpose in order to make (1, 2) negative
   pmns.transposeInPlace();
   BOOST_CHECK(is_unitary(pmns, 1.e-10));

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(pmns(1,2) >= 0.);
   BOOST_CHECK(pmns(2,2) >= 0.);

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve, 1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv, 1.e-10));

   {
      // check that converted matrices are consistent
      const Eigen::Matrix<double,3,3> pmns_check(ve*vv.adjoint());
      for (int i = 0; i < 3; i++) {
         for (int k = 0; k < 3; k++) {
            // check for transposed equality here, because the PMNS
            // matrix from above was transposed to generate a negative
            // (1,2) element
            BOOST_CHECK_CLOSE_FRACTION(Re(pmns(i,k)), Re(pmns_check(k,i)), 1.0e-10);
            BOOST_CHECK_CLOSE_FRACTION(Im(pmns(i,k)), Im(pmns_check(k,i)), 1.0e-10);
         }
      }
   }
}

BOOST_AUTO_TEST_CASE( test_complex_PMNS_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<std::complex<double>,3,3> me(
      Eigen::Matrix<std::complex<double>,3,3>::Random());
   Eigen::Matrix<std::complex<double>,3,3> mv(
      Eigen::Matrix<std::complex<double>,3,3>::Random());

   Symmetrize(mv);

   BOOST_REQUIRE(is_symmetric(mv));

   // mass eigenvalues
   Eigen::Array<double,3,1> se;
   Eigen::Array<double,3,1> sv;

   // mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> ve;
   Eigen::Matrix<std::complex<double>,3,3> ue;
   Eigen::Matrix<std::complex<double>,3,3> vv;

   fs_svd(me, se, ue, ve);
   fs_diagonalize_symmetric(mv, sv, vv);

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve, 1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv, 1.e-10));

   Eigen::Matrix<std::complex<double>,3,3> pmns(ve*vv.adjoint());

   PMNS_parameters::to_pdg_convention(pmns, vv, ve, ue);

   BOOST_CHECK(is_unitary(pmns, 1.e-10));

   BOOST_CHECK(is_equal(me, ue.transpose() * se.matrix().asDiagonal() * ve, 1.e-10));
   BOOST_CHECK(is_equal(mv, vv.transpose() * sv.matrix().asDiagonal() * vv, 1.e-10));

   BOOST_CHECK_LT(std::abs(std::arg(pmns(1,2))), 1e-15);
   BOOST_CHECK_LT(std::abs(std::arg(pmns(2,2))), 1e-15);
}
