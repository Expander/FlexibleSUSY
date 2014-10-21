
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ckm

#include <boost/test/unit_test.hpp>

#include "test.h"
#include "flavoursoft.h"
#include "ckm.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"

using namespace flexiblesusy;

bool is_orthogonal(const Eigen::Matrix<double,3,3>& m,
                   double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<double,3,3> unity(m * m.transpose());

   bool result = is_equal(unity, Eigen::Matrix<double,3,3>::Identity(), max_deviation);

   if (!result)
      BOOST_MESSAGE("matrix is not orthogonal!"
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
      BOOST_MESSAGE("matrix is not unitary!"
                    "\n   Matrix = " << m <<
                    "\n   Matrix * Matrix^+ = " << unity <<
                    "\n   Maximum allowed max_deviation from unity = " << max_deviation);

   return result;
}

BOOST_AUTO_TEST_CASE( test_CKM_unitarity_from_wolfenstein )
{
   CKM_parameters ckm_pars;

   const double lambda = 0.2272;
   const double A      = 0.818;
   const double rhobar = 0.221;
   const double etabar = 0.34;

   ckm_pars.set_from_wolfenstein(lambda, A, rhobar, etabar);

   const Eigen::Matrix<double,3,3> ckm(ckm_pars.get_real_ckm());

   BOOST_CHECK(is_unitary(ckm));
}

BOOST_AUTO_TEST_CASE( test_CKM_unitarity_from_angles )
{
   CKM_parameters ckm_pars;

   ckm_pars.theta_12 = 0.229206;
   ckm_pars.theta_13 = 0.003960;
   ckm_pars.theta_23 = 0.042223;
   ckm_pars.delta    = 0.1;

   const Eigen::Matrix<double,3,3> ckm_real(ckm_pars.get_real_ckm());
   BOOST_CHECK(is_unitary(ckm_real));

   const Eigen::Matrix<std::complex<double>,3,3> ckm_cmpl(ckm_pars.get_complex_ckm());
   BOOST_CHECK(is_unitary(ckm_cmpl));
}

BOOST_AUTO_TEST_CASE( test_CKM_pdg_convention )
{
   // fermion mass matrices
   const Eigen::Matrix<double,3,3>
      mu(Eigen::Matrix<double,3,3>::Random()),
      md(Eigen::Matrix<double,3,3>::Random());

   // mass eigenvalues
   Eigen::Array<double,3,1> su, sd;

   // mixing matrices
   Eigen::Matrix<double,3,3> vu, vd, uu, ud;

   fs_svd(mu, su, uu, vu);
   fs_svd(md, sd, ud, vd);

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));

   Eigen::Matrix<double,3,3> ckm(vu*vd.adjoint());

   ckm.transposeInPlace();
   BOOST_CHECK(is_unitary(ckm, 1.e-10));

   CKM_parameters::to_pdg_convention(ckm, vu, vd, uu, ud);

   BOOST_CHECK(ckm(0,0) > 0.);
   BOOST_CHECK(ckm(1,1) > 0.);
   BOOST_CHECK(ckm(2,2) > 0.);

   BOOST_CHECK(ckm(0,1) > 0.);
   // BOOST_CHECK(ckm(0,2) < 0.); // 13 element is not made positive
   BOOST_CHECK(ckm(1,2) > 0.);

   BOOST_CHECK(is_equal(mu, uu.transpose() * su.matrix().asDiagonal() * vu, 1.e-10));
   BOOST_CHECK(is_equal(md, ud.transpose() * sd.matrix().asDiagonal() * vd, 1.e-10));
}
