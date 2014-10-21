
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ckm

#include <boost/test/unit_test.hpp>

#include "flavoursoft.h"
#include "ckm.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

bool is_orthogonal(const Eigen::Matrix<double,3,3>& m,
                   double max_deviation = std::numeric_limits<double>::epsilon())
{
   const Eigen::Matrix<double,3,3> unity(m * m.transpose());

   const double max_dev =
      Abs((unity - Eigen::Matrix<double,3,3>::Identity()).maxCoeff());

   if (max_dev > max_deviation)
      BOOST_MESSAGE("matrix is not orthogonal!"
                    "\n   Matrix = " << m <<
                    "\n   Matrix * Matrix^T = " << unity <<
                    "\n   Maximum max_deviation from unity = " << max_dev <<
                    "\n   Maximum allowed max_deviation from unity = " << max_deviation);

   return max_dev <= max_deviation;
}

bool is_unitary(const Eigen::Matrix<double,3,3>& m)
{
   return is_orthogonal(m);
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
