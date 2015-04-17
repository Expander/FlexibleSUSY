
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_goldstones

#include <boost/test/unit_test.hpp>

#include "eigen_utils.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_reordering )
{
   const double mvz  = 90.;
   const double mvzp = 125.;

   const double mvz_pole  = mvz + 1.;
   const double mvzp_pole = mvzp + 1.;
   const double mh_pole   = 250.;

   Eigen::Array<double, 3, 1> higgs;
   higgs(0) = mh_pole;
   higgs(1) = mvz_pole;
   higgs(2) = mvzp_pole;

   Eigen::Matrix<double, 3, 3> Z;
   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         Z(i,k) = 10*i + k;

   // mass matrix
   const Eigen::Matrix<double, 3, 3> M(Z.transpose() * higgs.matrix().asDiagonal() * Z);

   move_goldstone_to(0, mvz , higgs, Z);
   move_goldstone_to(1, mvzp, higgs, Z);

   BOOST_CHECK_EQUAL(higgs(0), mvz_pole);
   BOOST_CHECK_EQUAL(higgs(1), mvzp_pole);
   BOOST_CHECK_EQUAL(higgs(2), mh_pole);

   // ensure that mass matrix is preserved
   const Eigen::Matrix<double, 3, 3> M_new(Z.transpose() * higgs.matrix().asDiagonal() * Z);

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         BOOST_CHECK_EQUAL(M(i,k), M_new(i,k));
}

BOOST_AUTO_TEST_CASE( test_preserve_order )
{
   Eigen::Matrix<double, 3, 3> Z;
   Eigen::Array<double, 3, 1> higgs;
   higgs(0) = 10.;
   higgs(1) = 80.;
   higgs(2) = 91.;

   move_goldstone_to(0, 91., higgs, Z);

   BOOST_CHECK_EQUAL(higgs(0), 91.);
   BOOST_CHECK_EQUAL(higgs(1), 10.);
   BOOST_CHECK_EQUAL(higgs(2), 80.);
}
