
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_goldstones

#include <boost/test/unit_test.hpp>

#include "wrappers.hpp"

using namespace flexiblesusy;

template <typename DerivedArray, typename DerivedMatrix>
void move_to(unsigned idx, double mass, Eigen::ArrayBase<DerivedArray>& v,
             Eigen::MatrixBase<DerivedMatrix>& z)
{
   const unsigned pos = closest_index(mass, v);

   v.row(idx).swap(v.row(pos));
   z.row(idx).swap(z.row(pos));
}

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

   move_to(0, mvz , higgs, Z);
   move_to(1, mvzp, higgs, Z);

   BOOST_CHECK_EQUAL(higgs(0), mvz_pole);
   BOOST_CHECK_EQUAL(higgs(1), mvzp_pole);
   BOOST_CHECK_EQUAL(higgs(2), mh_pole);

   // ensure that mass matrix is preserved
   const Eigen::Matrix<double, 3, 3> M_new(Z.transpose() * higgs.matrix().asDiagonal() * Z);

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         BOOST_CHECK_EQUAL(M(i,k), M_new(i,k));
}
