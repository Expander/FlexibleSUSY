#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_2L_limits

#include <boost/test/unit_test.hpp>
#include <cmath>

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"

using namespace std;
using namespace flexiblesusy;

constexpr double sqr(double x) { return x*x; }

struct Point {
   double mt2 = sqr(173.34);
   double mg = 1000;
   double mst12 = sqr(4000);
   double mst22 = sqr(5000);
   double st = 0.0000001;
   double ct = sqrt (1. - sqr(st));
   double q2 = sqr(3000);
   double mu = 2000;
   double tb = 5;
   double v2 = sqr(245);
   double g3 = sqrt(4*M_PI*0.1184);
};

// Pietro Slavich implementation
Eigen::Matrix<double, 2, 1> calc_tad_as_at_PS(Point p)
{
   double t1, t2;

   ewsb2loop_(
      &p.mt2, &p.mg, &p.mst12, &p.mst22, &p.st, &p.ct,
      &p.q2, &p.mu, &p.tb, &p.v2, &p.g3, &t1, &t2);

   Eigen::Matrix<double, 2, 1> result;
   result << t1, t2;

   return -result;
}

// FlexibleSUSY wrapper
Eigen::Matrix<double, 2, 1> calc_tad_as_at_FS(Point p)
{
   return mssm_twoloophiggs::tadpole_higgs_2loop_at_as_mssm(
      p.mt2, p.mg, p.mst12, p.mst22, p.st, p.ct, p.q2, p.mu, p.tb, p.v2, p.g3);
}

BOOST_AUTO_TEST_CASE( MSSM_tadpole_at_as_st_0 )
{
   Point p_close, p_exact;
   p_close.st = 0.0000001;
   p_exact.st = 0;

   const auto tad_ps       = calc_tad_as_at_PS(p_close);
   const auto tad_fs_close = calc_tad_as_at_FS(p_close);
   const auto tad_fs_exact = calc_tad_as_at_FS(p_exact);

   BOOST_CHECK_EQUAL(tad_ps(0), tad_fs_close(0));
   BOOST_CHECK_EQUAL(tad_ps(1), tad_fs_close(1));

   BOOST_CHECK_CLOSE(tad_ps(0), tad_fs_exact(0), 1e-3);
   BOOST_CHECK_CLOSE(tad_ps(1), tad_fs_exact(1), 1e-3);

   BOOST_MESSAGE("Pietro Slavich: " << tad_ps.transpose());
   BOOST_MESSAGE("Limit st -> 0 : " << tad_fs_exact.transpose());
}
