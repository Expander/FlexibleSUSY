#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_weinberg_angle_new

#include <boost/test/unit_test.hpp>

#include "CMSSM_two_scale_model.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include <Eigen/Core>

#define private public

#include "weinberg_angle.hpp"

using namespace flexiblesusy;
using namespace weinberg_angle;

void ensure_tree_level_ewsb(CMSSM<Two_scale>& m)
{
   // ensure that the EWSB eqs. are satisfied (Drees p.222)
   const double vu = m.get_vu();
   const double vd = m.get_vd();
   const double gY = m.get_g1() * sqrt(0.6);
   const double g2 = m.get_g2();
   const double Mu = m.get_Mu();
   const double BMu = m.get_BMu();
   const double mHd2 = BMu*vu/vd - (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double mHu2 = BMu*vd/vu + (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   m.set_mHd2(mHd2);
   m.set_mHu2(mHu2);
}

// template version of this function already included
// in test_CMSSMNoFV.hpp and test_CMSSM_two_loop_spectrum.cpp
void setup_CMSSM_const(CMSSM<Two_scale>& m, const CMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_input_parameters(input);
   m.set_scale(scale);
   m.set_loops(1);
   m.set_thresholds(3);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();
}

BOOST_AUTO_TEST_CASE( test_1 )
{
    CMSSM<Two_scale> model;
}
