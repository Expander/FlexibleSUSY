
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_beta_functions

#include <boost/test/unit_test.hpp>

#include "test.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void setup_NMSSM(NMSSM<Two_scale>& m, NmssmSoftsusy& s, const NMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double lambda = input.LambdaInput;
   const double kappa = 0.01;
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
   const double vS = 1000.;
   const double susyMu = 0;
   const double BMu = 0;
   const double scale = Electroweak_constants::MZ;

   DoubleMatrix Yu_SS(3,3), Yd_SS(3,3), Ye_SS(3,3);
   Yu_SS(3,3) = 165.0   * root2 / (vev * sinBeta);
   Yd_SS(3,3) = 2.9     * root2 / (vev * cosBeta);
   Ye_SS(3,3) = 1.77699 * root2 / (vev * cosBeta);
   DoubleMatrix ID(3, 3), mm0_SS(3, 3);
   for (int i=1; i<=3; i++) ID(i, i) = 1.0;
   mm0_SS = ID * sqr(m0);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
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
   m.set_ms2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_TLambdax(a0);
   m.set_TKappa(a0);
   // m.set_Mu(susyMu);
   // m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vS(vS);

   s.setMu(scale);
   s.setLoops(1);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setYukawaMatrix(YU, Yu_SS);
   s.setYukawaMatrix(YD, Yd_SS);
   s.setYukawaMatrix(YE, Ye_SS);
   s.setLambda(lambda);
   s.setKappa(kappa);
   s.setGauginoMass(1, M12);
   s.setGauginoMass(2, M12);
   s.setGauginoMass(3, M12);
   s.setSoftMassMatrix(mQl, mm0_SS);
   s.setSoftMassMatrix(mUr, mm0_SS);
   s.setSoftMassMatrix(mDr, mm0_SS);
   s.setSoftMassMatrix(mLl, mm0_SS);
   s.setSoftMassMatrix(mEr, mm0_SS);
   s.setMh1Squared(sqr(m0));
   s.setMh2Squared(sqr(m0));
   s.setMsSquared(sqr(m0));
   s.setTrilinearMatrix(UA, a0 * Yu_SS);
   s.setTrilinearMatrix(DA, a0 * Yd_SS);
   s.setTrilinearMatrix(EA, a0 * Ye_SS);
   s.setTrialambda(a0);
   s.setTriakappa(a0);
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);
   s.setHvev(vev);
   s.setTanb(tanBeta);
   s.setSvev(vS);
   s.setMw(s.displayMwRun());
}

void test_parameter_equality(const SoftParsNmssm& a, const NMSSM_soft_parameters& b)
{
   TEST_EQUALITY(a.displayLoops(), b.get_loops());
   TEST_EQUALITY(a.displayMu(), b.get_scale());
   TEST_EQUALITY(a.displayThresholds(), b.get_thresholds());

   TEST_EQUALITY(a.displayGaugeCoupling(1), b.get_g1());
   TEST_EQUALITY(a.displayGaugeCoupling(2), b.get_g2());
   TEST_EQUALITY(a.displayGaugeCoupling(3), b.get_g3());

   TEST_EQUALITY(a.displayYukawaMatrix(YU), b.get_Yu());
   TEST_EQUALITY(a.displayYukawaMatrix(YD), b.get_Yd());
   TEST_EQUALITY(a.displayYukawaMatrix(YE), b.get_Ye());

   TEST_EQUALITY(a.displayGaugino(1), b.get_MassB());
   TEST_EQUALITY(a.displayGaugino(2), b.get_MassWB());
   TEST_EQUALITY(a.displayGaugino(3), b.get_MassG());

   TEST_EQUALITY(a.displayMh1Squared(), b.get_mHd2());
   TEST_EQUALITY(a.displayMh2Squared(), b.get_mHu2());
   TEST_EQUALITY(a.displaySoftMassSquared(mQl), b.get_mq2());
   TEST_EQUALITY(a.displaySoftMassSquared(mUr), b.get_mu2());
   TEST_EQUALITY(a.displaySoftMassSquared(mDr), b.get_md2());
   TEST_EQUALITY(a.displaySoftMassSquared(mLl), b.get_ml2());
   TEST_EQUALITY(a.displaySoftMassSquared(mEr), b.get_me2());

   TEST_EQUALITY(a.displayTrilinear(UA), b.get_TYu());
   TEST_EQUALITY(a.displayTrilinear(DA), b.get_TYd());
   TEST_EQUALITY(a.displayTrilinear(EA), b.get_TYe());

   // TEST_EQUALITY(a.displaySusyMu(), b.get_Mu());
   // TEST_EQUALITY(a.displayM3Squared(), b.get_BMu());

   TEST_EQUALITY(a.displayLambda(), b.get_Lambdax());
   TEST_EQUALITY(a.displayKappa(), b.get_Kappa());
   TEST_EQUALITY(a.displaySvev(), b.get_vS());

   TEST_EQUALITY(a.displayTrialambda(), b.get_TLambdax());
   TEST_EQUALITY(a.displayTriakappa(), b.get_TKappa());
   TEST_EQUALITY(a.displayMsSquared(), b.get_ms2());

   TEST_EQUALITY(a.displayMspSquared(), 0.);
   TEST_EQUALITY(a.displayXiF(), 0.);
   TEST_EQUALITY(a.displayXiS(), 0.);
   TEST_EQUALITY(a.displayMupr(), 0.);

   const double vu = b.get_vu(), vd = b.get_vd();
   double tanBeta;
   if (is_zero(vu))
      tanBeta = 0.;
   else if (is_zero(vd))
      tanBeta = std::numeric_limits<double>::infinity();
   else
      tanBeta = vu/vd;
   const double vev = sqrt(sqr(b.get_vu()) + sqr(b.get_vd()));
   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_beta_functions )
{
   NMSSM_input_parameters input;
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   test_parameter_equality(s, m);
   if (gErrors) {
      BOOST_FAIL("Error: parameters are not equal");
      gErrors = 0;
   }
}
