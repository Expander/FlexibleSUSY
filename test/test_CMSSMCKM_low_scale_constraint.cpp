
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCKM_low_scale_constraint

#include <boost/test/unit_test.hpp>
#include <functional>
#include <Eigen/Dense>

#define private public

#include "CMSSMCKM_two_scale_model.hpp"
#include "CMSSMCKM_two_scale_low_scale_constraint.hpp"
#include "flavoursoft.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "ckm.hpp"
#include "conversion.hpp"
#include "test.h"

using namespace flexiblesusy;
using namespace softsusy;

void ensure_tree_level_ewsb(CMSSMCKM<Two_scale>& m)
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

void ensure_tree_level_ewsb(FlavourMssmSoftsusy& softSusy)
{
   const double Mu = softSusy.displaySusyMu();
   // const int signMu = Mu >= 0.0 ? 1 : -1;
   const double vev = softSusy.displayHvev();
   const double tanBeta = softSusy.displayTanb();
   const double beta = atan(tanBeta);
   const double sinBeta = sin(beta);
   const double cosBeta = cos(beta);
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double g1 = softSusy.displayGaugeCoupling(1);
   const double gY = g1 * sqrt(0.6);
   const double g2 = softSusy.displayGaugeCoupling(2);
   const double BMu = softSusy.displayM3Squared();
   const double mHd2 = BMu*vu/vd - (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double mHu2 = BMu*vd/vu + (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double MZrun = 0.5 * vev * sqrt(sqr(gY) + sqr(g2));

   softSusy.setMh1Squared(mHd2);
   softSusy.setMh2Squared(mHu2);

   TEST_CLOSE(MZrun, softSusy.displayMzRun(), 1.0e-10);
   TEST_CLOSE(-2 * BMu, (mHd2 - mHu2) * tan(2*beta) + sqr(MZrun) * sin(2*beta), 1.0e-10);
}

void setup_CMSSMCKM(CMSSMCKM<Two_scale>& m, FlavourMssmSoftsusy& s,
                    CMSSMCKM_input_parameters& input)
{
   input.TanBeta = 10.;

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
   const double M12 = 500.;
   const double m0 = 125.;
   const double a0 = 100.;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu, Yd, Ye;

   // Yu << 0.0023, 1.e-3, 1.e-3,
   //       0     , 1.275, 1.e-3,
   //       0     , 0    , 165.0;

   // Yd << 0.0048, 1.e-3, 1.e-3,
   //       0     , 0.095, 1.e-3,
   //       0     , 0    , 2.9;

   // Ye << 0.000511, 1.e-5, 1.e-5,
   //       0       , 0.105, 1.e-5,
   //       0       , 0    , 1.77699;

   Yu << 0, 0, 0,
         0, 0, 0,
         0, 0, 165.0;

   Yd << 0, 0, 0,
         0, 0, 0,
         0, 0, 2.9;

   Ye << 0, 0, 0,
         0, 0, 0,
         0, 0, 1.77699;

   Yu *= root2 / vu;
   Yd *= root2 / vd;
   Ye *= root2 / vd;

   Eigen::Matrix<double,3,3> mm0;

   mm0 << Sqr(m0), 0        , 0,
          0      , 2*Sqr(m0), 0,
          0      , 0        , 3*Sqr(m0);

   // mm0 << Sqr(130), 200     , 100,
   //        200     , Sqr(170), 300,
   //        100     , 300     , Sqr(200);

   Eigen::Matrix<double,3,3> TYu, TYd, TYe;
   TYu = a0 * Yu;
   TYd = a0 * Yd;
   TYe = a0 * Ye;

   // TYuInput << 100, 2  , 3,
   //             4  , 500, 6,
   //             7  , 8  , 900;

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
   m.set_TYu(TYu);
   m.set_TYd(TYd);
   m.set_TYe(TYe);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   s.setMu(scale);
   s.setLoops(1);
   s.setThresholds(3);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setYukawaMatrix(YU, ToDoubleMatrix(Yu));
   s.setYukawaMatrix(YD, ToDoubleMatrix(Yd));
   s.setYukawaMatrix(YE, ToDoubleMatrix(Ye));
   s.setGauginoMass(1, M12);
   s.setGauginoMass(2, M12);
   s.setGauginoMass(3, M12);
   s.setSoftMassMatrix(mQl, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mUr, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mDr, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mLl, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mEr, ToDoubleMatrix(mm0));
   s.setMh1Squared(sqr(m0));
   s.setMh2Squared(sqr(m0));
   s.setTrilinearMatrix(UA, ToDoubleMatrix(TYu));
   s.setTrilinearMatrix(DA, ToDoubleMatrix(TYd));
   s.setTrilinearMatrix(EA, ToDoubleMatrix(TYe));
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);
   s.setHvev(vev);
   s.setTanb(tanBeta);
   s.setMw(s.displayMwRun());

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_masses();

   ensure_tree_level_ewsb(s);
   s.calcDrBarPars();
}

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   CMSSMCKM<Two_scale> m; FlavourMssmSoftsusy s;
   CMSSMCKM_input_parameters input;
   QedQcd oneset;
   setup_CMSSMCKM(m, s, input);
   s.setData(oneset);

   CMSSMCKM_low_scale_constraint<Two_scale> constraint(&m, input, oneset);

   const double alpha_em = oneset.displayAlpha(ALPHA);
   const double alpha_s  = oneset.displayAlpha(ALPHAS);
   const double scale = m.get_scale();

   const double delta_alpha_em_fs = constraint.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_fs  = constraint.calculate_delta_alpha_s(alpha_s);

   const double delta_alpha_em_ss = 1.0 - alpha_em / s.qedSusythresh(alpha_em, scale);
   const double delta_alpha_s_ss  = 1.0 - alpha_s  / s.qcdSusythresh(alpha_s , scale);

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_fs, delta_alpha_em_ss, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_fs , delta_alpha_s_ss , 1.0e-12);
}

BOOST_AUTO_TEST_CASE( test_low_energy_constraint_with_flavour_mixing )
{
   CMSSMCKM_input_parameters input;
   QedQcd oneset;
   oneset.setPoleMt(175.);       // non-default
   oneset.setMass(mBottom, 4.3); // non-default
   CMSSMCKM<Two_scale> m; FlavourMssmSoftsusy s;

   // set non-diagonal CKM matrix
   CKM_parameters ckm_parameters;
   ckm_parameters.reset_to_observation();
   ckm_parameters.theta_12 = 0.;
   ckm_parameters.theta_13 = 0.;
   ckm_parameters.theta_23 = 0.;
   ckm_parameters.delta = 0.;
   oneset.setCKM(ckm_parameters);

   s.setData(oneset);
   setup_CMSSMCKM(m, s, input);

   s.setTheta12(ckm_parameters.theta_12);
   s.setTheta13(ckm_parameters.theta_13);
   s.setTheta23(ckm_parameters.theta_23);
   s.setDelta(ckm_parameters.delta);
   softsusy::MIXING = 3; // up-type mixing with only one CKM factor

   {
      // compare CKM matrices
      const Eigen::Matrix<double,3,3> ckm_fs(ckm_parameters.get_real_ckm());
      const DoubleMatrix ckm_ss(s.displayCkm());
      TEST_CLOSE(ckm_ss, ckm_fs, 1.0e-10);
      BOOST_REQUIRE(gErrors == 0);
   }

   CMSSMCKM_low_scale_constraint<Two_scale> constraint(&m, input, oneset);

   const double TanBeta = input.TanBeta;
   const double g1 = m.get_g1();
   const double g2 = m.get_g2();

   const double ss_mt = s.calcRunningMt();
   const double ss_mb = s.calcRunningMb();
   const double ss_me = s.calcRunningMtau();
   const double MZ    = s.displayMz();
   const double pizzt = s.piZZT(MZ, s.displayMu());
   const double ss_MZ = Sqrt(Sqr(MZ) + pizzt);
   const double ss_new_vev = s.getVev();

   const double fs_mt = m.calculate_MFu_DRbar(oneset.displayPoleMt(), 2);
   const double fs_mb = m.calculate_MFd_DRbar(oneset.displayMass(mBottom), 2);
   const double fs_me = m.calculate_MFe_DRbar(oneset.displayMass(mTau), 2);
   const double fs_MZ = m.calculate_MVZ_DRbar(Electroweak_constants::MZ);
   const double fs_old_vd = m.get_vd();
   const double fs_old_vu = m.get_vu();
   // const double fs_old_vev = Sqrt(Sqr(fs_old_vu) + Sqr(fs_old_vd));
   const double fs_new_vd = (2*fs_MZ)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vu = (2*fs_MZ*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vev = Sqrt(Sqr(fs_new_vu) + Sqr(fs_new_vd));

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, ss_mt, 9.5e-5);
   BOOST_CHECK_CLOSE_FRACTION(fs_mb, ss_mb, 9.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(fs_me, ss_me, 6.0e-7);
   BOOST_CHECK_CLOSE_FRACTION(fs_MZ, ss_MZ, 4.5e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_new_vev, ss_new_vev, 4.5e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_old_vu / fs_old_vd, s.displayTanb(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_new_vu / fs_new_vd, s.displayTanb(), 1.0e-10);

   // apply constraints
   constraint.apply();
   s.sparticleThresholdCorrections(input.TanBeta);

   BOOST_CHECK_CLOSE_FRACTION(m.get_g1(), s.displayGaugeCoupling(1), 0.0025);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g2(), s.displayGaugeCoupling(2), 0.0070);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g3(), s.displayGaugeCoupling(3), 1.0e-10);

   // test off-diagonal elements
   BOOST_MESSAGE("testing off-diagonal yukawa elements");
   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         if (i == k)
            continue;
         BOOST_MESSAGE("testing yukawa elements " << i << ", " << k);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(i-1,k-1), s.displayYukawaMatrix(YU)(i,k), 0.00001);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(i-1,k-1), s.displayYukawaMatrix(YD)(i,k), 0.00001);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(i-1,k-1), s.displayYukawaMatrix(YE)(i,k), 0.00001);
      }
   }

   // The following Yukawa couplings differ a lot from Softsusy,
   // because Ben uses the new vev (= the value of the vev after
   // sparticleThresholdCorrections() was called) to calculate the
   // Yukawa couplings.  We use the old vev (= combination of vu, vd
   // from the last run) to calculate the Yukawa couplings.

   BOOST_MESSAGE("testing diagonal yukawa elements");
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(0,0), s.displayYukawaMatrix(YU)(1,1), 0.011);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(0,0), s.displayYukawaMatrix(YD)(1,1), 0.011);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(0,0), s.displayYukawaMatrix(YE)(1,1), 0.011);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(1,1), s.displayYukawaMatrix(YU)(2,2), 0.011);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(1,1), s.displayYukawaMatrix(YD)(2,2), 0.011);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(1,1), s.displayYukawaMatrix(YE)(2,2), 0.011);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(2,2), s.displayYukawaMatrix(YU)(3,3), 0.011);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(2,2), s.displayYukawaMatrix(YD)(3,3), 0.011);
   BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(2,2), s.displayYukawaMatrix(YE)(3,3), 0.011);

   BOOST_MESSAGE("testing running VEV");
   const double running_vev = Sqrt(Sqr(m.get_vu()) +  Sqr(m.get_vd()));
   BOOST_CHECK_CLOSE_FRACTION(running_vev, s.displayHvev(), 1.0e-9);
}
