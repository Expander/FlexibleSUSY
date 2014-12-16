
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
                    const CMSSMCKM_input_parameters& input)
{
   Eigen::Matrix<double,3,3>
      Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero());
   Eigen::Matrix<double,3,3> mm0;
   Eigen::Matrix<double,3,3> TYuInput;

   const double scale = MZ;
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;

   Yu << 0.0023, 1.e-3, 1.e-3,
         0     , 1.275, 1.e-3,
         0     , 0    , 165.0;

   Yd << 0.0048, 1.e-3, 1.e-3,
         0     , 0.095, 1.e-3,
         0     , 0    , 2.9;

   Ye << 0.000511, 1.e-5, 1.e-5,
         0       , 0.105, 1.e-5,
         0       , 0    , 1.77699;

   Yu *= root2 / vu;
   Yd *= root2 / vd;
   Ye *= root2 / vd;

   mm0 << Sqr(130), 200     , 100,
          200     , Sqr(170), 300,
          100     , 300     , Sqr(200);

   TYuInput << 100, 2  , 3,
               4  , 500, 6,
               7  , 8  , 900;

   m.set_scale(scale);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(125));
   m.set_mHu2(Sqr(150));
   m.set_TYu(TYuInput);
   m.set_TYd(m.get_TYu());
   m.set_TYe(m.get_TYu());

   s.setMu(scale);
   s.setYukawaMatrix(YU, ToDoubleMatrix(Yu));
   s.setYukawaMatrix(YD, ToDoubleMatrix(Yd));
   s.setYukawaMatrix(YE, ToDoubleMatrix(Ye));
   s.setTanb(tanBeta);
   s.setHvev(vev);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);
   s.setGauginoMass(1, M12);
   s.setGauginoMass(2, M12);
   s.setGauginoMass(3, M12);
   s.setSoftMassMatrix(mQl, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mLl, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mUr, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mDr, ToDoubleMatrix(mm0));
   s.setSoftMassMatrix(mEr, ToDoubleMatrix(mm0));
   s.setMh1Squared(Sqr(125));
   s.setMh2Squared(Sqr(150));
   s.setTrilinearMatrix(UA, ToDoubleMatrix(m.get_TYu()));
   s.setTrilinearMatrix(DA, ToDoubleMatrix(m.get_TYu()));
   s.setTrilinearMatrix(EA, ToDoubleMatrix(m.get_TYu()));

   ensure_tree_level_ewsb(m);
   ensure_tree_level_ewsb(s);
}

BOOST_AUTO_TEST_CASE( test_low_energy_constraint_with_flavour_mixing )
{
   CMSSMCKM_input_parameters input;
   input.TanBeta = 10.;
   QedQcd oneset;
   oneset.setPoleMt(175.);       // non-default
   oneset.setMass(mBottom, 4.3); // non-default

   // set non-diagonal CKM matrix
   CKM_parameters ckm_parameters;
   ckm_parameters.reset_to_observation();
   // ckm_parameters.theta_12 = 0.;
   // ckm_parameters.theta_13 = 0.;
   // ckm_parameters.theta_23 = 0.;
   // ckm_parameters.delta = 0.;
   oneset.setCKM(ckm_parameters);

   CMSSMCKM<Two_scale> m(input);
   FlavourMssmSoftsusy s;

   m.set_loops(0);
   s.setLoops(0);

   s.setData(oneset);
   setup_CMSSMCKM(m, s, input);

   s.setTheta12(ckm_parameters.theta_12);
   s.setTheta13(ckm_parameters.theta_13);
   s.setTheta23(ckm_parameters.theta_23);
   s.setDelta(ckm_parameters.delta);
   softsusy::MIXING = 3; // up-type mixing with only one CKM factor

   m.calculate_DRbar_masses();
   s.calcDrBarPars();

   {
      // compare CKM matrices
      const Eigen::Matrix<double,3,3> ckm_fs(ckm_parameters.get_real_ckm());
      const DoubleMatrix ckm_ss(s.displayCkm());
      TEST_CLOSE(ckm_ss, ckm_fs, 1.0e-10);
      BOOST_REQUIRE(gErrors == 0);
   }

   BOOST_CHECK_CLOSE_FRACTION(m.get_scale(), s.displayMu(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MFu(2), s.displayDrBarPars().mt, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MFd(2), s.displayDrBarPars().mb, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_MFe(2), s.displayDrBarPars().mtau, 1.0e-10);
   BOOST_CHECK(is_equal(ToDoubleMatrix(m.get_Yu()), s.displayYukawaMatrix(YU), 1.0e-10));
   BOOST_CHECK(is_equal(ToDoubleMatrix(m.get_Yd()), s.displayYukawaMatrix(YD), 1.0e-10));
   BOOST_CHECK(is_equal(ToDoubleMatrix(m.get_Ye()), s.displayYukawaMatrix(YE), 1.0e-10));
   BOOST_CHECK_CLOSE_FRACTION(m.v(), s.displayHvev(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g1(), s.displayGaugeCoupling(1), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g2(), s.displayGaugeCoupling(2), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g3(), s.displayGaugeCoupling(3), 1.0e-10);

   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mDown)    , s.displayDataSet().displayMass(mDown)    , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mStrange) , s.displayDataSet().displayMass(mStrange) , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mUp)      , s.displayDataSet().displayMass(mUp)      , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mCharm)   , s.displayDataSet().displayMass(mCharm)   , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mElectron), s.displayDataSet().displayMass(mElectron), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mMuon)    , s.displayDataSet().displayMass(mMuon)    , 1.0e-10);

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
   const double fs_MZ = m.calculate_MVZ_DRbar(MZ);
   const double fs_old_vd = m.get_vd();
   const double fs_old_vu = m.get_vu();
   // const double fs_old_vev = Sqrt(Sqr(fs_old_vu) + Sqr(fs_old_vd));
   const double fs_new_vd = (2*fs_MZ)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vu = (2*fs_MZ*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vev = Sqrt(Sqr(fs_new_vu) + Sqr(fs_new_vd));

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, ss_mt, 9.5e-5);
   BOOST_CHECK_CLOSE_FRACTION(fs_mb, ss_mb, 3.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(fs_me, ss_me, 6.0e-7);

   BOOST_CHECK_CLOSE_FRACTION(Re(m.self_energy_VZ(MZ)), pizzt, 4.5e-10);
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

   // The Yukawa couplings differ a lot from Softsusy, because Ben
   // uses the new vev (= the value of the vev after
   // sparticleThresholdCorrections() was called) to calculate the
   // Yukawa couplings.  We use the old vev (= combination of vu, vd
   // from the last run) to calculate the Yukawa couplings.

   // test off-diagonal elements
   BOOST_MESSAGE("testing off-diagonal yukawa elements");
   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         if (i == k)
            continue;
         BOOST_MESSAGE("testing yukawa elements " << i << ", " << k);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(i-1,k-1), s.displayYukawaMatrix(YU)(k,i), 0.011);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(i-1,k-1), s.displayYukawaMatrix(YD)(k,i), 0.00001);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(i-1,k-1), s.displayYukawaMatrix(YE)(k,i), 0.00001);
      }
   }

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
