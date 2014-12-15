
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCKM_low_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_CMSSM.hpp"
#include <functional>
#include <Eigen/Dense>

#define private public

#include "CMSSMCKM_two_scale_model.hpp"
#include "CMSSMCKM_two_scale_low_scale_constraint.hpp"
#include "softsusy.h"
#include "flavoursoft.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "ckm.hpp"

BOOST_AUTO_TEST_CASE( test_low_energy_constraint_with_flavour_mixing )
{
   CMSSMCKM_input_parameters input;
   QedQcd oneset;
   oneset.setPoleMt(175.);       // non-default
   oneset.setMass(mBottom, 4.3); // non-default

   // set non-diagonal CKM matrix
   CKM_parameters ckm_parameters;
   ckm_parameters.reset_to_observation();
   oneset.setCKM(ckm_parameters);

   CMSSMCKM<Two_scale> m;
   FlavourMssmSoftsusy s;

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

   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mDown)    , s.displayDataSet().displayMass(mDown)    , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mStrange) , s.displayDataSet().displayMass(mStrange) , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mUp)      , s.displayDataSet().displayMass(mUp)      , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mCharm)   , s.displayDataSet().displayMass(mCharm)   , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mElectron), s.displayDataSet().displayMass(mElectron), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(oneset.displayMass(mMuon)    , s.displayDataSet().displayMass(mMuon)    , 1.0e-10);

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, ss_mt, 9.5e-5);
   BOOST_CHECK_CLOSE_FRACTION(fs_mb, ss_mb, 3.0e-15);
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
