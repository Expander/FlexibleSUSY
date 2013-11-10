
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_low_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_MSSM.hpp"
#include <functional>
#include <Eigen/Dense>

#define private public

#include "MSSM_two_scale_model.hpp"
#include "MSSM_two_scale_low_scale_constraint.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"

DoubleVector calculate_gauge_couplings(MSSM<Two_scale> model,
                                       MSSM_low_scale_constraint<Two_scale> constraint,
                                       double scale)
{
   model.set_scale(scale);
   constraint.set_model(&model);
   constraint.apply();

   DoubleVector g(3);
   g(1) = model.get_g1();
   g(2) = model.get_g2();
   g(3) = model.get_g3();

   return g;
}

BOOST_AUTO_TEST_CASE( test_threshold_corrections )
{
   MSSM<Two_scale> m; MssmSoftsusy s;
   MSSM_input_parameters input;
   QedQcd oneset;
   setup_MSSM(m, s, input);

   MSSM_low_scale_constraint<Two_scale> constraint(input, oneset);

   const double Q1 = constraint.get_scale();
   const double Q2 = 2. * Q1;
   const double gut_normalization = 3./5.;
   DoubleVector g_old(3);
   g_old(1) = Electroweak_constants::g1;
   g_old(2) = Electroweak_constants::g2;
   g_old(3) = Electroweak_constants::g3;
   DoubleVector prefactor(3);
   for (int i = 1; i <= 3; i++)
      prefactor(i) = 1. / (oneOver16PiSqr * Power(g_old(i),3));

   const DoubleVector g_Q1(calculate_gauge_couplings(m, constraint, Q1));
   const DoubleVector g_Q2(calculate_gauge_couplings(m, constraint, Q2));

   const DoubleVector beta_numeric((g_Q1 - g_Q2) * prefactor * (1. / log(Q1/Q2)));
   DoubleVector beta_SM(3);
   beta_SM(1) = 41./6. * gut_normalization;
   beta_SM(2) = -19./6.;
   beta_SM(3) = -7. - 2./3;
   DoubleVector beta_MSSM(3);
   beta_MSSM(1) = 11. * gut_normalization;
   beta_MSSM(2) = 1.;
   beta_MSSM(3) = -3.;

   // BOOST_CHECK_CLOSE_FRACTION(beta_numeric(1), beta_MSSM(1) - beta_SM(1), 0.04);
   // BOOST_CHECK_CLOSE_FRACTION(beta_numeric(2), beta_MSSM(2) - beta_SM(2), 0.05);
   BOOST_CHECK_CLOSE_FRACTION(beta_numeric(3), beta_MSSM(3) - beta_SM(3), 0.075);
}

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   MSSM<Two_scale> m; MssmSoftsusy s;
   MSSM_input_parameters input;
   QedQcd oneset;
   setup_MSSM(m, s, input);
   s.setData(oneset);

   MSSM_low_scale_constraint<Two_scale> constraint(input, oneset);
   constraint.set_model(&m);

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

BOOST_AUTO_TEST_CASE( test_low_energy_constraint )
{
   MSSM_input_parameters input;
   QedQcd oneset;
   oneset.setPoleMt(175.);       // non-default
   oneset.setMass(mBottom, 4.3); // non-default
   MSSM<Two_scale> m; MssmSoftsusy s;
   s.setData(oneset);
   setup_MSSM(m, s, input);

   MSSM_low_scale_constraint<Two_scale> constraint(input, oneset);
   constraint.set_model(&m);

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

   const double fs_mt = m.calculate_MFu_DRbar_1loop(oneset.displayPoleMt(), 2);
   const double fs_mb = m.calculate_MFd_DRbar_1loop(oneset.displayMass(mBottom), 2);
   const double fs_me = m.calculate_MFe_DRbar_1loop(oneset.displayMass(mTau), 2);
   const double fs_MZ = m.calculate_MVZ_DRbar_1loop(Electroweak_constants::MZ);
   const double fs_old_vd = m.get_vd();
   const double fs_old_vu = m.get_vu();
   // const double fs_old_vev = Sqrt(Sqr(fs_old_vu) + Sqr(fs_old_vd));
   const double fs_new_vd = (2*fs_MZ)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vu = (2*fs_MZ*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)));
   const double fs_new_vev = Sqrt(Sqr(fs_new_vu) + Sqr(fs_new_vd));

   BOOST_CHECK_CLOSE_FRACTION(fs_mt, ss_mt, 9.5e-5);
   BOOST_CHECK_CLOSE_FRACTION(fs_mb, ss_mb, 3.0e-15);
   BOOST_CHECK_CLOSE_FRACTION(fs_me, ss_me, 4.3e-7);
   BOOST_CHECK_CLOSE_FRACTION(fs_MZ, ss_MZ, 4.5e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_new_vev, ss_new_vev, 4.5e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_old_vu / fs_old_vd, s.displayTanb(), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(fs_new_vu / fs_new_vd, s.displayTanb(), 1.0e-10);

   // apply constraints
   constraint.apply();
   s.sparticleThresholdCorrections(input.TanBeta);

   BOOST_CHECK_CLOSE_FRACTION(m.get_g1(), s.displayGaugeCoupling(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g2(), s.displayGaugeCoupling(2), 0.0075);
   BOOST_CHECK_CLOSE_FRACTION(m.get_g3(), s.displayGaugeCoupling(3), 1.0e-8);

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
