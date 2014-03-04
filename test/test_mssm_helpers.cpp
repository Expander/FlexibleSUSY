
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_helpers

#include <boost/test/unit_test.hpp>

#include "MSSM_two_scale_model.hpp"
#include "softsusy.h"
#include "mssm_helpers.hpp"
#include "wrappers.hpp"
#include "test_MSSM.hpp"

using namespace flexiblesusy;

void setup_mssm_models(MSSM<Two_scale>& m, MssmSoftsusy& softSusy)
{
   const int loopLevel = 1;
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * PI * alpha1);
   const double g2 = sqrt(4 * PI * alpha2);
   const double g3 = sqrt(4 * PI * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = sqr(2.0 * susyMu);
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
   mm0 = sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_scale(91);
   m.set_loops(loopLevel);
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
   m.set_mHd2(sqr(m0));
   m.set_mHu2(sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   softSusy.setMu(91);
   softSusy.setLoops(loopLevel);
   softSusy.setGaugeCoupling(1, g1);
   softSusy.setGaugeCoupling(2, g2);
   softSusy.setGaugeCoupling(3, g3);
   softSusy.setYukawaMatrix(YU, Yu_SS);
   softSusy.setYukawaMatrix(YD, Yd_SS);
   softSusy.setYukawaMatrix(YE, Ye_SS);
   softSusy.setGauginoMass(1, M12);
   softSusy.setGauginoMass(2, M12);
   softSusy.setGauginoMass(3, M12);
   softSusy.setSoftMassMatrix(mQl, mm0_SS);
   softSusy.setSoftMassMatrix(mUr, mm0_SS);
   softSusy.setSoftMassMatrix(mDr, mm0_SS);
   softSusy.setSoftMassMatrix(mLl, mm0_SS);
   softSusy.setSoftMassMatrix(mEr, mm0_SS);
   softSusy.setMh1Squared(sqr(m0));
   softSusy.setMh2Squared(sqr(m0));
   softSusy.setTrilinearMatrix(UA, a0 * Yu_SS);
   softSusy.setTrilinearMatrix(DA, a0 * Yd_SS);
   softSusy.setTrilinearMatrix(EA, a0 * Ye_SS);
   softSusy.setSusyMu(susyMu);
   softSusy.setM3Squared(BMu);
   softSusy.setHvev(vev);
   softSusy.setTanb(tanBeta);

   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_parameters();
   softSusy.calcDrBarPars();
}

BOOST_AUTO_TEST_CASE( test_stop )
{
   MSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().mu(1,3);
   mf(2) = s.displayDrBarPars().mu(2,3);
   const double thetat = s.displayDrBarPars().thetat;

   const double vev = s.displayHvev();
   const double tanBeta = s.displayTanb();
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;

   mssm_helpers::Sfermion_mass_data stop_data;
   stop_data.ml2 = s.displaySoftMassSquared(mQl, 3, 3);
   stop_data.mr2 = s.displaySoftMassSquared(mUr, 3, 3);
   stop_data.yf  = s.displayYukawaElement(YU, 3, 3);
   stop_data.vd  = vd;
   stop_data.vu  = vu;
   stop_data.gY  = sqrt(0.6) * s.displayGaugeCoupling(1);
   stop_data.g2  = s.displayGaugeCoupling(2);
   stop_data.Tyf = s.displayTrilinear(UA, 3, 3);
   stop_data.mu  = s.displaySusyMu();
   stop_data.T3  = mssm_helpers::Isospin[mssm_helpers::up];
   stop_data.Yl  = mssm_helpers::Hypercharge_left[mssm_helpers::up];
   stop_data.Yr  = mssm_helpers::Hypercharge_right[mssm_helpers::up];

   Eigen::Array<double,2,1> mstop;
   double theta_stop;

   theta_stop = mssm_helpers::diagonalize_sfermions_2x2(stop_data, mstop);

   BOOST_CHECK_CLOSE_FRACTION(mstop(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(mstop(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(Abs(theta_stop), Abs(thetat), 1.0e-9);
}

BOOST_AUTO_TEST_CASE( test_sbottom )
{
   MSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().md(1,3);
   mf(2) = s.displayDrBarPars().md(2,3);
   double thetab = s.displayDrBarPars().thetab;
   DoubleMatrix Zf(rot2d(thetab));
   Zf.associateOrderAbs(mf);
   thetab = ArcCos(Abs(Zf(1,1)));

   const double vev = s.displayHvev();
   const double tanBeta = s.displayTanb();
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;

   mssm_helpers::Sfermion_mass_data sbottom_data;
   sbottom_data.ml2 = s.displaySoftMassSquared(mQl, 3, 3);
   sbottom_data.mr2 = s.displaySoftMassSquared(mDr, 3, 3);
   sbottom_data.yf  = s.displayYukawaElement(YD, 3, 3);
   sbottom_data.vd  = vd;
   sbottom_data.vu  = vu;
   sbottom_data.gY  = sqrt(0.6) * s.displayGaugeCoupling(1);
   sbottom_data.g2  = s.displayGaugeCoupling(2);
   sbottom_data.Tyf = s.displayTrilinear(DA, 3, 3);
   sbottom_data.mu  = s.displaySusyMu();
   sbottom_data.T3  = mssm_helpers::Isospin[mssm_helpers::down];
   sbottom_data.Yl  = mssm_helpers::Hypercharge_left[mssm_helpers::down];
   sbottom_data.Yr  = mssm_helpers::Hypercharge_right[mssm_helpers::down];

   Eigen::Array<double,2,1> msbottom;
   double theta_sbottom;

   theta_sbottom = mssm_helpers::diagonalize_sfermions_2x2(sbottom_data, msbottom);

   BOOST_CHECK_CLOSE_FRACTION(msbottom(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(msbottom(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_sbottom, thetab, 1.0e-9);
}

BOOST_AUTO_TEST_CASE( test_stau )
{
   MSSM<Two_scale> m;
   MssmSoftsusy s;
   setup_mssm_models(m, s);

   DoubleVector mf(2);
   mf(1) = s.displayDrBarPars().me(1,3);
   mf(2) = s.displayDrBarPars().me(2,3);
   double thetatau = s.displayDrBarPars().thetatau;
   DoubleMatrix Zf(rot2d(thetatau));
   Zf.associateOrderAbs(mf);
   thetatau = ArcCos(Abs(Zf(1,1)));

   const double vev = s.displayHvev();
   const double tanBeta = s.displayTanb();
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;

   mssm_helpers::Sfermion_mass_data stau_data;
   stau_data.ml2 = s.displaySoftMassSquared(mLl, 3, 3);
   stau_data.mr2 = s.displaySoftMassSquared(mEr, 3, 3);
   stau_data.yf  = s.displayYukawaElement(YE, 3, 3);
   stau_data.vd  = vd;
   stau_data.vu  = vu;
   stau_data.gY  = sqrt(0.6) * s.displayGaugeCoupling(1);
   stau_data.g2  = s.displayGaugeCoupling(2);
   stau_data.Tyf = s.displayTrilinear(EA, 3, 3);
   stau_data.mu  = s.displaySusyMu();
   stau_data.T3  = mssm_helpers::Isospin[mssm_helpers::electron];
   stau_data.Yl  = mssm_helpers::Hypercharge_left[mssm_helpers::electron];
   stau_data.Yr  = mssm_helpers::Hypercharge_right[mssm_helpers::electron];

   Eigen::Array<double,2,1> mstau;
   double theta_stau;

   theta_stau = mssm_helpers::diagonalize_sfermions_2x2(stau_data, mstau);

   BOOST_CHECK_CLOSE_FRACTION(mstau(0), mf(1), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(mstau(1), mf(2), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(theta_stau, thetatau, 1.0e-9);
}
