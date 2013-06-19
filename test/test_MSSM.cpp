
#include "MSSM_model.hpp"
#include "MSSM_highScaleConstraint.hpp"
#include "test.h"
#include "softsusy.h"
#include "wrappers.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

void test_high_scale_constraint(MSSM m)
{
   double diff;
   int iteration = 0;
   MSSM_input_parameters input;
   input.m0 = 100;
   input.m12 = 200;
   input.TanBeta = m.get_vd() / m.get_vu();
   input.Azero = 300;
   MSSM_high_scale_constraint constraint(input);
   constraint.set_model(&m);
   double mgut = constraint.get_scale(); // initial guess
   double mgut_new = mgut;

   do {
      m.runto(mgut);
      constraint.apply();
      mgut_new = constraint.get_scale();
      diff = std::fabs((mgut - mgut_new)/mgut);
      mgut = mgut_new;
      iteration++;
   } while (diff > 1.0e-7 && iteration < 20);

   TEST_GREATER(1.0e-7, diff);
   TEST_GREATER(20, iteration);

   m.runto(mgut_new);
   TEST_CLOSE(m.get_g1(), m.get_g2(), 1.0e-6);
   TEST_CLOSE(m.get_mHd2(), sqr(input.m0), 1.0e-6);
   TEST_CLOSE(m.get_mHu2(), sqr(input.m0), 1.0e-6);
   TEST_CLOSE(m.get_mq2() , sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_ml2() , sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_md2() , sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_mu2() , sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_me2() , sqr(input.m0) * UNITMATRIX(3), 1.0e-6);
   TEST_CLOSE(m.get_MassB() , input.m12, 1.0e-6);
   TEST_CLOSE(m.get_MassG() , input.m12, 1.0e-6);
   TEST_CLOSE(m.get_MassWB(), input.m12, 1.0e-6);
   TEST_CLOSE(m.get_TYu(), input.Azero * m.get_Yu(), 1.0e-6);
   TEST_CLOSE(m.get_TYd(), input.Azero * m.get_Yd(), 1.0e-6);
   TEST_CLOSE(m.get_TYe(), input.Azero * m.get_Ye(), 1.0e-6);
}

void compare_anomalous_dimensions(const SoftParsMssm& a, const MSSM_soft_parameters& b)
{
  DoubleMatrix gEE(3,3),gLL(3,3),gQQ(3,3),gDD(3,3),gUU(3,3);
  double gH1H1 = 0.0, gH2H2 = 0.0;
  DoubleVector dg(1,3);
  sBrevity brevity;
  a.anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, brevity);

  TEST_EQUALITY(a.displayLoops(), b.displayLoops());
  TEST_EQUALITY(a.displayMu(), b.displayMu());
  TEST_EQUALITY(a.displayThresholds(), b.displayThresholds());

  TEST_EQUALITY(gEE, b.get_SeRSeR());
  TEST_EQUALITY(gLL, b.get_SlSl());
  TEST_EQUALITY(gQQ, b.get_SqSq());
  TEST_EQUALITY(gUU, b.get_SuRSuR());
  TEST_EQUALITY(gDD, b.get_SdRSdR());
  TEST_EQUALITY(gH1H1, b.get_SHdSHd());
  TEST_EQUALITY(gH2H2, b.get_SHuSHu());
}

void test_parameter_equality(const SoftParsMssm& a, const MSSM_soft_parameters& b)
{
   TEST_EQUALITY(a.displayLoops(), b.displayLoops());
   TEST_EQUALITY(a.displayMu(), b.displayMu());
   TEST_EQUALITY(a.displayThresholds(), b.displayThresholds());

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

   TEST_EQUALITY(a.displaySusyMu(), b.get_Mu());
   TEST_EQUALITY(a.displayM3Squared(), b.get_BMu());

   const double tanBeta = b.get_vu() / b.get_vd();
   const double vev = sqrt(sqr(b.get_vu()) + sqr(b.get_vd()));
   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
}

void test_beta_function_equality(const SoftParsMssm& a, const MSSM_soft_parameters& b)
{
   SoftParsMssm beta_a(a.beta2());
   MSSM_soft_parameters beta_b(b.calcBeta());

   TEST_EQUALITY(beta_a.displayLoops(), beta_b.displayLoops());
   TEST_EQUALITY(beta_a.displayMu(), beta_b.displayMu());
   TEST_EQUALITY(beta_a.displayThresholds(), beta_b.displayThresholds());

   TEST_EQUALITY(beta_a.displayGaugeCoupling(1), beta_b.get_g1());
   TEST_EQUALITY(beta_a.displayGaugeCoupling(2), beta_b.get_g2());
   TEST_EQUALITY(beta_a.displayGaugeCoupling(3), beta_b.get_g3());

   TEST_EQUALITY(beta_a.displayYukawaMatrix(YU), beta_b.get_Yu());
   TEST_EQUALITY(beta_a.displayYukawaMatrix(YD), beta_b.get_Yd());
   TEST_EQUALITY(beta_a.displayYukawaMatrix(YE), beta_b.get_Ye());

   TEST_EQUALITY(beta_a.displayGaugino(1), beta_b.get_MassB());
   TEST_EQUALITY(beta_a.displayGaugino(2), beta_b.get_MassWB());
   TEST_EQUALITY(beta_a.displayGaugino(3), beta_b.get_MassG());

   TEST_EQUALITY(beta_a.displayMh1Squared(), beta_b.get_mHd2());
   TEST_EQUALITY(beta_a.displayMh2Squared(), beta_b.get_mHu2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mQl), beta_b.get_mq2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mUr), beta_b.get_mu2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mDr), beta_b.get_md2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mLl), beta_b.get_ml2());
   TEST_EQUALITY(beta_a.displaySoftMassSquared(mEr), beta_b.get_me2());

   TEST_EQUALITY(beta_a.displayTrilinear(UA), beta_b.get_TYu());
   TEST_EQUALITY(beta_a.displayTrilinear(DA), beta_b.get_TYd());
   TEST_EQUALITY(beta_a.displayTrilinear(EA), beta_b.get_TYe());

   TEST_EQUALITY(beta_a.displaySusyMu(), beta_b.get_Mu());
   TEST_EQUALITY(beta_a.displayM3Squared(), beta_b.get_BMu());

   const double vu = b.get_vu();
   const double vd = b.get_vd();
   const double tanBeta = vu / vd;
   const double beta_tanBeta = tanBeta * (beta_b.get_vu()/vu - beta_b.get_vd() / vd);
   const double vev = sqrt(sqr(vu) + sqr(vd));
   const double beta_vev = (vu * beta_b.get_vu() + vd * beta_b.get_vd()) / vev;

   TEST_EQUALITY(a.displayTanb(), tanBeta);
   TEST_EQUALITY(beta_a.displayTanb(), beta_tanBeta);
   TEST_EQUALITY(a.displayHvev(), vev);
   TEST_EQUALITY(beta_a.displayHvev(), beta_vev);
}

void ensure_tree_level_ewsb(MSSM& m)
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

void ensure_tree_level_ewsb(MssmSoftsusy& softSusy)
{
   const double Mu = softSusy.displaySusyMu();
   // const int signMu = Mu >= 0.0 ? 1 : -1;
   const double vev = softSusy.displayHvev();
   const double tanBeta = softSusy.displayTanb();
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double g1 = softSusy.displayGaugeCoupling(1);
   const double gY = g1 * sqrt(0.6);
   const double g2 = softSusy.displayGaugeCoupling(2);
   const double BMu = softSusy.displayM3Squared();
   const double mHd2 = BMu*vu/vd - (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);
   const double mHu2 = BMu*vd/vu + (sqr(gY) + sqr(g2))*(sqr(vd) - sqr(vu))/8. - sqr(Mu);

   // softSusy.rewsbTreeLevel(signMu);

   softSusy.setMh2Squared(mHd2);
   softSusy.setMh1Squared(mHu2);
}

void ensure_one_loop_ewsb(MSSM& m)
{
   ensure_tree_level_ewsb(m);

   const double precision = m.get_ewsb_iteration_precision();
   m.set_ewsb_loop_order(1);
   m.solve_ewsb();
   TEST_CLOSE(m.get_tadpole_vd() - m.tadpole_hh(1).real(), 0.0, precision);
   TEST_CLOSE(m.get_tadpole_vu() - m.tadpole_hh(2).real(), 0.0, precision);
}

void compare_tree_level_masses(MssmSoftsusy s, MSSM m)
{
   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_parameters();
   s.calcDrBarPars();

   const double beta = atan(s.displayTanb());
   const double alpha = s.displayDrBarPars().thetaH;

   // check that tadpole eqs. are fulfilled
   TEST_CLOSE(m.get_tadpole_vd(), 0.0, 1.0e-8);
   TEST_CLOSE(m.get_tadpole_vu(), 0.0, 1.0e-9);

   // neutral CP even Higgs
   DoubleVector hh(m.get_Masshh().sort());
   TEST_EQUALITY(hh(1), s.displayDrBarPars().mh0);
   TEST_EQUALITY(hh(2), s.displayDrBarPars().mH0);
   TEST_CLOSE(m.get_ZH()(1,1), cos(alpha), 1.0e-12);
   TEST_CLOSE(m.get_ZH()(1,2), sin(alpha), 1.0e-12);
   TEST_CLOSE(m.get_ZH()(2,1), -sin(alpha), 1.0e-12);
   TEST_CLOSE(m.get_ZH()(2,2), cos(alpha), 1.0e-12);

   // neutral CP odd Higgs
   DoubleVector Ah(m.get_MassAh().sort());
   TEST_CLOSE(Ah(1), s.displayMzRun(), 1.0e-11);
   TEST_EQUALITY(Ah(2), s.displayDrBarPars().mA0);
   TEST_CLOSE(m.get_ZA()(1,1), sin(beta), 1.0e-12);
   TEST_CLOSE(m.get_ZA()(1,2), cos(beta), 1.0e-12);
   TEST_CLOSE(m.get_ZA()(2,1), -cos(beta), 1.0e-12);
   TEST_CLOSE(m.get_ZA()(2,2), sin(beta), 1.0e-12);

   // charged Higgs
   DoubleVector Hpm(m.get_MassHpm().sort());
   TEST_EQUALITY(Hpm(1), s.displayMwRun()); // for RXi(Wm) == 1
   TEST_EQUALITY(Hpm(2), s.displayDrBarPars().mHpm);
   // This test assumes that we have a twisted rotation using beta a
   // mixing angle.  But in our general approach ZP is an ordinary 2
   // by 2 rotation matrix.
   // TEST_CLOSE(m.get_ZP()(1,1), -cos(beta), 1.0e-12);
   // TEST_CLOSE(m.get_ZP()(1,2), sin(beta), 1.0e-12);
   // TEST_CLOSE(m.get_ZP()(2,1), sin(beta), 1.0e-12);
   // TEST_CLOSE(m.get_ZP()(2,2), cos(beta), 1.0e-12);

   // neutralinos
   DoubleVector mneut(m.get_MassChi());
   TEST_EQUALITY(mneut(1), s.displayDrBarPars().mnBpmz(1));
   TEST_EQUALITY(mneut(2), s.displayDrBarPars().mnBpmz(2));
   TEST_EQUALITY(mneut(3), s.displayDrBarPars().mnBpmz(3));
   TEST_EQUALITY(mneut(4), s.displayDrBarPars().mnBpmz(4));
   const ComplexMatrix ZN(m.get_ZN());
   TEST_CLOSE(ZN, s.displayDrBarPars().nBpmz, 1.0e-12);
   // check neutralino diagonalization convention
   // diag = ZN^* M ZN^\dagger
   const DoubleMatrix m_chi(m.get_mass_matrix_Chi());
   TEST_CLOSE(DoubleMatrix(mneut), ZN.complexConjugate() * m_chi *
              ZN.hermitianConjugate(), 1.0e-10);

   // charginos
   DoubleVector mch(m.get_MassCha());
   TEST_EQUALITY(mch(1), s.displayDrBarPars().mchBpmz(1));
   TEST_EQUALITY(mch(2), s.displayDrBarPars().mchBpmz(2));
   TEST_CLOSE(m.get_UM(), s.displayDrBarPars().uBpmz, 1.0e-12);
   TEST_CLOSE(m.get_UP(), s.displayDrBarPars().vBpmz, 1.0e-12);

   // photon, W and Z mass
   const double vp = m.get_MassVP();
   const double vz = m.get_MassVZ();
   const double vw = m.get_MassVWm();
   TEST_CLOSE(vp, 0.0, 1.0e-6);
   TEST_EQUALITY(vz, s.displayMzRun());
   TEST_EQUALITY(vw, s.displayMwRun());

   // test MSSM tree-level mass relations
   TEST_CLOSE(sqr(hh(1)) + sqr(hh(2)), sqr(Ah(2)) + sqr(vz), 1.0e-9);
   const double vev = sqrt(sqr(m.get_vu()) + sqr(m.get_vd()));
   const double MW = 0.5 * m.get_g2() * vev;
   TEST_CLOSE(sqr(Hpm(2)), sqr(Ah(2)) + MW*MW, 1.0e-9);

   // Compare the sfermion masses

   // Note: In the diagonalization the eigenvalues are ordered
   // automatically.

   // down-type squarks
   const DoubleVector Sd(m.get_MassSd());
   const DoubleMatrix md(s.displayDrBarPars().md);
   const DoubleMatrix ZD(m.get_ZD());
   const double thetab = s.displayDrBarPars().thetab;
   TEST_EQUALITY(Sd(1), md(1,1));
   TEST_EQUALITY(Sd(2), md(1,2));
   TEST_EQUALITY(Sd(3), md(1,3));
   TEST_EQUALITY(Sd(4), md(2,1));
   TEST_EQUALITY(Sd(5), md(2,2));
   TEST_EQUALITY(Sd(6), md(2,3));
   TEST_CLOSE(ZD(1,1), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(1,4), 0.0, 1.0e-12);
   TEST_CLOSE(ZD(4,1), 0.0, 1.0e-12);
   TEST_CLOSE(ZD(4,4), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(2,2), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(2,5), 0.0, 1.0e-12);
   TEST_CLOSE(ZD(5,2), 0.0, 1.0e-12);
   TEST_CLOSE(ZD(5,5), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(3,3), cos(thetab), 1.0e-12);
   TEST_CLOSE(ZD(3,6), sin(thetab), 1.0e-12);
   TEST_CLOSE(ZD(6,3), -sin(thetab), 1.0e-12);
   TEST_CLOSE(ZD(6,6), cos(thetab), 1.0e-12);

   // up-type squarks
   const DoubleVector Su(m.get_MassSu());
   const DoubleMatrix mu(s.displayDrBarPars().mu);
   const DoubleMatrix ZU(m.get_ZU());
   const double thetat = s.displayDrBarPars().thetat;
   TEST_EQUALITY(Su(1), mu(1,1));
   TEST_EQUALITY(Su(2), mu(1,2));
   TEST_EQUALITY(Su(3), mu(1,3));
   TEST_EQUALITY(Su(4), mu(2,1));
   TEST_EQUALITY(Su(5), mu(2,2));
   TEST_EQUALITY(Su(6), mu(2,3));
   TEST_CLOSE(ZU(1,1), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(1,4), 0.0, 1.0e-12);
   TEST_CLOSE(ZU(4,1), 0.0, 1.0e-12);
   TEST_CLOSE(ZU(4,4), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(2,2), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(2,5), 0.0, 1.0e-12);
   TEST_CLOSE(ZU(5,2), 0.0, 1.0e-12);
   TEST_CLOSE(ZU(5,5), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(3,3), cos(thetat), 1.0e-12);
   TEST_CLOSE(ZU(3,6), sin(thetat), 1.0e-12);
   TEST_CLOSE(ZU(6,3), -sin(thetat), 1.0e-12);
   TEST_CLOSE(ZU(6,6), cos(thetat), 1.0e-12);

   // sleptons
   const DoubleVector Se(m.get_MassSe());
   const DoubleMatrix me(s.displayDrBarPars().me);
   const DoubleMatrix ZE(m.get_ZE());
   const double thetatau = s.displayDrBarPars().thetatau;
   TEST_EQUALITY(Se(1), me(1,1));
   TEST_EQUALITY(Se(2), me(1,2));
   TEST_EQUALITY(Se(3), me(1,3));
   TEST_EQUALITY(Se(4), me(2,1));
   TEST_EQUALITY(Se(5), me(2,2));
   TEST_EQUALITY(Se(6), me(2,3));
   TEST_CLOSE(ZE(1,1), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(1,4), 0.0, 1.0e-12);
   TEST_CLOSE(ZE(4,1), 0.0, 1.0e-12);
   TEST_CLOSE(ZE(4,4), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(2,2), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(2,5), 0.0, 1.0e-12);
   TEST_CLOSE(ZE(5,2), 0.0, 1.0e-12);
   TEST_CLOSE(ZE(5,5), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(3,3), cos(thetatau), 1.0e-12);
   TEST_CLOSE(ZE(3,6), sin(thetatau), 1.0e-12);
   TEST_CLOSE(ZE(6,3), -sin(thetatau), 1.0e-12);
   TEST_CLOSE(ZE(6,6), cos(thetatau), 1.0e-12);

   // sneutrinos
   DoubleVector msnu(s.displayDrBarPars().msnu.sort());
   DoubleVector Snu(m.get_MassSv());

   TEST_EQUALITY(Snu(1), msnu(1));
   TEST_EQUALITY(Snu(2), msnu(2));
   TEST_EQUALITY(Snu(3), msnu(3));

   // gluons
   TEST_EQUALITY(m.get_MassVG(), 0.0);

   // gluinos
   TEST_EQUALITY(m.get_MassGlu(), s.displayDrBarPars().mGluino);

   // neutrinos
   TEST_EQUALITY(m.get_MassFv()(1), 0.0);
   TEST_EQUALITY(m.get_MassFv()(2), 0.0);
   TEST_EQUALITY(m.get_MassFv()(3), 0.0);

   // leptons
   TEST_EQUALITY(m.get_MassFe()(1), 0.0);
   TEST_EQUALITY(m.get_MassFe()(2), 0.0);
   TEST_EQUALITY(m.get_MassFe()(3), s.displayDrBarPars().mtau);
   DoubleMatrix unity(3,3);
   unity(1,1) = -1.0; // why is this chosen to be negative?
   unity(2,2) = -1.0; // why is this chosen to be negative?
   unity(3,3) = 1.0;
   TEST_EQUALITY(m.get_ZEL(), unity);
   TEST_EQUALITY(m.get_ZER(), unity);

   // ups
   TEST_EQUALITY(m.get_MassFu()(1), 0.0);
   TEST_EQUALITY(m.get_MassFu()(2), 0.0);
   TEST_EQUALITY(m.get_MassFu()(3), s.displayDrBarPars().mt);
   TEST_EQUALITY(m.get_ZUL(), unity);
   TEST_EQUALITY(m.get_ZUR(), unity);

   // downs
   TEST_EQUALITY(m.get_MassFd()(1), 0.0);
   TEST_EQUALITY(m.get_MassFd()(2), 0.0);
   TEST_EQUALITY(m.get_MassFd()(3), s.displayDrBarPars().mb);
   TEST_EQUALITY(m.get_ZDL(), unity);
   TEST_EQUALITY(m.get_ZDR(), unity);
}

void compare_gluino_self_energy(MssmSoftsusy s, MSSM m)
{
   // tree-level
   s.gluino(0);
   const double m3 = s.displayGaugino(3);
   const double g3 = s.displayGaugeCoupling(3);
   const double softsusy_gluino_tree = s.displayPhys().mGluino;

   TEST_CLOSE(softsusy_gluino_tree, m3, 1.0e-12);
   TEST_CLOSE(softsusy_gluino_tree, m.get_MassGlu(), 1.0e-12);
   TEST_EQUALITY(g3, m.get_g3());

   // one-loop
   s.gluino(1);
   const double softsusy_gluino_se = s.displayPhys().mGluino - softsusy_gluino_tree;
   const double p = std::fabs(m3);
   const double glu_scalar = m.self_energy_Glu_1(p).real();
   const double glu_left   = m.self_energy_Glu_PL(p).real();
   const double glu_right  = m.self_energy_Glu_PR(p).real();
   const double glu_se     = - (glu_scalar + m3 * (glu_left + glu_right));

   TEST_CLOSE(softsusy_gluino_se, glu_se, 1.0e-4);
}

void compare_neutralino_self_energy(MssmSoftsusy s, MSSM m)
{
   const double p = s.displayDrBarPars().mneut(1);

   const double tanb = s.displayTanb();
   const double cosb = cos(atan(tanb));
   const double m1 = s.displayGaugino(1);
   const double m2 = s.displayGaugino(2);
   const double smu = s.displaySusyMu();
   DoubleMatrix softsusy_tree(4, 4); // tree-level mass matrix

   // fill tree-level mass matrix
   softsusy_tree(1, 1) = m1;
   softsusy_tree(2, 2) = m2;
   softsusy_tree(1, 3) = - s.displayMzRun() * cosb * s.calcSinthdrbar();
   softsusy_tree(1, 4) = - softsusy_tree(1, 3) * tanb;
   softsusy_tree(2, 3) = s.displayMwRun() * cosb;
   softsusy_tree(2, 4) = - softsusy_tree(2, 3) * tanb;
   softsusy_tree(3, 4) = - smu;
   softsusy_tree.symmetrise();

   TEST_EQUALITY(softsusy_tree, softsusy_tree.transpose());

   ComplexMatrix sarah_sigma_L(4,4), sarah_sigma_R(4,4), sarah_sigma_S(4,4);
   for (unsigned i = 1; i <= 4; ++i) {
      for (unsigned k = 1; k <= 4; ++k) {
         sarah_sigma_L(i,k) = m.self_energy_Chi_PL(p,i,k);
         sarah_sigma_R(i,k) = m.self_energy_Chi_PR(p,i,k);
         sarah_sigma_S(i,k) = m.self_energy_Chi_1(p,i,k);
      }
   }

   DoubleVector Chi(m.get_MassChi());
   ComplexMatrix M_tree(m.get_mass_matrix_Chi());

   // check that tree-level mass matrix is real
   TEST_EQUALITY(M_tree.imag(), DoubleMatrix(4,4));

   // check that SoftSusy and SARAH give the same tree-level mass
   // matrix
   TEST_EQUALITY(M_tree.real(), softsusy_tree);

   // calculate SoftSusy self-energy
   DoubleMatrix softsusy_sigma(softsusy_tree);
   s.addNeutralinoLoop(p, softsusy_sigma);
   softsusy_sigma = softsusy_sigma - softsusy_tree;

   ComplexMatrix sarah_delta_M(- sarah_sigma_R * M_tree
                               - M_tree * sarah_sigma_L
                               - sarah_sigma_S);
   ComplexMatrix sarah_sigma(0.5 * (sarah_delta_M + sarah_delta_M.transpose()));

   TEST_EQUALITY(sarah_sigma.imag(), DoubleMatrix(4,4));
   TEST_EQUALITY(sarah_sigma.real(), sarah_sigma.real().transpose());
   TEST_CLOSE(softsusy_sigma, sarah_sigma.real(), 1.0e-10);
}

void compare_chargino_self_energy(MssmSoftsusy s, MSSM m)
{
   const double p = std::fabs(s.displayDrBarPars().mch(1));

   const double tanb = s.displayTanb();
   const double cosb = cos(atan(tanb));
   const double m2 = s.displayGaugino(2);
   const double smu = s.displaySusyMu();
   const double mwOneLarg = sqr(s.displayMwRun());
   DoubleMatrix softsusy_tree(2, 2); // tree-level mass matrix

   // fill tree-level mass matrix
   softsusy_tree(1, 1) = m2;
   softsusy_tree(2, 1) = root2 * sqrt(fabs(mwOneLarg)) * cosb;
   softsusy_tree(1, 2) = softsusy_tree(2, 1) * tanb;
   softsusy_tree(2, 2) = smu;

   ComplexMatrix sarah_sigma_L(2,2), sarah_sigma_R(2,2), sarah_sigma_S(2,2);
   for (unsigned i = 1; i <= 2; ++i) {
      for (unsigned k = 1; k <= 2; ++k) {
         sarah_sigma_L(i,k) = m.self_energy_Cha_PL(p,i,k);
         sarah_sigma_R(i,k) = m.self_energy_Cha_PR(p,i,k);
         sarah_sigma_S(i,k) = m.self_energy_Cha_1(p,i,k);
      }
   }

   DoubleVector Cha(m.get_MassCha());
   ComplexMatrix M_tree(m.get_mass_matrix_Cha());

   // check that tree-level mass matrix is real
   TEST_EQUALITY(M_tree.imag(), DoubleMatrix(2,2));

   // check that SoftSusy and SARAH give the same tree-level mass
   // matrix
   TEST_EQUALITY(M_tree.real(), softsusy_tree);

   // calculate SoftSusy self-energy
   DoubleMatrix softsusy_sigma(softsusy_tree);
   s.addCharginoLoop(p, softsusy_sigma);
   softsusy_sigma = softsusy_sigma - softsusy_tree;

   ComplexMatrix sarah_sigma(- sarah_sigma_R * M_tree
                             - M_tree * sarah_sigma_L
                             - sarah_sigma_S);

   TEST_EQUALITY(sarah_sigma.imag(), DoubleMatrix(2,2));
   TEST_CLOSE(softsusy_sigma, sarah_sigma.real(), 1.0e-10);
}

void compare_sneutrino_self_energy(MssmSoftsusy s, MSSM m)
{
   // tree-level
   s.doSnu(0.0, 0);
   const DoubleVector Snu_softsusy_tree(s.displayPhys().msnu);
   DoubleVector Snu_sarah_tree(m.get_MassSv());
   TEST_CLOSE(Snu_softsusy_tree(1), Snu_sarah_tree(1), 1.0e-10);
   TEST_CLOSE(Snu_softsusy_tree(2), Snu_sarah_tree(2), 1.0e-10);
   TEST_CLOSE(Snu_softsusy_tree(3), Snu_sarah_tree(3), 1.0e-10);

   // one-loop
   s.doSnu(0.0, 1);
   const DoubleVector Snu_softsusy_1loop(s.displayPhys().msnu);
   ComplexMatrix Snu_sarah_se(3,3);
   for (unsigned i1 = 1; i1 <= 3; ++i1) {
      const double p = Snu_sarah_tree(i1);
      for (unsigned i2 = 1; i2 <= 3; ++i2) {
         Snu_sarah_se(i1, i2) = m.self_energy_Sv(p, i1, i2);
      }
   }

   // check self-energy matrix is hermitian
   TEST_CLOSE(Snu_sarah_se, Snu_sarah_se.hermitianConjugate(), 1.0e-10);

   // check that off-diagonal self-energies are zero
   TEST_CLOSE(Snu_sarah_se(1,2), Complex(0,0), 1.0e-10);
   TEST_CLOSE(Snu_sarah_se(1,3), Complex(0,0), 1.0e-10);
   TEST_CLOSE(Snu_sarah_se(2,3), Complex(0,0), 1.0e-10);

   DoubleVector Snu_sarah_1loop(3);
   for (unsigned i = 1; i <= 3; ++i) {
      Snu_sarah_1loop(i) = zeroSqrt(sqr(Snu_sarah_tree(i))
                                     - Snu_sarah_se(i,i).real());
   }
   TEST_CLOSE(Snu_softsusy_1loop(1), Snu_sarah_1loop(1), 1.0e-10);
   TEST_CLOSE(Snu_softsusy_1loop(2), Snu_sarah_1loop(2), 1.0e-10);
   TEST_CLOSE(Snu_softsusy_1loop(3), Snu_sarah_1loop(3), 1.0e-10);
}

void compare_selectron_self_energy(MssmSoftsusy s, MSSM m)
{
   const double mtau = s.displayDrBarPars().mtau;
   const double sinthDRbar = s.calcSinthdrbar();

   // tree-level
   s.doChargedSleptons(mtau, 0.0, sinthDRbar, 0);
   const DoubleMatrix Se_softsusy_tree(s.displayPhys().me);
   const DoubleVector Se_sarah_tree(m.get_MassSe());
   TEST_EQUALITY(Se_sarah_tree(1), Se_softsusy_tree(1,1));
   TEST_EQUALITY(Se_sarah_tree(2), Se_softsusy_tree(1,2));
   TEST_EQUALITY(Se_sarah_tree(3), Se_softsusy_tree(1,3));
   TEST_EQUALITY(Se_sarah_tree(4), Se_softsusy_tree(2,1));
   TEST_EQUALITY(Se_sarah_tree(5), Se_softsusy_tree(2,2));
   TEST_EQUALITY(Se_sarah_tree(6), Se_softsusy_tree(2,3));

   // one-loop
   s.doChargedSleptons(mtau, 0.0, sinthDRbar, 1);
   const DoubleMatrix Se_softsusy_1loop(s.displayPhys().me);
   ComplexMatrix Se_sarah_se(6,6);
   for (unsigned i1 = 1; i1 <= 6; ++i1) {
      const double p = Se_sarah_tree(i1);
      for (unsigned i2 = 1; i2 <= 6; ++i2) {
         Se_sarah_se(i1, i2) = m.self_energy_Se(p, i1, i2);
      }
   }

   DoubleMatrix Se_softsusy_se(6,6);
   for (unsigned i = 1; i <= 2; ++i) {
      DoubleMatrix mat(2,2);
      s.addSlepCorrection(mat, i);
      Se_softsusy_se(i,i)     = -mat(1,1);
      Se_softsusy_se(i+3,i+3) = -mat(2,2);
      Se_softsusy_se(i+3,i)   = -mat(2,1);
      Se_softsusy_se(i,i+3)   = -mat(1,2);
   }
   {
      DoubleMatrix mSlepSquared_1(2,2), mSlepSquared_2(2,2);
      s.addStauCorrection(Se_softsusy_tree(1,3), mSlepSquared_1, mtau);
      s.addStauCorrection(Se_softsusy_tree(2,3), mSlepSquared_2, mtau);
      Se_softsusy_se(3,3)     = -mSlepSquared_1(1,1);
      Se_softsusy_se(3,6)     = -mSlepSquared_1(1,2);
      Se_softsusy_se(6,3)     = -mSlepSquared_2(2,1);
      Se_softsusy_se(6,6)     = -mSlepSquared_2(2,2);
   }

   // compare families 1 and 2
   TEST_CLOSE(Se_softsusy_se(1,1), Se_sarah_se(1,1), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(1,2), Se_sarah_se(1,2), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(2,1), Se_sarah_se(2,1), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(2,2), Se_sarah_se(2,2), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(4,4), Se_sarah_se(4,4), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(4,5), Se_sarah_se(4,5), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(5,4), Se_sarah_se(5,4), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(5,5), Se_sarah_se(5,5), 1.0e-10);
   // compare 3rd family
   TEST_CLOSE(Se_softsusy_se(3,3), Se_sarah_se(3,3), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(3,6), Se_sarah_se(3,6), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(6,3), Se_sarah_se(6,3), 1.0e-10);
   TEST_CLOSE(Se_softsusy_se(6,6), Se_sarah_se(6,6), 1.0e-10);

   // one-loop diagonalization
   DoubleMatrix MSe_sarah_1loop(m.get_mass_matrix_Se()
                                - Se_sarah_se.real());
   DoubleMatrix ZSe_sarah_1loop(6,6);
   DoubleVector Se_sarah_1loop(6);
   DiagonalizeUnsorted(MSe_sarah_1loop, ZSe_sarah_1loop, Se_sarah_1loop);
   Se_sarah_1loop = Se_sarah_1loop.apply(zeroSqrt);

   // The differences of 0.02 and 0.03 GeV come from the different way
   // we do the diagnoalization.
   TEST_CLOSE(Se_softsusy_1loop(1,1), Se_sarah_1loop(1), 1.0e-10);
   TEST_CLOSE(Se_softsusy_1loop(1,2), Se_sarah_1loop(2), 1.0e-10);
   TEST_CLOSE(Se_softsusy_1loop(1,3), Se_sarah_1loop(3), 0.02);
   TEST_CLOSE(Se_softsusy_1loop(2,1), Se_sarah_1loop(4), 1.0e-10);
   TEST_CLOSE(Se_softsusy_1loop(2,2), Se_sarah_1loop(5), 1.0e-10);
   TEST_CLOSE(Se_softsusy_1loop(2,3), Se_sarah_1loop(6), 0.03);
}

void compare_sup_self_energy(MssmSoftsusy s, MSSM m)
{
   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();

   // tree-level
   s.doUpSquarks(mt, 0.0, sinthDRbar, 0);
   const DoubleMatrix Su_softsusy_tree(s.displayPhys().mu);
   const DoubleVector Su_sarah_tree(m.get_MassSu());
   TEST_EQUALITY(Su_sarah_tree(1), Su_softsusy_tree(1,1));
   TEST_EQUALITY(Su_sarah_tree(2), Su_softsusy_tree(1,2));
   TEST_EQUALITY(Su_sarah_tree(3), Su_softsusy_tree(1,3));
   TEST_EQUALITY(Su_sarah_tree(4), Su_softsusy_tree(2,1));
   TEST_EQUALITY(Su_sarah_tree(5), Su_softsusy_tree(2,2));
   TEST_EQUALITY(Su_sarah_tree(6), Su_softsusy_tree(2,3));

   // one-loop
   s.doUpSquarks(mt, 0.0, sinthDRbar, 1);
   const DoubleMatrix Su_softsusy_1loop(s.displayPhys().mu);
   ComplexMatrix Su_sarah_se(6,6);
   for (unsigned i1 = 1; i1 <= 6; ++i1) {
      for (unsigned i2 = 1; i2 <= 6; ++i2) {
         const double p = sqrt(Su_sarah_tree(i1) * Su_sarah_tree(i2));
         Su_sarah_se(i1, i2) = m.self_energy_Su(p, i1, i2);
      }
   }

   DoubleMatrix Su_softsusy_se(6,6);
   for (unsigned i = 1; i <= 2; ++i) {
      DoubleMatrix mat(2,2);
      s.addSupCorrection(mat, i);
      Su_softsusy_se(i,i)     = -mat(1,1);
      Su_softsusy_se(i+3,i+3) = -mat(2,2);
      Su_softsusy_se(i+3,i)   = -mat(2,1);
      Su_softsusy_se(i,i+3)   = -mat(1,2);
   }

   // compare families 1 and 2
   TEST_CLOSE(Su_softsusy_se(1,1), Su_sarah_se(1,1), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(1,2), Su_sarah_se(1,2), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(2,1), Su_sarah_se(2,1), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(2,2), Su_sarah_se(2,2), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(4,4), Su_sarah_se(4,4), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(4,5), Su_sarah_se(4,5), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(5,4), Su_sarah_se(5,4), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(5,5), Su_sarah_se(5,5), 1.0e-10);

   // compare 3rd family
   {
      // for simplicity we take only one fixed momentum
      const double p = Su_softsusy_tree(1,3);

      DoubleMatrix mStopSquared(2,2);
      s.addStopCorrection(p, mStopSquared, mt);

      Su_softsusy_se(3,3) = -mStopSquared(1,1);
      Su_softsusy_se(3,6) = -mStopSquared(1,2);
      Su_softsusy_se(6,3) = -mStopSquared(2,1);
      Su_softsusy_se(6,6) = -mStopSquared(2,2);

      Su_sarah_se(3,3) = m.self_energy_Su(p,3,3);
      Su_sarah_se(3,6) = m.self_energy_Su(p,3,6);
      Su_sarah_se(6,3) = m.self_energy_Su(p,6,3);
      Su_sarah_se(6,6) = m.self_energy_Su(p,6,6);
   }

   TEST_CLOSE(Su_softsusy_se(3,3), Su_sarah_se(3,3), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(3,6), Su_sarah_se(3,6), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(6,3), Su_sarah_se(6,3), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(6,6), Su_sarah_se(6,6), 1.0e-9);
}

void compare_sdown_self_energy(MssmSoftsusy s, MSSM m)
{
   const double mt = s.displayDrBarPars().mt;
   const double mb = s.displayDrBarPars().mb;
   const double sinthDRbar = s.calcSinthdrbar();

   // tree-level
   s.doDownSquarks(mb, 0.0, sinthDRbar, 0, mt);
   const DoubleMatrix Sd_softsusy_tree(s.displayPhys().md);
   const DoubleVector Sd_sarah_tree(m.get_MassSd());
   TEST_EQUALITY(Sd_sarah_tree(1), Sd_softsusy_tree(1,1));
   TEST_EQUALITY(Sd_sarah_tree(2), Sd_softsusy_tree(1,2));
   TEST_EQUALITY(Sd_sarah_tree(3), Sd_softsusy_tree(1,3));
   TEST_EQUALITY(Sd_sarah_tree(4), Sd_softsusy_tree(2,1));
   TEST_EQUALITY(Sd_sarah_tree(5), Sd_softsusy_tree(2,2));
   TEST_EQUALITY(Sd_sarah_tree(6), Sd_softsusy_tree(2,3));

   // one-loop
   s.doDownSquarks(mb, 0.0, sinthDRbar, 1, mt);
   const DoubleMatrix Sd_softsusy_1loop(s.displayPhys().md);
   ComplexMatrix Sd_sarah_se(6,6);
   for (unsigned i1 = 1; i1 <= 6; ++i1) {
      for (unsigned i2 = 1; i2 <= 6; ++i2) {
         const double p = sqrt(Sd_sarah_tree(i1) * Sd_sarah_tree(i2));
         Sd_sarah_se(i1, i2) = m.self_energy_Sd(p, i1, i2);
      }
   }

   DoubleMatrix Sd_softsusy_se(6,6);
   for (unsigned i = 1; i <= 2; ++i) {
      DoubleMatrix mat(2,2);
      s.addSdownCorrection(mat, i);
      Sd_softsusy_se(i,i)     = -mat(1,1);
      Sd_softsusy_se(i+3,i+3) = -mat(2,2);
      Sd_softsusy_se(i+3,i)   = -mat(2,1);
      Sd_softsusy_se(i,i+3)   = -mat(1,2);
   }

   // compare families 1 and 2
   TEST_CLOSE(Sd_softsusy_se(1,1), Sd_sarah_se(1,1), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(1,2), Sd_sarah_se(1,2), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(2,1), Sd_sarah_se(2,1), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(2,2), Sd_sarah_se(2,2), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(4,4), Sd_sarah_se(4,4), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(4,5), Sd_sarah_se(4,5), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(5,4), Sd_sarah_se(5,4), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(5,5), Sd_sarah_se(5,5), 1.0e-10);

   // compare 3rd family
   {
      // for simplicity we take only one fixed momentum
      const double p = Sd_softsusy_tree(1,3);

      DoubleMatrix mSbotSquared(2,2);
      s.addSbotCorrection(p, mSbotSquared, mt);

      Sd_softsusy_se(3,3) = -mSbotSquared(1,1);
      Sd_softsusy_se(3,6) = -mSbotSquared(1,2);
      Sd_softsusy_se(6,3) = -mSbotSquared(2,1);
      Sd_softsusy_se(6,6) = -mSbotSquared(2,2);

      Sd_sarah_se(3,3) = m.self_energy_Sd(p,3,3);
      Sd_sarah_se(3,6) = m.self_energy_Sd(p,3,6);
      Sd_sarah_se(6,3) = m.self_energy_Sd(p,6,3);
      Sd_sarah_se(6,6) = m.self_energy_Sd(p,6,6);
   }

   TEST_CLOSE(Sd_softsusy_se(3,3), Sd_sarah_se(3,3), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(3,6), Sd_sarah_se(3,6), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(6,3), Sd_sarah_se(6,3), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(6,6), Sd_sarah_se(6,6), 1.2e-10);
}

void compare_CP_even_higgs_self_energy(MssmSoftsusy s, MSSM m)
{
   const double mh0 = s.displayDrBarPars().mh0;
   const double mH0 = s.displayDrBarPars().mH0;
   const double scale = s.displayMu();
   const DoubleVector hh(m.get_Masshh());

   TEST_EQUALITY(s.displayMu(), m.displayMu());

   if (hh(1) <= hh(2)) {
      TEST_CLOSE(mh0, hh(1), 1.0e-10);
      TEST_CLOSE(mH0, hh(2), 1.0e-10);
   } else {
      TEST_CLOSE(mh0, hh(2), 1.0e-10);
      TEST_CLOSE(mH0, hh(1), 1.0e-10);
   }

   DoubleMatrix softsusy_sigma_light(2,2);
   softsusy_sigma_light(1,1) = s.pis1s1(mh0, scale);
   softsusy_sigma_light(1,2) = s.pis1s2(mh0, scale);
   softsusy_sigma_light(2,1) = softsusy_sigma_light(1,2);
   softsusy_sigma_light(2,2) = s.pis2s2(mh0, scale);

   DoubleMatrix softsusy_sigma_heavy(2,2);
   softsusy_sigma_heavy(1,1) = s.pis1s1(mH0, scale);
   softsusy_sigma_heavy(1,2) = s.pis1s2(mH0, scale);
   softsusy_sigma_heavy(2,1) = softsusy_sigma_heavy(1,2);
   softsusy_sigma_heavy(2,2) = s.pis2s2(mH0, scale);

   ComplexMatrix sarah_sigma_light(2,2), sarah_sigma_heavy(2,2);
   for (unsigned i1 = 1; i1 <= 2; ++i1) {
      for (unsigned i2 = 1; i2 <= 2; ++i2) {
         sarah_sigma_light(i1,i2) = m.self_energy_hh(mh0,i1,i2);
         sarah_sigma_heavy(i1,i2) = m.self_energy_hh(mH0,i1,i2);
      }
   }

   TEST_CLOSE(sarah_sigma_light.imag(), DoubleMatrix(2,2), 1.0e-10);
   TEST_CLOSE(sarah_sigma_light.real(), sarah_sigma_light.real().transpose(), 1.0e-10);
   TEST_CLOSE(sarah_sigma_light.real(), softsusy_sigma_light, 1.0e-10);

   TEST_CLOSE(sarah_sigma_heavy.imag(), DoubleMatrix(2,2), 1.0e-10);
   TEST_CLOSE(sarah_sigma_heavy.real(), sarah_sigma_heavy.real().transpose(), 1.0e-10);
   TEST_CLOSE(sarah_sigma_heavy.real(), softsusy_sigma_heavy, 1.0e-10);
}

void compare_CP_odd_higgs_self_energy(MssmSoftsusy s, MSSM m)
{
   const double mA0 = s.displayDrBarPars().mA0;
   const double mZrun = s.displayMzRun();
   const double scale = s.displayMu();
   const DoubleVector Ah(m.get_MassAh());

   TEST_EQUALITY(s.displayMu(), m.displayMu());

   bool twisted;    // true if Ah(1) == A0, false otherwise
   const double p = mA0;

   if (Ah(1) <= Ah(2)) {
      TEST_CLOSE(mZrun, Ah(1), 1.0e-10);
      TEST_CLOSE(mA0  , Ah(2), 1.0e-10);
      twisted = false;
   } else {
      TEST_CLOSE(mZrun, Ah(2), 1.0e-10);
      TEST_CLOSE(mA0  , Ah(1), 1.0e-10);
      twisted = true;
   }

   const double softsusy_sigma_AA = s.piAA(mA0, scale);

   ComplexMatrix sarah_sigma_AA(2,2);
   for (unsigned i1 = 1; i1 <= 2; ++i1) {
      for (unsigned i2 = 1; i2 <= 2; ++i2) {
         sarah_sigma_AA(i1,i2) = m.self_energy_Ah(p,i1,i2);
      }
   }

   // do tree-level rotation to compare with SoftSusy
   const DoubleMatrix tree_tevel_rotation(m.get_ZA());
   sarah_sigma_AA = tree_tevel_rotation * sarah_sigma_AA
      * tree_tevel_rotation.transpose();

   TEST_CLOSE(sarah_sigma_AA.imag(), DoubleMatrix(2,2), 1.0e-10);
   TEST_CLOSE(sarah_sigma_AA.real(), sarah_sigma_AA.real().transpose(), 1.0e-10);

   if (twisted) {
      TEST_CLOSE(sarah_sigma_AA.real()(1,1), softsusy_sigma_AA, 1.0e-10);
   } else {
      TEST_CLOSE(sarah_sigma_AA.real()(2,2), softsusy_sigma_AA, 1.0e-10);
   }

   // TEST_CLOSE(sarah_sigma_AA.real()(1,2), 0.0, 1.0e-10);
   // TEST_CLOSE(sarah_sigma_AA.real()(2,1), 0.0, 1.0e-10);
}

void compare_charged_higgs_self_energy(MssmSoftsusy s, MSSM m)
{
   const double mHpm = s.displayDrBarPars().mHpm;
   const double mWrun = s.displayMwRun();
   const double scale = s.displayMu();
   const DoubleVector Hpm(m.get_MassHpm());

   TEST_EQUALITY(s.displayMu(), m.displayMu());

   bool twisted;    // true if Hpm(1) == Hpm, false otherwise
   const double p = mHpm;

   if (Hpm(1) <= Hpm(2)) {
      TEST_CLOSE(mWrun, Hpm(1), 1.0e-10);
      TEST_CLOSE(mHpm , Hpm(2), 1.0e-10);
      twisted = false;
   } else {
      TEST_CLOSE(mWrun, Hpm(2), 1.0e-10);
      TEST_CLOSE(mHpm , Hpm(1), 1.0e-10);
      twisted = true;
   }

   const double softsusy_sigma_HpHm = s.piHpHm(mHpm, scale);

   ComplexMatrix sarah_sigma_Hpm(2,2);
   for (unsigned i1 = 1; i1 <= 2; ++i1) {
      for (unsigned i2 = 1; i2 <= 2; ++i2) {
         sarah_sigma_Hpm(i1,i2) = m.self_energy_Hpm(p,i1,i2);
      }
   }

   // do tree-level rotation to compare with SoftSusy
   const DoubleMatrix tree_tevel_rotation(m.get_ZP());
   sarah_sigma_Hpm = tree_tevel_rotation * sarah_sigma_Hpm
      * tree_tevel_rotation.transpose();

   TEST_CLOSE(sarah_sigma_Hpm.imag(), DoubleMatrix(2,2), 1.0e-10);
   TEST_CLOSE(sarah_sigma_Hpm.real(), sarah_sigma_Hpm.real().transpose(), 1.0e-10);

   if (twisted) {
      TEST_CLOSE(sarah_sigma_Hpm.real()(1,1), softsusy_sigma_HpHm, 1.0e-10);
   } else {
      TEST_CLOSE(sarah_sigma_Hpm.real()(2,2), softsusy_sigma_HpHm, 1.0e-10);
   }

   // TEST_CLOSE(sarah_sigma_Hpm.real()(1,2), 0.0, 1.0e-10);
   // TEST_CLOSE(sarah_sigma_Hpm.real()(2,1), 0.0, 1.0e-10);
}

void compare_z_self_energy(MssmSoftsusy s, MSSM m)
{
   const double p = m.get_MassVZ();
   const double scale = m.displayMu();
   const Complex sarah_z_se(m.self_energy_VZ(p));
   const double softsusy_z_se = s.piZZT(p, scale, false);

   TEST_EQUALITY(scale, m.displayMu());
   TEST_EQUALITY(scale, s.displayMu());
   TEST_EQUALITY(m.get_MassVZ(), s.displayMzRun());
   TEST_EQUALITY(sarah_z_se.imag(), 0.0);
   // Note: Softsusy uses on-shell masses for the 1st and 2nd
   // generation fermions.  FlexibleSUSY allways uses running DRbar
   // masses.  This leads to some deviation.
   TEST_CLOSE(sarah_z_se.real(), softsusy_z_se, 1.0e-4);
}

void compare_w_self_energy(MssmSoftsusy s, MSSM m)
{
   const double p = m.get_MassVWm();
   const double scale = m.displayMu();
   const Complex sarah_w_se(m.self_energy_VWm(p));
   const double softsusy_w_se = s.piWWT(p, scale, false);

   TEST_EQUALITY(scale, m.displayMu());
   TEST_EQUALITY(scale, s.displayMu());
   TEST_EQUALITY(m.get_MassVZ(), s.displayMzRun());
   TEST_EQUALITY(sarah_w_se.imag(), 0.0);
   // Note: Softsusy uses on-shell masses for the 1st and 2nd
   // generation fermions.  FlexibleSUSY allways uses running DRbar
   // masses.  This leads to some deviation.
   TEST_CLOSE(sarah_w_se.real(), softsusy_w_se, 2.0e-3);
}

void compare_top_self_energy(MssmSoftsusy s, MSSM m)
{
   // Note: in calcRunningMt() the running and the pole top mass are
   // used simultaneously, which leads to some deviations

   s.setMu(MZ);
   m.setMu(MZ);
   m.calculate_DRbar_parameters();

   const double mtpole = s.displayDataSet().displayPoleMt();
   const double softsusy_mtop = s.calcRunningMt();
   const double sarah_mtop = m.calculate_MassFu_DRbar_1loop(mtpole, 3);

   TEST_CLOSE(sarah_mtop, softsusy_mtop, 1.0e-10);
}

void compare_bot_self_energy(MssmSoftsusy s, MSSM m)
{
   // Note: in calcRunningMb() the running and the pole bottom mass
   // are used simultaneously, which leads to some deviations

   s.setMu(MZ);
   m.setMu(MZ);
   m.calculate_DRbar_parameters();

   const double mbpole = s.displayDataSet().displayPoleMb();
   const double softsusy_mbot = s.calcRunningMb();
   const double sarah_mbot = m.calculate_MassFd_DRbar_1loop(mbpole, 3);

   TEST_CLOSE(sarah_mbot, softsusy_mbot, 1.0e-10);
}

void compare_tau_self_energy(MssmSoftsusy s, MSSM m)
{
   s.setMu(MZ);
   m.setMu(MZ);
   m.calculate_DRbar_parameters();

   const double mtaupole = s.displayDataSet().displayPoleMtau();
   const double softsusy_mtau = s.calcRunningMtau();
   const double sarah_mtau = m.calculate_MassFe_DRbar_1loop(mtaupole, 3);

   TEST_CLOSE(sarah_mtau, softsusy_mtau, 1.0e-10);
}

void compare_self_energies(MssmSoftsusy s, MSSM m)
{
   ensure_tree_level_ewsb(m);
   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   TEST_EQUALITY(s.displayMu(), m.displayMu());

   compare_gluino_self_energy(s, m);
   compare_neutralino_self_energy(s, m);
   compare_chargino_self_energy(s, m);
   compare_sneutrino_self_energy(s, m);
   compare_selectron_self_energy(s, m);
   compare_sup_self_energy(s, m);
   compare_sdown_self_energy(s, m);
   compare_CP_even_higgs_self_energy(s, m);
   compare_CP_odd_higgs_self_energy(s, m);
   compare_charged_higgs_self_energy(s, m);
   compare_z_self_energy(s, m);
   compare_w_self_energy(s, m);
   compare_top_self_energy(s, m);
   compare_bot_self_energy(s, m);
   compare_tau_self_energy(s, m);
}

void compare_tadpoles(MssmSoftsusy s, MSSM m)
{
   ensure_tree_level_ewsb(m);
   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();
   const double vd = m.get_vd();
   const double vu = m.get_vu();

   const double td = m.tadpole_hh(1).real();
   const double tu = m.tadpole_hh(2).real();

   // check equality of tadpoles
   TEST_CLOSE(td / vd, s.doCalcTadpole1oneLoop(mt, sinthDRbar), 1.0e-11);
   TEST_CLOSE(tu / vu, s.doCalcTadpole2oneLoop(mt, sinthDRbar), 1.0e-11);
}

void compare_loop_masses(MssmSoftsusy s, MSSM m)
{
   ensure_tree_level_ewsb(m);
   softsusy::numHiggsMassLoops = 1;
   s.physical(1);
   m.calculate_1loop_masses();

   TEST_EQUALITY(s.displayMu(), m.displayMu());

   TEST_CLOSE(s.displayPhys().msnu, m.get_physical().MassSv, 1.0e-10);
   TEST_CLOSE(s.displayPhys().mGluino, m.get_physical().MassGlu, 1.0e-4);
   TEST_CLOSE(s.displayPhys().mneut.apply(fabs), m.get_physical().MassChi, 1.0e-10);
   TEST_CLOSE(s.displayPhys().mch.apply(fabs), m.get_physical().MassCha, 1.0e-10);

   TEST_CLOSE(s.displayPhys().me(1,1), m.get_physical().MassSe(1), 1.0e-10);
   TEST_CLOSE(s.displayPhys().me(1,2), m.get_physical().MassSe(2), 1.0e-10);
   TEST_CLOSE(s.displayPhys().me(1,3), m.get_physical().MassSe(3), 1.0e-10);
   TEST_CLOSE(s.displayPhys().me(2,1), m.get_physical().MassSe(4), 1.0e-10);
   TEST_CLOSE(s.displayPhys().me(2,2), m.get_physical().MassSe(5), 1.0e-10);
   TEST_CLOSE(s.displayPhys().me(2,3), m.get_physical().MassSe(6), 1.0e-10);

   TEST_CLOSE(s.displayPhys().mu(1,1), m.get_physical().MassSu(1), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mu(1,2), m.get_physical().MassSu(2), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mu(1,3), m.get_physical().MassSu(3), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mu(2,1), m.get_physical().MassSu(4), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mu(2,2), m.get_physical().MassSu(5), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mu(2,3), m.get_physical().MassSu(6), 1.0e-10);

   TEST_CLOSE(s.displayPhys().md(1,1), m.get_physical().MassSd(1), 1.0e-10);
   TEST_CLOSE(s.displayPhys().md(1,2), m.get_physical().MassSd(2), 1.0e-10);
   TEST_CLOSE(s.displayPhys().md(1,3), m.get_physical().MassSd(3), 1.0e-10);
   TEST_CLOSE(s.displayPhys().md(2,1), m.get_physical().MassSd(4), 1.0e-10);
   TEST_CLOSE(s.displayPhys().md(2,2), m.get_physical().MassSd(5), 1.0e-10);
   TEST_CLOSE(s.displayPhys().md(2,3), m.get_physical().MassSd(6), 1.0e-10);

   TEST_EQUALITY(0.0, m.get_physical().MassVG);
   TEST_EQUALITY(0.0, m.get_physical().MassVP);

   ensure_one_loop_ewsb(m);
   ensure_tree_level_ewsb(s);
   s.physical(1);
   m.calculate_1loop_masses();

   DoubleVector hh(m.get_physical().Masshh);
   if (hh(1) <= hh(2)) {
      TEST_CLOSE(s.displayPhys().mh0, m.get_physical().Masshh(1), 1.0e-10);
      TEST_CLOSE(s.displayPhys().mH0, m.get_physical().Masshh(2), 1.0e-10);
   } else {
      TEST_CLOSE(s.displayPhys().mh0, m.get_physical().Masshh(2), 1.0e-10);
      TEST_CLOSE(s.displayPhys().mH0, m.get_physical().Masshh(1), 1.0e-10);
   }
   TEST_CLOSE(m.get_physical().MassVZ , m.get_physical().MassAh(1), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mA0     , m.get_physical().MassAh(2), 1.0e-10);
   TEST_CLOSE(m.get_physical().MassVWm, m.get_physical().MassHpm(1), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mHpm    , m.get_physical().MassHpm(2), 1.0e-10);
}

void test_ewsb_tree(MSSM model, MssmSoftsusy softSusy)
{
   softSusy.calcDrBarPars();
   model.calculate_DRbar_parameters();

   const double BMu = model.get_BMu();
   const double Mu  = model.get_Mu();
   const double Mu2 = sqr(Mu);
   const double m1sq = model.get_mHd2();
   const double m2sq = -model.get_mHu2();
   const int signMu = Mu >= 0.0 ? 1 : -1;
   const DoubleVector pars(3); // unused
   const double precision = model.get_ewsb_iteration_precision();
   model.set_mHd2(m1sq);
   model.set_mHu2(m2sq);
   softSusy.setMh1Squared(m1sq);
   softSusy.setMh2Squared(m2sq);

   // these conditions must be fulfilled to have EWSB
   // see Drees p. 221 and 222
   TEST_GREATER(sqr(BMu), (m2sq + Mu2)*(m1sq + Mu2));
   // TEST_GREATER(m1sq + m2sq + 2*Mu2, 2*std::abs(BMu));

   // tree-level
   model.set_ewsb_loop_order(0);
   model.solve_ewsb();
   TEST_CLOSE(model.get_tadpole_vd(), 0.0, precision);
   TEST_CLOSE(model.get_tadpole_vu(), 0.0, precision);

   softSusy.rewsbTreeLevel(signMu);
   TEST_CLOSE(softSusy.displayM3Squared(), model.get_BMu(), 5.0);
   TEST_CLOSE(softSusy.displaySusyMu(), model.get_Mu(), 0.1);
}

void test_ewsb_1loop(MSSM model, MssmSoftsusy softSusy)
{
   softSusy.calcDrBarPars();
   model.calculate_DRbar_parameters();

   const double BMu = model.get_BMu();
   const double Mu  = model.get_Mu();
   const double Mu2 = sqr(Mu);
   const double m1sq = model.get_mHd2();
   const double m2sq = -model.get_mHu2();
   const int signMu = Mu >= 0.0 ? 1 : -1;
   const DoubleVector pars(3); // unused
   const double precision = model.get_ewsb_iteration_precision();
   model.set_mHd2(m1sq);
   model.set_mHu2(m2sq);
   softSusy.setMh1Squared(m1sq);
   softSusy.setMh2Squared(m2sq);

   // these conditions must be fulfilled to have EWSB
   // see Drees p. 221 and 222
   TEST_GREATER(sqr(BMu), (m2sq + Mu2)*(m1sq + Mu2));

   // one-loop
   model.set_ewsb_loop_order(1);
   model.solve_ewsb();
   TEST_CLOSE(model.get_tadpole_vd() - model.tadpole_hh(1).real(), 0.0, precision);
   TEST_CLOSE(model.get_tadpole_vu() - model.tadpole_hh(2).real(), 0.0, precision);

   softsusy::numRewsbLoops = 1;
   softSusy.rewsb(signMu, softSusy.displayDrBarPars().mt, pars);
   TEST_CLOSE(softSusy.displayM3Squared(), model.get_BMu(), 5.0);
   TEST_CLOSE(softSusy.displaySusyMu(), model.get_Mu(), 0.1);
}

void compare_models(int loopLevel)
{
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
   DoubleMatrix Yu(3,3), Yd(3,3), Ye(3,3);
   Yu(3,3) = 165.0   * root2 / (vev * sinBeta);
   Yd(3,3) = 2.9     * root2 / (vev * cosBeta);
   Ye(3,3) = 1.77699 * root2 / (vev * cosBeta);
   DoubleMatrix ID(3, 3), mm0(3, 3);
   for (int i=1; i<=3; i++) ID(i, i) = 1.0;
   mm0 = ID * sqr(m0);

   MSSM m;
   m.setMu(91);
   m.setLoops(loopLevel);
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

   MssmSoftsusy softSusy;
   softSusy.setMu(91);
   softSusy.setLoops(loopLevel);
   softSusy.setGaugeCoupling(1, g1);
   softSusy.setGaugeCoupling(2, g2);
   softSusy.setGaugeCoupling(3, g3);
   softSusy.setYukawaMatrix(YU, Yu);
   softSusy.setYukawaMatrix(YD, Yd);
   softSusy.setYukawaMatrix(YE, Ye);
   softSusy.setGauginoMass(1, M12);
   softSusy.setGauginoMass(2, M12);
   softSusy.setGauginoMass(3, M12);
   softSusy.setSoftMassMatrix(mQl, mm0);
   softSusy.setSoftMassMatrix(mUr, mm0);
   softSusy.setSoftMassMatrix(mDr, mm0);
   softSusy.setSoftMassMatrix(mLl, mm0);
   softSusy.setSoftMassMatrix(mEr, mm0);
   softSusy.setMh1Squared(sqr(m0));
   softSusy.setMh2Squared(sqr(m0));
   softSusy.setTrilinearMatrix(UA, a0 * Yu);
   softSusy.setTrilinearMatrix(DA, a0 * Yd);
   softSusy.setTrilinearMatrix(EA, a0 * Ye);
   softSusy.setSusyMu(susyMu);
   softSusy.setM3Squared(BMu);
   softSusy.setHvev(vev);
   softSusy.setTanb(tanBeta);

   std::cout << "comparing parameters ... ";
   test_parameter_equality(softSusy, m);
   std::cout << "done\n";
   std::cout << "comparing beta functions ... ";
   test_beta_function_equality(softSusy, m);
   std::cout << "done\n";
   std::cout << "comparing anomalous dimensions ... ";
   compare_anomalous_dimensions(softSusy, m);
   std::cout << "done\n";

   if (loopLevel == 1) {
      std::cout << "test tree-level ewsb ... ";
      test_ewsb_tree(m, softSusy);
      std::cout << "done\n";

      std::cout << "test one-loop ewsb ... ";
      test_ewsb_1loop(m, softSusy);
      std::cout << "done\n";

      std::cout << "comparing tree level masses ... ";
      compare_tree_level_masses(softSusy, m);
      std::cout << "done\n";

      std::cout << "comparing self-energies ... ";
      compare_self_energies(softSusy, m);
      std::cout << "done\n";

      std::cout << "comparing tadpoles ... ";
      compare_tadpoles(softSusy, m);
      std::cout << "done\n";

      std::cout << "comparing loop-masses ... ";
      compare_loop_masses(softSusy, m);
      std::cout << "done\n";

      std::cout << "testing high scale constraint ... ";
      test_high_scale_constraint(m);
      std::cout << "done\n";
   }
}

int main()
{
   std::cout << "====================\n";
   std::cout << "compare 1-loop level\n";
   std::cout << "====================\n";
   compare_models(1);

   std::cout << "====================\n";
   std::cout << "compare 2-loop level\n";
   std::cout << "====================\n";
   compare_models(2);

   return gErrors;
}
