
#include "MSSM_two_scale_model.hpp"
#include "test.h"
#include "test_MSSM.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "conversion.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

void OrderAccordingTo(DoubleVector& m, DoubleMatrix& z, const DoubleMatrix& ref)
{
   const int cols = ref.displayCols();
   const int rows = ref.displayRows();
   const int size = rows * cols;

   if (cols != 3) {
      cout << "<OrderAccordingTo> Error: reference vector dose not have"
         " 2 columns" << endl;
      return;
   }
   if (rows != 2) {
      cout << "<OrderAccordingTo> Error: reference vector dose not have"
         " 3 rows" << endl;
      return;
   }
   if (m.displayStart() != 1) {
      cout << "<OrderAccordingTo> Error: mass vector dose not begin"
         " at index 1" << endl;
      return;
   }
   if (m.displayEnd() != size) {
      cout << "<OrderAccordingTo> Error: mass vector dose not end"
         " at index " << size << endl;
      return;
   }
   if (z.displayCols() != size || z.displayCols() != size) {
      cout << "<OrderAccordingTo> Error: mixing matrix dose not have"
         " " << size << " rows or cols" << endl;
      return;
   }

   for (int i = 1; i <= cols; i++) {
      for (int k = 1; k <= rows; k++) {
         const double m1 = ref(k,i);
         const int idx = m.closest(m1);
         m.swap(idx, (k-1) * cols + i);
         z.swaprows(idx, (k-1) * cols + i);
      }
   }
}

void compare_anomalous_dimensions(const SoftParsMssm& a, const MSSM_soft_parameters& b)
{
  DoubleMatrix gEE(3,3),gLL(3,3),gQQ(3,3),gDD(3,3),gUU(3,3);
  double gH1H1 = 0.0, gH2H2 = 0.0;
  DoubleVector dg(1,3);
  sBrevity brevity;
  a.anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, gH1H1, gH2H2, brevity);

  TEST_EQUALITY(a.displayLoops(), b.get_loops());
  TEST_EQUALITY(a.displayMu(), b.get_scale());
  TEST_EQUALITY(a.displayThresholds(), b.get_thresholds());

  TEST_EQUALITY(gEE, b.get_SeRSeR());
  TEST_EQUALITY(gLL, b.get_SlSl());
  TEST_EQUALITY(gQQ, b.get_SqSq());
  TEST_EQUALITY(gUU, b.get_SuRSuR());
  TEST_EQUALITY(gDD, b.get_SdRSdR());
  TEST_EQUALITY(gH1H1, b.get_SHdSHd());
  TEST_EQUALITY(gH2H2, b.get_SHuSHu());
}

void test_beta_function_equality(const SoftParsMssm& a, const MSSM_soft_parameters& b)
{
   SoftParsMssm beta_a(a.beta2());
   MSSM_soft_parameters beta_b(b.calc_beta());

   TEST_EQUALITY(beta_a.displayLoops(), beta_b.get_loops());
   TEST_EQUALITY(beta_a.displayMu(), beta_b.get_scale());
   TEST_EQUALITY(beta_a.displayThresholds(), beta_b.get_thresholds());

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

void compare_tree_level_masses(MssmSoftsusy s, MSSM<Two_scale> m)
{
   ensure_tree_level_ewsb(m);
   m.calculate_DRbar_parameters();
   s.calcDrBarPars();

   const double beta = atan(s.displayTanb());
   const double alpha = s.displayDrBarPars().thetaH;

   // check that tadpole eqs. are fulfilled
   TEST_CLOSE(m.get_ewsb_eq_vd(), 0.0, 1.0e-8);
   TEST_CLOSE(m.get_ewsb_eq_vu(), 0.0, 1.0e-9);

   // neutral CP even Higgs
   DoubleVector hh(ToDoubleVector(m.get_Mhh()));
   TEST_EQUALITY(hh, s.displayDrBarPars().mh0);
   TEST_CLOSE(m.get_ZH()(0,0), -sin(alpha), 1.0e-12);
   TEST_CLOSE(m.get_ZH()(0,1), cos(alpha), 1.0e-12);
   TEST_CLOSE(m.get_ZH()(1,0), cos(alpha), 1.0e-12);
   TEST_CLOSE(m.get_ZH()(1,1), sin(alpha), 1.0e-12);

   // neutral CP odd Higgs
   DoubleVector Ah(ToDoubleVector(m.get_MAh()));
   TEST_CLOSE(Ah(1), s.displayMzRun(), 1.0e-11);
   TEST_EQUALITY(Ah(2), s.displayDrBarPars().mA0(1));
   TEST_CLOSE(m.get_ZA()(0,0), -cos(beta), 1.0e-12);
   TEST_CLOSE(m.get_ZA()(0,1), sin(beta), 1.0e-12);
   TEST_CLOSE(m.get_ZA()(1,0), sin(beta), 1.0e-12);
   TEST_CLOSE(m.get_ZA()(1,1), cos(beta), 1.0e-12);

   // charged Higgs
   DoubleVector Hpm(ToDoubleVector(m.get_MHpm()));
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
   DoubleVector mneut(ToDoubleVector(m.get_MChi()));
   TEST_EQUALITY(mneut(1), s.displayDrBarPars().mnBpmz(1));
   TEST_EQUALITY(mneut(2), s.displayDrBarPars().mnBpmz(2));
   TEST_EQUALITY(mneut(3), s.displayDrBarPars().mnBpmz(3));
   TEST_EQUALITY(mneut(4), s.displayDrBarPars().mnBpmz(4));
   const ComplexMatrix ZN(ToComplexMatrix(m.get_ZN()));
   TEST_CLOSE(ZN, s.displayDrBarPars().nBpmz, 1.0e-12);
   // check neutralino diagonalization convention
   // diag = ZN^* M ZN^\dagger
   const DoubleMatrix m_chi(ToDoubleMatrix(m.get_mass_matrix_Chi()));
   TEST_CLOSE(DoubleMatrix(mneut), ZN.complexConjugate() * m_chi *
              ZN.hermitianConjugate(), 1.0e-10);

   // charginos
   DoubleVector mch(ToDoubleVector(m.get_MCha()));
   TEST_EQUALITY(mch(1), s.displayDrBarPars().mchBpmz(1));
   TEST_EQUALITY(mch(2), s.displayDrBarPars().mchBpmz(2));
   const ComplexMatrix UM(ToComplexMatrix(m.get_UM()) * Complex(0.,-1.));
   const ComplexMatrix UP(ToComplexMatrix(m.get_UP()) * Complex(0.,1.));
   TEST_CLOSE(UM, s.displayDrBarPars().uBpmz, 1.0e-12);
   TEST_CLOSE(UP, s.displayDrBarPars().vBpmz, 1.0e-12);

   // photon, W and Z mass
   const double vp = m.get_MVP();
   const double vz = m.get_MVZ();
   const double vw = m.get_MVWm();
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
   DoubleVector Sd(ToDoubleVector(m.get_MSd()));
   DoubleMatrix ZD(ToDoubleMatrix(m.get_ZD()));
   const DoubleMatrix md(s.displayDrBarPars().md);
   const double thetab = s.displayDrBarPars().thetab;
   OrderAccordingTo(Sd, ZD, md);
   TEST_EQUALITY(Sd(1), md(1,1));
   TEST_EQUALITY(Sd(2), md(1,2));
   TEST_EQUALITY(Sd(3), md(1,3));
   TEST_EQUALITY(Sd(4), md(2,1));
   TEST_EQUALITY(Sd(5), md(2,2));
   TEST_EQUALITY(Sd(6), md(2,3));
   // TEST_CLOSE(ZD(1,1), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(1,4), 0.0, 1.0e-12);
   TEST_CLOSE(ZD(4,1), 0.0, 1.0e-12);
   // TEST_CLOSE(ZD(4,4), 1.0, 1.0e-12);
   // TEST_CLOSE(ZD(2,2), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(2,5), 0.0, 1.0e-12);
   TEST_CLOSE(ZD(5,2), 0.0, 1.0e-12);
   // TEST_CLOSE(ZD(5,5), 1.0, 1.0e-12);
   TEST_CLOSE(ZD(3,3), cos(thetab), 1.0e-12);
   TEST_CLOSE(ZD(3,6), sin(thetab), 1.0e-12);
   TEST_CLOSE(ZD(6,3), -sin(thetab), 1.0e-12);
   TEST_CLOSE(ZD(6,6), cos(thetab), 1.0e-12);

   // up-type squarks
   DoubleVector Su(ToDoubleVector(m.get_MSu()));
   DoubleMatrix ZU(ToDoubleMatrix(m.get_ZU()));
   const DoubleMatrix mu(s.displayDrBarPars().mu);
   const double thetat = s.displayDrBarPars().thetat;
   OrderAccordingTo(Su, ZU, mu);
   TEST_EQUALITY(Su(1), mu(1,1));
   TEST_EQUALITY(Su(2), mu(1,2));
   TEST_EQUALITY(Su(3), mu(1,3));
   TEST_EQUALITY(Su(4), mu(2,1));
   TEST_EQUALITY(Su(5), mu(2,2));
   TEST_EQUALITY(Su(6), mu(2,3));
   // TEST_CLOSE(ZU(1,1), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(1,4), 0.0, 1.0e-12);
   TEST_CLOSE(ZU(4,1), 0.0, 1.0e-12);
   // TEST_CLOSE(ZU(4,4), 1.0, 1.0e-12);
   // TEST_CLOSE(ZU(2,2), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(2,5), 0.0, 1.0e-12);
   TEST_CLOSE(ZU(5,2), 0.0, 1.0e-12);
   // TEST_CLOSE(ZU(5,5), 1.0, 1.0e-12);
   TEST_CLOSE(ZU(3,3), -cos(thetat), 1.0e-12);
   TEST_CLOSE(ZU(3,6), -sin(thetat), 1.0e-12);
   TEST_CLOSE(ZU(6,3), -sin(thetat), 1.0e-12);
   TEST_CLOSE(ZU(6,6), cos(thetat), 1.0e-12);

   // sleptons
   DoubleVector Se(ToDoubleVector(m.get_MSe()));
   DoubleMatrix ZE(ToDoubleMatrix(m.get_ZE()));
   const DoubleMatrix me(s.displayDrBarPars().me);
   const double thetatau = s.displayDrBarPars().thetatau;
   OrderAccordingTo(Se, ZE, me);
   TEST_EQUALITY(Se(1), me(1,1));
   TEST_EQUALITY(Se(2), me(1,2));
   TEST_EQUALITY(Se(3), me(1,3));
   TEST_EQUALITY(Se(4), me(2,1));
   TEST_EQUALITY(Se(5), me(2,2));
   TEST_EQUALITY(Se(6), me(2,3));
   // TEST_CLOSE(ZE(1,1), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(1,4), 0.0, 1.0e-12);
   TEST_CLOSE(ZE(4,1), 0.0, 1.0e-12);
   // TEST_CLOSE(ZE(4,4), 1.0, 1.0e-12);
   // TEST_CLOSE(ZE(2,2), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(2,5), 0.0, 1.0e-12);
   TEST_CLOSE(ZE(5,2), 0.0, 1.0e-12);
   // TEST_CLOSE(ZE(5,5), 1.0, 1.0e-12);
   TEST_CLOSE(ZE(3,3), cos(thetatau), 1.0e-12);
   TEST_CLOSE(ZE(3,6), sin(thetatau), 1.0e-12);
   TEST_CLOSE(ZE(6,3), -sin(thetatau), 1.0e-12);
   TEST_CLOSE(ZE(6,6), cos(thetatau), 1.0e-12);

   // sneutrinos
   DoubleVector msnu(s.displayDrBarPars().msnu);
   DoubleVector Snu(ToDoubleVector(m.get_MSv()));

   TEST_EQUALITY(Snu(1), msnu(1));
   TEST_EQUALITY(Snu(2), msnu(2));
   TEST_EQUALITY(Snu(3), msnu(3));

   // gluons
   TEST_EQUALITY(m.get_MVG(), 0.0);

   // gluinos
   TEST_EQUALITY(m.get_MGlu(), s.displayDrBarPars().mGluino);

   // neutrinos
   TEST_EQUALITY(m.get_MFv()(0), 0.0);
   TEST_EQUALITY(m.get_MFv()(1), 0.0);
   TEST_EQUALITY(m.get_MFv()(2), 0.0);

   // leptons
   TEST_EQUALITY(m.get_MFe()(0), 0.0);
   TEST_EQUALITY(m.get_MFe()(1), 0.0);
   TEST_EQUALITY(m.get_MFe()(2), s.displayDrBarPars().mtau);
   DoubleMatrix unity(3,3);
   unity(1,1) = 1.0;
   unity(2,2) = 1.0;
   unity(3,3) = 1.0;
   TEST_EQUALITY(ToComplexMatrix(m.get_ZEL()), unity);
   TEST_EQUALITY(ToComplexMatrix(m.get_ZER()), unity);

   // ups
   TEST_EQUALITY(m.get_MFu()(0), 0.0);
   TEST_EQUALITY(m.get_MFu()(1), 0.0);
   TEST_EQUALITY(m.get_MFu()(2), s.displayDrBarPars().mt);
   TEST_EQUALITY(ToComplexMatrix(m.get_ZUL()), unity);
   TEST_EQUALITY(ToComplexMatrix(m.get_ZUR()), unity);

   // downs
   TEST_EQUALITY(m.get_MFd()(0), 0.0);
   TEST_EQUALITY(m.get_MFd()(1), 0.0);
   TEST_EQUALITY(m.get_MFd()(2), s.displayDrBarPars().mb);
   TEST_EQUALITY(ToComplexMatrix(m.get_ZDL()), unity);
   TEST_EQUALITY(ToComplexMatrix(m.get_ZDR()), unity);
}

void compare_gluino_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   // tree-level
   s.gluino(0);
   const double m3 = s.displayGaugino(3);
   const double g3 = s.displayGaugeCoupling(3);
   const double softsusy_gluino_tree = s.displayPhys().mGluino;

   TEST_CLOSE(softsusy_gluino_tree, m3, 1.0e-12);
   TEST_CLOSE(softsusy_gluino_tree, m.get_MGlu(), 1.0e-12);
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

void compare_neutralino_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
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
         sarah_sigma_L(i,k) = m.self_energy_Chi_PL(p,i-1,k-1);
         sarah_sigma_R(i,k) = m.self_energy_Chi_PR(p,i-1,k-1);
         sarah_sigma_S(i,k) = m.self_energy_Chi_1(p,i-1,k-1);
      }
   }

   DoubleVector Chi(ToDoubleVector(m.get_MChi()));
   ComplexMatrix M_tree(ToComplexMatrix(m.get_mass_matrix_Chi()));

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

void compare_chargino_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
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
         sarah_sigma_L(i,k) = m.self_energy_Cha_PL(p,i-1,k-1);
         sarah_sigma_R(i,k) = m.self_energy_Cha_PR(p,i-1,k-1);
         sarah_sigma_S(i,k) = m.self_energy_Cha_1(p,i-1,k-1);
      }
   }

   DoubleVector Cha(ToDoubleVector(m.get_MCha()));
   ComplexMatrix M_tree(ToComplexMatrix(m.get_mass_matrix_Cha()));

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

void compare_sneutrino_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   // tree-level
   s.doSnu(0.0, 0);
   const DoubleVector Snu_softsusy_tree(s.displayPhys().msnu);
   DoubleVector Snu_sarah_tree(ToDoubleVector(m.get_MSv()));
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
         Snu_sarah_se(i1, i2) = m.self_energy_Sv(p, i1-1, i2-1);
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

void compare_selectron_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double mtau = s.displayDrBarPars().mtau;
   const double sinthDRbar = s.calcSinthdrbar();

   // tree-level
   s.doChargedSleptons(mtau, 0.0, sinthDRbar, 0);
   const DoubleMatrix Se_softsusy_tree(s.displayPhys().me);
   DoubleVector Se_sarah_tree(ToDoubleVector(m.get_MSe()));
   DoubleMatrix ZE(ToDoubleMatrix(m.get_ZE()));
   OrderAccordingTo(Se_sarah_tree, ZE, Se_softsusy_tree);
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
         Se_sarah_se(i1, i2) = m.self_energy_Se(p, i1-1, i2-1);
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
}

void compare_sup_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();

   // tree-level
   s.doUpSquarks(mt, 0.0, sinthDRbar, 0);
   const DoubleMatrix Su_softsusy_tree(s.displayPhys().mu);
   DoubleVector Su_sarah_tree(ToDoubleVector(m.get_MSu()));
   DoubleMatrix ZU(ToDoubleMatrix(m.get_ZU()));
   OrderAccordingTo(Su_sarah_tree, ZU, Su_softsusy_tree);
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
         Su_sarah_se(i1, i2) = m.self_energy_Su(p, i1-1, i2-1);
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

      Su_sarah_se(3,3) = m.self_energy_Su(p,2,2);
      Su_sarah_se(3,6) = m.self_energy_Su(p,2,5);
      Su_sarah_se(6,3) = m.self_energy_Su(p,5,2);
      Su_sarah_se(6,6) = m.self_energy_Su(p,5,5);
   }

   TEST_CLOSE(Su_softsusy_se(3,3), Su_sarah_se(3,3), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(3,6), Su_sarah_se(3,6), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(6,3), Su_sarah_se(6,3), 1.0e-10);
   TEST_CLOSE(Su_softsusy_se(6,6), Su_sarah_se(6,6), 1.0e-9);
}

void compare_sdown_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double mt = s.displayDrBarPars().mt;
   const double mb = s.displayDrBarPars().mb;
   const double sinthDRbar = s.calcSinthdrbar();

   // tree-level
   s.doDownSquarks(mb, 0.0, sinthDRbar, 0, mt);
   const DoubleMatrix Sd_softsusy_tree(s.displayPhys().md);
   DoubleVector Sd_sarah_tree(ToDoubleVector(m.get_MSd()));
   DoubleMatrix ZD(ToDoubleMatrix(m.get_ZD()));
   OrderAccordingTo(Sd_sarah_tree, ZD, Sd_softsusy_tree);
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
         Sd_sarah_se(i1, i2) = m.self_energy_Sd(p, i1-1, i2-1);
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

      Sd_sarah_se(3,3) = m.self_energy_Sd(p,2,2);
      Sd_sarah_se(3,6) = m.self_energy_Sd(p,2,5);
      Sd_sarah_se(6,3) = m.self_energy_Sd(p,5,2);
      Sd_sarah_se(6,6) = m.self_energy_Sd(p,5,5);
   }

   TEST_CLOSE(Sd_softsusy_se(3,3), Sd_sarah_se(3,3), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(3,6), Sd_sarah_se(3,6), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(6,3), Sd_sarah_se(6,3), 1.0e-10);
   TEST_CLOSE(Sd_softsusy_se(6,6), Sd_sarah_se(6,6), 1.2e-10);
}

void compare_CP_even_higgs_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double mh0 = s.displayDrBarPars().mh0(1);
   const double mH0 = s.displayDrBarPars().mh0(2);
   const double scale = s.displayMu();
   const DoubleVector hh(ToDoubleVector(m.get_Mhh()));

   TEST_EQUALITY(s.displayMu(), m.get_scale());

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
         sarah_sigma_light(i1,i2) = m.self_energy_hh(mh0,i1-1,i2-1);
         sarah_sigma_heavy(i1,i2) = m.self_energy_hh(mH0,i1-1,i2-1);
      }
   }

   TEST_CLOSE(sarah_sigma_light.imag(), DoubleMatrix(2,2), 1.0e-10);
   TEST_CLOSE(sarah_sigma_light.real(), sarah_sigma_light.real().transpose(), 1.0e-10);
   TEST_CLOSE(sarah_sigma_light.real(), softsusy_sigma_light, 1.0e-10);

   TEST_CLOSE(sarah_sigma_heavy.imag(), DoubleMatrix(2,2), 1.0e-10);
   TEST_CLOSE(sarah_sigma_heavy.real(), sarah_sigma_heavy.real().transpose(), 1.0e-10);
   TEST_CLOSE(sarah_sigma_heavy.real(), softsusy_sigma_heavy, 1.0e-10);
}

void compare_CP_odd_higgs_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double mA0 = s.displayDrBarPars().mA0(1);
   const double mZrun = s.displayMzRun();
   const double scale = s.displayMu();
   const DoubleVector Ah(ToDoubleVector(m.get_MAh()));

   TEST_EQUALITY(s.displayMu(), m.get_scale());

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
         sarah_sigma_AA(i1,i2) = m.self_energy_Ah(p,i1-1,i2-1);
      }
   }

   // do tree-level rotation to compare with SoftSusy
   const DoubleMatrix tree_tevel_rotation(ToDoubleMatrix(m.get_ZA()));
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

void compare_charged_higgs_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double mHpm = s.displayDrBarPars().mHpm;
   const double mWrun = s.displayMwRun();
   const double scale = s.displayMu();
   const DoubleVector Hpm(ToDoubleVector(m.get_MHpm()));

   TEST_EQUALITY(s.displayMu(), m.get_scale());

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
         sarah_sigma_Hpm(i1,i2) = m.self_energy_Hpm(p,i1-1,i2-1);
      }
   }

   // do tree-level rotation to compare with SoftSusy
   const DoubleMatrix tree_tevel_rotation(ToDoubleMatrix(m.get_ZP()));
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

void compare_z_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double p = m.get_MVZ();
   const double scale = m.get_scale();
   const Complex sarah_z_se(m.self_energy_VZ(p));
   const double softsusy_z_se = s.piZZT(p, scale, false);

   TEST_EQUALITY(scale, m.get_scale());
   TEST_EQUALITY(scale, s.displayMu());
   TEST_EQUALITY(m.get_MVZ(), s.displayMzRun());
   TEST_EQUALITY(sarah_z_se.imag(), 0.0);
   // Note: Softsusy uses on-shell masses for the 1st and 2nd
   // generation fermions.  FlexibleSUSY allways uses running DRbar
   // masses.  This leads to some deviation.
   TEST_CLOSE(sarah_z_se.real(), softsusy_z_se, 1.0e-4);
}

void compare_w_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   const double p = m.get_MVWm();
   const double scale = m.get_scale();
   const Complex sarah_w_se(m.self_energy_VWm(p));
   const double softsusy_w_se = s.piWWT(p, scale, false);

   TEST_EQUALITY(scale, m.get_scale());
   TEST_EQUALITY(scale, s.displayMu());
   TEST_EQUALITY(m.get_MVZ(), s.displayMzRun());
   TEST_EQUALITY(sarah_w_se.imag(), 0.0);
   // Note: Softsusy uses on-shell masses for the 1st and 2nd
   // generation fermions.  FlexibleSUSY allways uses running DRbar
   // masses.  This leads to some deviation.
   TEST_CLOSE(sarah_w_se.real(), softsusy_w_se, 2.0e-3);
}

void compare_top_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   // Note: in calcRunningMt() the running and the pole top mass are
   // used simultaneously, which leads to some deviations

   s.setMu(MZ);
   s.calcDrBarPars();
   m.set_scale(MZ);
   m.calculate_DRbar_parameters();

   const double mtpole = s.displayDataSet().displayPoleMt();
   const double softsusy_mtop = s.calcRunningMt();
   const double sarah_mtop = m.calculate_MFu_DRbar_1loop(mtpole, 2);

   TEST_CLOSE(sarah_mtop, softsusy_mtop, 0.14);
}

void compare_bot_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   // Note: in calcRunningMb() the running and the pole bottom mass
   // are used simultaneously, which leads to some deviations

   s.setMu(MZ);
   s.calcDrBarPars();
   m.set_scale(MZ);
   m.calculate_DRbar_parameters();

   const double mb_ms_bar = s.displayDataSet().displayMass(mBottom);
   const double softsusy_mbot = s.calcRunningMb();
   const double sarah_mbot = m.calculate_MFd_DRbar_1loop(mb_ms_bar, 2);

   TEST_CLOSE(sarah_mbot, softsusy_mbot, 0.0013);
}

void compare_tau_self_energy(MssmSoftsusy s, MSSM<Two_scale> m)
{
   s.setMu(MZ);
   s.calcDrBarPars();
   m.set_scale(MZ);
   m.calculate_DRbar_parameters();

   const double mtaupole = s.displayDataSet().displayPoleMtau();
   const double softsusy_mtau = s.calcRunningMtau();
   const double sarah_mtau = m.calculate_MFe_DRbar_1loop(mtaupole, 2);

   // Softsusy:
   // * no photon contribution
   // * sometimes on-shell masses in the loop functions
   // FlexibleSUSY:
   // * allways running DRbar masses in the loop functions
   TEST_CLOSE(sarah_mtau, softsusy_mtau, 0.04);
}

void compare_self_energies(MssmSoftsusy s, MSSM<Two_scale> m)
{
   ensure_tree_level_ewsb(m);
   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   TEST_EQUALITY(s.displayMu(), m.get_scale());

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

void compare_tadpoles(MssmSoftsusy s, MSSM<Two_scale> m)
{
   ensure_tree_level_ewsb(m);
   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();
   const double vd = m.get_vd();
   const double vu = m.get_vu();

   const double td = m.tadpole_hh(0).real();
   const double tu = m.tadpole_hh(1).real();

   // check equality of tadpoles
   TEST_CLOSE(td / vd, s.doCalcTadpole1oneLoop(mt, sinthDRbar), 1.0e-11);
   TEST_CLOSE(tu / vu, s.doCalcTadpole2oneLoop(mt, sinthDRbar), 1.0e-11);
}

void compare_loop_masses(MssmSoftsusy s, MSSM<Two_scale> m)
{
   ensure_tree_level_ewsb(m);
   ensure_tree_level_ewsb(s);
   softsusy::numHiggsMassLoops = 1;
   s.physical(1);
   m.calculate_DRbar_parameters();
   m.calculate_pole_masses();

   TEST_EQUALITY(s.displayMu(), m.get_scale());

   TEST_CLOSE(s.displayPhys().msnu             , ToDoubleVector(m.get_physical().MSv), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mGluino          , m.get_physical().MGlu, 1.0e-4);
   TEST_CLOSE(s.displayPhys().mneut.apply(fabs), ToDoubleVector(m.get_physical().MChi), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mch.apply(fabs)  , ToDoubleVector(m.get_physical().MCha), 1.0e-11);

   TEST_CLOSE(s.displayPhys().me.flatten().sort(), ToDoubleVector(m.get_physical().MSe), 1.0e-10);
   TEST_CLOSE(s.displayPhys().mu.flatten().sort(), ToDoubleVector(m.get_physical().MSu), 1.0e-10);
   TEST_CLOSE(s.displayPhys().md.flatten().sort(), ToDoubleVector(m.get_physical().MSd), 1.0e-10);

   TEST_EQUALITY(0.0, m.get_physical().MVG);
   TEST_EQUALITY(0.0, m.get_physical().MVP);

   // ensure that the important scalar potential parameters are equal
   // before solvin the 1-loop EWSB
   TEST_CLOSE(m.get_mHd2(), s.displayMh1Squared(), 1.0e-10);
   TEST_CLOSE(m.get_mHu2(), s.displayMh2Squared(), 1.0e-10);
   TEST_CLOSE(m.get_Mu(), s.displaySusyMu(), 1.0e-10);
   TEST_CLOSE(m.get_BMu(), s.displayM3Squared(), 1.0e-10);

   ensure_one_loop_ewsb(m);
   ensure_one_loop_ewsb(s);
   // check that the important scalar potential parameters are equal
   TEST_CLOSE(m.get_mHd2(), s.displayMh1Squared(), 1.0e-10);
   TEST_CLOSE(m.get_mHu2(), s.displayMh2Squared(), 1.0e-10);
   TEST_CLOSE_REL(m.get_Mu(), s.displaySusyMu(), 0.0007);
   TEST_CLOSE_REL(m.get_BMu(), s.displayM3Squared(), 0.0003);
   m.calculate_pole_masses();
   s.physical(1);

   const DoubleVector hh(ToDoubleVector(m.get_physical().Mhh));
   TEST_CLOSE(s.displayPhys().mh0(1), hh(1), 0.114);
   TEST_CLOSE(s.displayPhys().mh0(2), hh(2), 0.04);
   TEST_CLOSE_REL(s.displayPhys().mh0(1), hh(1), 0.00115);
   TEST_CLOSE_REL(s.displayPhys().mh0(2), hh(2), 6.0e-5);

   const DoubleVector Ah(ToDoubleVector(m.get_physical().MAh));
   TEST_CLOSE(s.displayPhys().mA0(1), Ah(2), 0.05);
   TEST_CLOSE_REL(s.displayPhys().mA0(1), Ah(2), 6.0e-5);

   const DoubleVector Hpm(ToDoubleVector(m.get_physical().MHpm));
   TEST_CLOSE(s.displayPhys().mHpm, Hpm(2), 0.09);
   TEST_CLOSE_REL(s.displayPhys().mHpm, Hpm(2), 1.3e-4);

   TEST_EQUALITY(0.0, m.get_physical().MVG);
   TEST_EQUALITY(0.0, m.get_physical().MVP);

   TEST_CLOSE(s.displayPhys().mGluino, m.get_physical().MGlu, 1.0e-4);
   TEST_CLOSE_REL(s.displayPhys().msnu(1), m.get_physical().MSv(0), 0.00013);
   TEST_CLOSE_REL(s.displayPhys().msnu(2), m.get_physical().MSv(1), 0.00013);
   TEST_CLOSE_REL(s.displayPhys().msnu(3), m.get_physical().MSv(2), 0.00014);
   TEST_CLOSE_REL(s.displayPhys().mneut.apply(fabs)(1), m.get_physical().MChi(0), 0.0012);
   TEST_CLOSE_REL(s.displayPhys().mneut.apply(fabs)(2), m.get_physical().MChi(1), 0.0004);
   TEST_CLOSE_REL(s.displayPhys().mneut.apply(fabs)(3), m.get_physical().MChi(2), 0.0017);
   TEST_CLOSE_REL(s.displayPhys().mneut.apply(fabs)(4), m.get_physical().MChi(3), 0.0005);
   TEST_CLOSE_REL(s.displayPhys().mch.apply(fabs)(1), m.get_physical().MCha(0), 0.0011);
   TEST_CLOSE_REL(s.displayPhys().mch.apply(fabs)(2), m.get_physical().MCha(1), 0.0008);

   TEST_CLOSE_REL(s.displayPhys().me.flatten().sort()(1), m.get_physical().MSe(0), 0.00011);
   TEST_CLOSE_REL(s.displayPhys().me.flatten().sort()(2), m.get_physical().MSe(1), 0.0002);
   TEST_CLOSE_REL(s.displayPhys().me.flatten().sort()(3), m.get_physical().MSe(2), 0.0002);
   TEST_CLOSE_REL(s.displayPhys().me.flatten().sort()(4), m.get_physical().MSe(3), 0.00012);
   TEST_CLOSE_REL(s.displayPhys().me.flatten().sort()(5), m.get_physical().MSe(4), 0.00012);
   TEST_CLOSE_REL(s.displayPhys().me.flatten().sort()(6), m.get_physical().MSe(5), 0.0001);

   TEST_CLOSE_REL(s.displayPhys().mu.flatten().sort()(1), m.get_physical().MSu(0), 0.0001);
   TEST_CLOSE_REL(s.displayPhys().mu.flatten().sort()(2), m.get_physical().MSu(1), 0.0001);
   TEST_CLOSE_REL(s.displayPhys().mu.flatten().sort()(3), m.get_physical().MSu(2), 0.00014);
   TEST_CLOSE_REL(s.displayPhys().mu.flatten().sort()(4), m.get_physical().MSu(3), 0.00014);
   TEST_CLOSE_REL(s.displayPhys().mu.flatten().sort()(5), m.get_physical().MSu(4), 0.00025);
   TEST_CLOSE_REL(s.displayPhys().mu.flatten().sort()(6), m.get_physical().MSu(5), 0.00044);

   TEST_CLOSE_REL(s.displayPhys().md.flatten().sort()(1), m.get_physical().MSd(0), 0.00003);
   TEST_CLOSE_REL(s.displayPhys().md.flatten().sort()(2), m.get_physical().MSd(1), 0.00007);
   TEST_CLOSE_REL(s.displayPhys().md.flatten().sort()(3), m.get_physical().MSd(2), 0.00007);
   TEST_CLOSE_REL(s.displayPhys().md.flatten().sort()(4), m.get_physical().MSd(3), 0.000026);
   TEST_CLOSE_REL(s.displayPhys().md.flatten().sort()(5), m.get_physical().MSd(4), 0.000026);
   TEST_CLOSE_REL(s.displayPhys().md.flatten().sort()(6), m.get_physical().MSd(5), 0.00005);
}

void test_ewsb_tree(MSSM<Two_scale> model, MssmSoftsusy softSusy)
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
   model.solve_ewsb_tree_level();
   TEST_CLOSE(model.get_ewsb_eq_vd(), 0.0, precision);
   TEST_CLOSE(model.get_ewsb_eq_vu(), 0.0, precision);

   softSusy.rewsbTreeLevel(signMu);
   TEST_CLOSE(softSusy.displayM3Squared(), model.get_BMu(), 5.0);
   TEST_CLOSE(softSusy.displaySusyMu(), model.get_Mu(), 0.1);
}

void test_ewsb_1loop(MSSM<Two_scale> model, MssmSoftsusy softSusy)
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
   model.solve_ewsb_one_loop();
   TEST_CLOSE(model.get_ewsb_eq_vd() - model.tadpole_hh(0).real(), 0.0, precision);
   TEST_CLOSE(model.get_ewsb_eq_vu() - model.tadpole_hh(1).real(), 0.0, precision);

   softsusy::numRewsbLoops = 1;
   softSusy.rewsb(signMu, softSusy.displayDrBarPars().mt, pars);
   TEST_CLOSE_REL(softSusy.displayM3Squared(), model.get_BMu(), 0.0004);
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

   MSSM<Two_scale> m;
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

   MssmSoftsusy softSusy;
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
