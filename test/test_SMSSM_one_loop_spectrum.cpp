
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SMSSM.hpp"
#include "wrappers.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void ensure_tree_level_ewsb(SMSSM<Two_scale>& m, NmssmSoftsusy& s,
                            const SMSSM_input_parameters& input)
{
   const double precision = 1.0e-5;

   // initial guess
   m.set_Kappa(0.1);
   m.set_vS(5000.);
   m.set_ms2(-Sqr(input.m0));
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   s.setKappa(m.get_Kappa());
   s.setSvev(m.get_vS());
   s.setMsSquared(m.get_ms2());
   s.setMh1Squared(m.get_mHd2());
   s.setMh2Squared(m.get_mHu2());

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_CHECK_EQUAL(error, 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vd(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vu(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vS(), precision);

   softsusy::Z3 = false;
   s.rewsbTreeLevel(1);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Mu() , s.displaySusyMu()   , precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_BMu(), s.displayM3Squared(), precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_LL1(), s.displayXiS()      , precision * 5.);

   m.set_Mu(s.displaySusyMu());
   m.set_BMu(s.displayM3Squared());
   m.set_LL1(s.displayXiS());
}

void ensure_one_loop_ewsb(SMSSM<Two_scale>& m, NmssmSoftsusy& s)
{
   const double precision = 1.0e-5;

   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   const double mt = s.displayDrBarPars().mt;
   const int signMu = 1;

   m.set_ewsb_iteration_precision(precision);
   m.solve_ewsb_one_loop();

   softsusy::Z3 = false;
   softsusy::numRewsbLoops = 1;
   s.rewsb(signMu, mt);

   const double mu_ss  = s.displaySusyMu();
   const double bmu_ss = s.displayM3Squared();
   const double xis_ss = s.displayXiS();

   const double mu_fs  = m.get_Mu();
   const double bmu_fs = m.get_BMu();
   const double xis_fs = m.get_LL1();

   BOOST_CHECK_CLOSE_FRACTION(mu_ss , mu_fs , 5.0e-11);
   BOOST_CHECK_CLOSE_FRACTION(bmu_ss, bmu_fs, 2.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(xis_ss, xis_fs, 1.0e-10);

   m.set_Mu(s.displaySusyMu());
   m.set_BMu(s.displayM3Squared());
   m.set_LL1(s.displayXiS());
}

BOOST_AUTO_TEST_CASE( test_SMSSM_pole_masses )
{
   SMSSM_input_parameters input;
   input.m0       = 540.;
   input.Azero    = -350.;
   input.MSInput  = 290.;
   input.BMSInput = 400.;
   input.L1Input  = 300.;
   SMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   setup_SMSSM(m, s, input);

   softsusy::Z3 = false;
   ensure_tree_level_ewsb(m, s, input);

   softsusy::numHiggsMassLoops = 1;
   s.physical(1);
   m.calculate_DRbar_parameters();
   m.calculate_1loop_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print(ostr);
      BOOST_FAIL(ostr.str());
   }

   BOOST_CHECK_EQUAL(s.displayMu(), m.get_scale());

   // neutralinos
   const DoubleVector MChi(m.get_physical().MChi);
   const DoubleVector mneut = s.displayPhys().mneut.apply(fabs);
   BOOST_CHECK_CLOSE(mneut(1), MChi(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(2), MChi(2), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(3), MChi(3), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(4), MChi(4), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(5), MChi(5), 1.0e-12);

   // charginos
   const DoubleVector MCha(m.get_physical().MCha);
   const DoubleVector mch = s.displayPhys().mch.apply(fabs);
   BOOST_CHECK_CLOSE(mch(1), MCha(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mch(2), MCha(2), 1.0e-12);

   // photon, gluon mass
   const double vp = m.get_physical().MVP;
   const double vg = m.get_physical().MVG;
   BOOST_CHECK_EQUAL(vp, 0.0);
   BOOST_CHECK_EQUAL(vg, 0.0);

   // gluinos
   const double MGlu = m.get_physical().MGlu;
   const double mGluino = s.displayPhys().mGluino;
   BOOST_CHECK_CLOSE(MGlu, mGluino, 1.0e-04);

   // down-type squarks
   const DoubleVector Sd(m.get_physical().MSd);
   const DoubleVector md(s.displayPhys().md.flatten().sort());
   BOOST_CHECK_CLOSE(Sd(1), md(1), 0.01); // @todo increase precision
   BOOST_CHECK_CLOSE(Sd(2), md(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(3), md(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(4), md(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(5), md(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(6), md(6), 0.01); // @todo increase precision

   // up-type squarks
   const DoubleVector Su(m.get_physical().MSu);
   const DoubleVector mu(s.displayPhys().mu.flatten().sort());
   BOOST_CHECK_CLOSE(Su(1), mu(1), 0.8); // @todo increase precision
   BOOST_CHECK_CLOSE(Su(2), mu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(3), mu(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(4), mu(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(5), mu(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(6), mu(6), 0.5); // @todo increase precision

   // down-type sleptons
   const DoubleVector Se(m.get_physical().MSe);
   const DoubleVector me(s.displayPhys().me.flatten().sort());
   BOOST_CHECK_CLOSE(Se(1), me(1), 0.003); // @todo increase precision
   BOOST_CHECK_CLOSE(Se(2), me(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(3), me(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(4), me(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(5), me(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(6), me(6), 0.005); // @todo increase precision

   // up-type sleptons
   const DoubleVector Sv(m.get_physical().MSv);
   const DoubleVector msnu(s.displayPhys().msnu.sort());
   BOOST_CHECK_CLOSE(Sv(1), msnu(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Sv(2), msnu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sv(3), msnu(3), 1.0e-12);

   // neutrinos
   const DoubleVector MFv(m.get_physical().MFv);
   BOOST_CHECK_EQUAL(MFv(1), 0.0);
   BOOST_CHECK_EQUAL(MFv(2), 0.0);
   BOOST_CHECK_EQUAL(MFv(3), 0.0);

   // leptons
   const DoubleVector MFe(m.get_physical().MFe);
   BOOST_CHECK_EQUAL(MFe(1), 0.0);
   BOOST_CHECK_EQUAL(MFe(2), 0.0);
   // BOOST_CHECK_CLOSE(MFe(3), s.displayPhys().mtau, 1.0e-12);

   // ups
   const DoubleVector MFu(m.get_physical().MFu);
   BOOST_CHECK_EQUAL(MFu(1), 0.0);
   BOOST_CHECK_EQUAL(MFu(2), 0.0);
   // BOOST_CHECK_CLOSE(MFu(3), s.displayPhys().mt, 1.0e-12);

   // downs
   const DoubleVector MFd(m.get_physical().MFd);
   BOOST_CHECK_EQUAL(MFd(1), 0.0);
   BOOST_CHECK_EQUAL(MFd(2), 0.0);
   // BOOST_CHECK_CLOSE(MFd(3), s.displayPhys().mb, 1.0e-12);

   ensure_one_loop_ewsb(m, s);
   m.calculate_1loop_masses();
   s.physical(1);

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print(ostr);
      BOOST_FAIL(ostr.str());
   }

   // neutral CP even Higgs
   const DoubleVector hh(m.get_physical().Mhh);
   const DoubleVector mh0(s.displayPhys().mh0);
   BOOST_CHECK_CLOSE(hh(1), mh0(1), 0.0025);
   BOOST_CHECK_CLOSE(hh(2), mh0(2), 0.015);
   BOOST_CHECK_CLOSE(hh(3), mh0(3), 0.01);

   // neutral CP odd Higgs
   const DoubleVector Ah(m.get_physical().MAh);
   const DoubleVector mA0(s.displayPhys().mA0);
   BOOST_CHECK_CLOSE(Ah(2), mA0(1), 0.015);
   BOOST_CHECK_CLOSE(Ah(3), mA0(2), 0.008);

   // charged Higgs
   const DoubleVector Hpm(m.get_physical().MHpm);
   const double mHpm = s.displayPhys().mHpm;
   BOOST_CHECK_CLOSE(Hpm(2), mHpm, 0.02);
}
