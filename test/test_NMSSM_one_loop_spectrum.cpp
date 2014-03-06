
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "test_NMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void ensure_tree_level_ewsb(NMSSM<Two_scale>& m, NmssmSoftsusy& s,
                            const NMSSM_input_parameters& input)
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

   softsusy::Z3 = true;
   s.rewsbTreeLevel(1);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Kappa(), s.displayKappa()    , precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_vS()   , s.displaySvev()     , precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_ms2()  , s.displayMsSquared(), precision * 5.);

   m.set_Kappa(s.displayKappa());
   m.set_vS(s.displaySvev());
   m.set_ms2(s.displayMsSquared());
}

void ensure_one_loop_ewsb(NMSSM<Two_scale>& m, NmssmSoftsusy& s)
{
   const double precision = 1.0e-5;

   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   const double mt = s.displayDrBarPars().mt;
   const int signMu = 1;

   m.set_ewsb_iteration_precision(precision);
   m.solve_ewsb_one_loop();

   softsusy::numRewsbLoops = 1;
   s.rewsb(signMu, mt);

   const double kappa_ss = s.displayKappa();
   const double vS_ss    = s.displaySvev();
   const double ms2_ss   = s.displayMsSquared();

   const double kappa_fs = m.get_Kappa();
   const double vS_fs    = m.get_vS();
   const double ms2_fs   = m.get_ms2();

   BOOST_CHECK_CLOSE_FRACTION(kappa_ss, kappa_fs, 1.0e-11);
   BOOST_CHECK_CLOSE_FRACTION(vS_ss   , vS_fs   , 2.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(ms2_ss  , ms2_fs  , 2.0e-10);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_pole_masses )
{
   NMSSM_input_parameters input;
   input.m0 = 250.;
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   ensure_tree_level_ewsb(m, s, input);

   softsusy::numHiggsMassLoops = 1;
   s.physical(1);
   m.calculate_DRbar_parameters();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print(ostr);
      BOOST_FAIL(ostr.str());
   }

   BOOST_CHECK_EQUAL(s.displayMu(), m.get_scale());

   // neutralinos
   const DoubleVector MChi(ToDoubleVector(m.get_physical().MChi));
   const DoubleVector mneut = s.displayPhys().mneut.apply(fabs);
   BOOST_CHECK_CLOSE(mneut(1), MChi(1), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(2), MChi(2), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(3), MChi(3), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(4), MChi(4), 1.0e-12);
   BOOST_CHECK_CLOSE(mneut(5), MChi(5), 1.0e-12);

   // charginos
   const DoubleVector MCha(ToDoubleVector(m.get_physical().MCha));
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
   const DoubleVector Sd(ToDoubleVector(m.get_physical().MSd));
   const DoubleVector md(s.displayPhys().md.flatten().sort());
   BOOST_CHECK_CLOSE(Sd(1), md(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(2), md(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(3), md(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(4), md(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(5), md(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Sd(6), md(6), 1.0e-12);

   // up-type squarks
   const DoubleVector Su(ToDoubleVector(m.get_physical().MSu));
   const DoubleVector mu(s.displayPhys().mu.flatten().sort());
   BOOST_CHECK_CLOSE(Su(1), mu(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(2), mu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(3), mu(3), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(4), mu(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(5), mu(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Su(6), mu(6), 1.0e-12);

   // down-type sleptons
   const DoubleVector Se(ToDoubleVector(m.get_physical().MSe));
   const DoubleVector me(s.displayPhys().me.flatten().sort());
   BOOST_CHECK_CLOSE(Se(1), me(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(2), me(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(3), me(3), 4.0e-12);
   BOOST_CHECK_CLOSE(Se(4), me(4), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(5), me(5), 1.0e-12);
   BOOST_CHECK_CLOSE(Se(6), me(6), 1.0e-12);

   // up-type sleptons
   const DoubleVector Sv(ToDoubleVector(m.get_physical().MSv));
   const DoubleVector msnu(s.displayPhys().msnu.sort());
   BOOST_CHECK_CLOSE(Sv(1), msnu(1), 1.0e-12);
   BOOST_CHECK_CLOSE(Sv(2), msnu(2), 1.0e-12);
   BOOST_CHECK_CLOSE(Sv(3), msnu(3), 1.0e-12);

   // neutrinos
   const DoubleVector MFv(ToDoubleVector(m.get_physical().MFv));
   BOOST_CHECK_EQUAL(MFv(1), 0.0);
   BOOST_CHECK_EQUAL(MFv(2), 0.0);
   BOOST_CHECK_EQUAL(MFv(3), 0.0);

   // leptons
   const DoubleVector MFe(ToDoubleVector(m.get_physical().MFe));
   BOOST_CHECK_EQUAL(MFe(1), 0.0);
   BOOST_CHECK_EQUAL(MFe(2), 0.0);
   // BOOST_CHECK_CLOSE(MFe(3), s.displayPhys().mtau, 1.0e-12);

   // ups
   const DoubleVector MFu(ToDoubleVector(m.get_physical().MFu));
   BOOST_CHECK_EQUAL(MFu(1), 0.0);
   BOOST_CHECK_EQUAL(MFu(2), 0.0);
   // BOOST_CHECK_CLOSE(MFu(3), s.displayPhys().mt, 1.0e-12);

   // downs
   const DoubleVector MFd(ToDoubleVector(m.get_physical().MFd));
   BOOST_CHECK_EQUAL(MFd(1), 0.0);
   BOOST_CHECK_EQUAL(MFd(2), 0.0);
   // BOOST_CHECK_CLOSE(MFd(3), s.displayPhys().mb, 1.0e-12);

   ensure_one_loop_ewsb(m, s);
   m.calculate_pole_masses();
   s.physical(1);

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print(ostr);
      BOOST_FAIL(ostr.str());
   }

   // neutral CP even Higgs
   const DoubleVector hh(ToDoubleVector(m.get_physical().Mhh));
   const DoubleVector mh0(s.displayPhys().mh0);
   BOOST_CHECK_CLOSE(hh(1), mh0(1), 0.002);
   BOOST_CHECK_CLOSE(hh(2), mh0(2), 0.0001);
   BOOST_CHECK_CLOSE(hh(3), mh0(3), 0.1);

   // neutral CP odd Higgs
   const DoubleVector Ah(ToDoubleVector(m.get_physical().MAh));
   const DoubleVector mA0(s.displayPhys().mA0);
   BOOST_CHECK_CLOSE(Ah(2), mA0(1), 0.0002);
   BOOST_CHECK_CLOSE(Ah(3), mA0(2), 0.0006);

   // charged Higgs
   const DoubleVector Hpm(ToDoubleVector(m.get_physical().MHpm));
   const double mHpm = s.displayPhys().mHpm;
   BOOST_CHECK_CLOSE(Hpm(2), mHpm, 0.003);
}
