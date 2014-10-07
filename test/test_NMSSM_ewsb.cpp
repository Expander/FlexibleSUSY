
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_ewsb

#include <boost/test/unit_test.hpp>

#define private public

#include "test_NMSSM.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_NMSSM_ewsb_tree_level )
{
   NMSSM_input_parameters input;
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   const double precision = 1.0e-5;
   setup_NMSSM(m, s, input);

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
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3(), precision);

   softsusy::Z3 = true;
   s.rewsbTreeLevel(1);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Kappa(), s.displayKappa()    , precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_vS()   , s.displaySvev()     , precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_ms2()  , s.displayMsSquared(), precision * 5.);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_ewsb_tree_level_via_soft_higgs_masses )
{
   NMSSM_input_parameters input;
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   const int error = m.solve_ewsb_tree_level_via_soft_higgs_masses();

   BOOST_CHECK_EQUAL(error, 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1(), 5.0e-09);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2(), 1.0e-09);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3(), 1.0e-09);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_one_loop_tadpoles )
{
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   NMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   setup_NMSSM_const(m, s, input);

   s.calcDrBarPars();
   m.calculate_DRbar_masses();

   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();
   const double vd = m.get_vd();
   const double vu = m.get_vu();
   const double vS = m.get_vS();

   const Complex tadpole_hh_1(m.tadpole_hh(0));
   const Complex tadpole_hh_2(m.tadpole_hh(1));
   const Complex tadpole_hh_3(m.tadpole_hh(2));

   const double tadpole_ss_1 = s.doCalcTadpole1oneLoop(mt, sinthDRbar);
   const double tadpole_ss_2 = s.doCalcTadpole2oneLoop(mt, sinthDRbar);
   const double tadpole_ss_3 = s.doCalcTadpoleSoneLoop(mt, sinthDRbar);

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);

   BOOST_CHECK_CLOSE(Re(tadpole_hh_1) / vd, tadpole_ss_1, 1.0e-11);
   BOOST_CHECK_CLOSE(Re(tadpole_hh_2) / vu, tadpole_ss_2, 1.0e-11);
   BOOST_CHECK_CLOSE(Re(tadpole_hh_3) / vS, tadpole_ss_3, 1.0e-11);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_one_loop_ewsb )
{
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   NMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   const double precision = 1.0e-5;
   setup_NMSSM_const(m, s, input);

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

   s.calcDrBarPars();
   m.calculate_DRbar_masses();

   const double mt = s.displayDrBarPars().mt;
   const int signMu = 1;

   m.set_ewsb_iteration_precision(precision);
   m.solve_ewsb_one_loop();

   const Complex tadpole_hh_1(m.tadpole_hh(0));
   const Complex tadpole_hh_2(m.tadpole_hh(1));
   const Complex tadpole_hh_3(m.tadpole_hh(2));

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);

   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1() - Re(tadpole_hh_1), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2() - Re(tadpole_hh_2), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3() - Re(tadpole_hh_3), precision);

   softsusy::numRewsbLoops = 1;
   s.rewsb(signMu, mt);

   const double kappa_ss = s.displayKappa();
   const double vS_ss    = s.displaySvev();
   const double ms2_ss   = s.displayMsSquared();

   const double kappa_fs = m.get_Kappa();
   const double vS_fs    = m.get_vS();
   const double ms2_fs   = m.get_ms2();

   BOOST_CHECK_CLOSE_FRACTION(kappa_ss, kappa_fs, 1.0e-11);
   BOOST_CHECK_CLOSE_FRACTION(vS_ss   , vS_fs   , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(ms2_ss  , ms2_fs  , 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_two_loop_tadpoles )
{
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   NMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   setup_NMSSM_const(m, s, input);

   s.calcDrBarPars();
   m.calculate_DRbar_masses();

   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();
   const double vd = m.get_vd();
   const double vu = m.get_vu();
   const double vS = m.get_vS();

   // compare one-loop tadpoles
   const double tadpole_hh_1(Re(m.tadpole_hh(0)));
   const double tadpole_hh_2(Re(m.tadpole_hh(1)));
   const double tadpole_hh_3(Re(m.tadpole_hh(2)));

   const double tadpole_ss_1 = s.doCalcTadpole1oneLoop(mt, sinthDRbar);
   const double tadpole_ss_2 = s.doCalcTadpole2oneLoop(mt, sinthDRbar);
   const double tadpole_ss_3 = s.doCalcTadpoleSoneLoop(mt, sinthDRbar);

   BOOST_CHECK_CLOSE(tadpole_hh_1 / vd, tadpole_ss_1, 1.0e-11);
   BOOST_CHECK_CLOSE(tadpole_hh_2 / vu, tadpole_ss_2, 1.0e-11);
   BOOST_CHECK_CLOSE(tadpole_hh_3 / vS, tadpole_ss_3, 1.0e-11);

   // compare two-loop tadpoles
   double two_loop_tadpole[3];
   m.tadpole_hh_2loop(two_loop_tadpole);

   softsusy::numRewsbLoops = 2;
   s.doTadpoles(mt, sinthDRbar);

   const double td_1_and_2loop_ss = s.displayTadpole1Ms();
   const double tu_1_and_2loop_ss = s.displayTadpole2Ms();
   const double ts_1_and_2loop_ss = s.displayTadpoleSMs();

   BOOST_CHECK_CLOSE(two_loop_tadpole[0] / vd, td_1_and_2loop_ss - tadpole_ss_1, 1.0e-10);
   BOOST_CHECK_CLOSE(two_loop_tadpole[1] / vu, tu_1_and_2loop_ss - tadpole_ss_2, 3.0e-11);
   BOOST_CHECK_CLOSE(two_loop_tadpole[2] / vS, ts_1_and_2loop_ss - tadpole_ss_3, 1.0e-11);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_two_loop_ewsb )
{
   NMSSM_input_parameters input;
   input.m0 = 300.; // avoids tree-level tachyons
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;
   NMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   const double precision = 1.0e-5;
   setup_NMSSM_const(m, s, input);

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

   s.calcDrBarPars();
   m.calculate_DRbar_masses();

   const double mt = s.displayDrBarPars().mt;
   const int signMu = 1;

   m.set_ewsb_iteration_precision(precision);
   m.set_ewsb_loop_order(2);
   m.solve_ewsb();

   const Complex tadpole_hh_1(m.tadpole_hh(0));
   const Complex tadpole_hh_2(m.tadpole_hh(1));
   const Complex tadpole_hh_3(m.tadpole_hh(2));

   double two_loop_tadpole[3];
   m.tadpole_hh_2loop(two_loop_tadpole);

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);

   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1() - Re(tadpole_hh_1) - two_loop_tadpole[0], precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2() - Re(tadpole_hh_2) - two_loop_tadpole[1], precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3() - Re(tadpole_hh_3) - two_loop_tadpole[2], precision);

   softsusy::numRewsbLoops = 2;
   s.rewsb(signMu, mt);

   const double kappa_ss = s.displayKappa();
   const double vS_ss    = s.displaySvev();
   const double ms2_ss   = s.displayMsSquared();

   const double kappa_fs = m.get_Kappa();
   const double vS_fs    = m.get_vS();
   const double ms2_fs   = m.get_ms2();

   BOOST_CHECK_CLOSE_FRACTION(kappa_ss, kappa_fs, 1.0e-11);
   BOOST_CHECK_CLOSE_FRACTION(vS_ss   , vS_fs   , 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(ms2_ss  , ms2_fs  , 1.0e-10);
}
