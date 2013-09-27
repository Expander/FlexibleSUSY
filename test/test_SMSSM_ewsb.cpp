
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_beta_functions

#include <boost/test/unit_test.hpp>

#define private public

#include "test_SMSSM.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"
#include "logger.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void test_tree_level_ewsb(const SMSSM_input_parameters& input)
{
   SMSSM<Two_scale> m(input);
   NmssmSoftsusy s;
   const double precision = 1.0e-5;
   setup_SMSSM(m, s, input);

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
   s.rewsbTreeLevel(input.SignMu);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Mu() , s.displaySusyMu()   , precision * 40.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_BMu(), s.displayM3Squared(), precision * 4.);
   BOOST_CHECK_CLOSE_FRACTION(m.get_LL1(), s.displayXiS()      , precision);
}

BOOST_AUTO_TEST_CASE( test_SMSSM_ewsb_tree_level )
{
   SMSSM_input_parameters input;

   input.SignMu = 1;
   test_tree_level_ewsb(input);

   input.SignMu = -1;
   test_tree_level_ewsb(input);
}

BOOST_AUTO_TEST_CASE( test_SMSSM_ewsb_tree_level_via_soft_higgs_masses )
{
   SMSSM_input_parameters input;
   SMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_SMSSM(m, s, input);

   const int error = m.solve_ewsb_tree_level_via_soft_higgs_masses();

   BOOST_CHECK_EQUAL(error, 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vd(), 2.0e-09);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vu(), 3.0e-09);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vS(), 2.0e-09);
}

BOOST_AUTO_TEST_CASE( test_SMSSM_one_loop_tadpoles )
{
   // set-up non-tachyonic point
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
   softsusy::numRewsbLoops = 1;
   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

   // chech that there this is a valid point
   if (m.get_problems().have_problem()) {
      INFO("FlexibleSUSY problem detected: ");
      m.get_problems().print();
      INFO("");
   }
   BOOST_REQUIRE(!m.get_problems().have_problem());
   BOOST_REQUIRE(!s.displayProblem().test());

   const double mt = s.displayDrBarPars().mt;
   const double sinthDRbar = s.calcSinthdrbar();
   const double vd = m.get_vd();
   const double vu = m.get_vu();
   const double vS = m.get_vS();

   const Complex tadpole_hh_1(m.tadpole_hh(1));
   const Complex tadpole_hh_2(m.tadpole_hh(2));
   const Complex tadpole_hh_3(m.tadpole_hh(3));

   const double tadpole_ss_1 = s.doCalcTadpole1oneLoop(mt, sinthDRbar);
   const double tadpole_ss_2 = s.doCalcTadpole2oneLoop(mt, sinthDRbar);
   const double tadpole_ss_3 = s.doCalcTadpoleSoneLoop(mt, sinthDRbar);

   BOOST_CHECK_SMALL(Im(tadpole_hh_1), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_2), 1.0e-12);
   BOOST_CHECK_SMALL(Im(tadpole_hh_3), 1.0e-12);

   // TODO: increase the test precision here
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole_hh_1) / vd, tadpole_ss_1, 0.0003);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole_hh_2) / vu, tadpole_ss_2, 0.00007);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole_hh_3) / vS, tadpole_ss_3, 0.095);
}
