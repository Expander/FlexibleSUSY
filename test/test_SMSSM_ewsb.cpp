
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SMSSM_beta_functions

#include <boost/test/unit_test.hpp>

#define private public

#include "test_SMSSM.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "SMSSM_two_scale_model.hpp"

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
