
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_beta_functions

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
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vd(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vu(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vS(), precision);

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
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vd(), 2.0e-09);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vu(), 1.0e-09);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_vS(), 1.0e-09);
}

BOOST_AUTO_TEST_CASE( test_NMSSM_one_loop_tadpoles )
{
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   s.calcDrBarPars();
   m.calculate_DRbar_parameters();

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

   BOOST_CHECK_CLOSE(Re(tadpole_hh_1) / vd, tadpole_ss_1, 1.0e-11);
   BOOST_CHECK_CLOSE(Re(tadpole_hh_2) / vu, tadpole_ss_2, 1.0e-11);
   BOOST_CHECK_CLOSE(Re(tadpole_hh_3) / vS, tadpole_ss_3, 1.0e-11);
}
