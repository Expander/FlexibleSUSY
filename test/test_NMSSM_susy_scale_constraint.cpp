
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_susy_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "test_NMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_ewsb_solver.hpp"
#include "NMSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_susy_scale_constraint.hpp"

#include <cmath>

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_susy_scale_constraint )
{
   QedQcd qedqcd;
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   setup_NMSSM(m, s, input);

   NMSSM_ewsb_solver<Two_scale> ewsb_solver;
   m.set_ewsb_solver(
      std::make_shared<NMSSM_ewsb_solver<Two_scale> >(ewsb_solver));

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

   const double precision = m.get_ewsb_iteration_precision();

   NMSSM_susy_scale_constraint<Two_scale> constraint(&m, qedqcd);
   constraint.apply();

   double tadpole[3];

   tadpole[0] = m.get_ewsb_eq_hh_1() - Re(m.tadpole_hh(0));
   tadpole[1] = m.get_ewsb_eq_hh_2() - Re(m.tadpole_hh(1));
   tadpole[2] = m.get_ewsb_eq_hh_3() - Re(m.tadpole_hh(2));

   if (m.get_ewsb_loop_order() > 1) {
      const auto two_loop_tadpole(m.tadpole_hh_2loop());
      tadpole[0] -= two_loop_tadpole[0];
      tadpole[1] -= two_loop_tadpole[1];
      tadpole[2] -= two_loop_tadpole[2];
   }

   // check that tree-level EWSB eqs. are fulfilled
   BOOST_CHECK_LT(std::fabs(tadpole[0]), 0.66);
   BOOST_CHECK_LT(std::fabs(tadpole[1]), 0.3);
   BOOST_CHECK_LT(std::fabs(tadpole[2]), 23.0);
}
