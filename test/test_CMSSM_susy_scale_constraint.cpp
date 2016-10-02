
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSM_susy_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_CMSSM.hpp"

#define private public

#include "CMSSM_two_scale_model.hpp"
#include "CMSSM_two_scale_susy_scale_constraint.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include <cmath>

BOOST_AUTO_TEST_CASE( test_susy_scale_constraint )
{
   QedQcd qedqcd;
   CMSSM_input_parameters input;
   CMSSM<Two_scale> m; MssmSoftsusy s;
   setup_CMSSM(m, s, input);

   CMSSM_susy_scale_constraint<Two_scale> constraint(&m, qedqcd);
   constraint.apply();

   double tadpole[2];

   tadpole[0] = m.get_ewsb_eq_hh_1() - Re(m.tadpole_hh(0));
   tadpole[1] = m.get_ewsb_eq_hh_2() - Re(m.tadpole_hh(1));

   if (m.get_ewsb_loop_order() > 1) {
      const auto two_loop_tadpole(m.tadpole_hh_2loop());
      tadpole[0] -= two_loop_tadpole[0];
      tadpole[1] -= two_loop_tadpole[1];
   }

   // check that EWSB eqs. are fulfilled
   BOOST_CHECK_SMALL(std::fabs(tadpole[0]), 0.3);
   BOOST_CHECK_SMALL(std::fabs(tadpole[1]), 0.1);
}
