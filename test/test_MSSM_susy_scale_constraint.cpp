
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_susy_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_MSSM.hpp"

#define private public

#include "MSSM_two_scale_model.hpp"
#include "MSSM_two_scale_susy_scale_constraint.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include <cmath>

BOOST_AUTO_TEST_CASE( test_susy_scale_constraint )
{
   MSSM_input_parameters input;
   MSSM<Two_scale> m; MssmSoftsusy s;
   setup_MSSM(m, s, input);

   MSSM_susy_scale_constraint<Two_scale> constraint(input);
   constraint.set_model(&m);
   constraint.apply();

   double tadpole[2];

   tadpole[0] = m.get_ewsb_eq_hh_1() - Re(m.tadpole_hh(0));
   tadpole[1] = m.get_ewsb_eq_hh_2() - Re(m.tadpole_hh(1));

   if (m.get_ewsb_loop_order() > 1) {
      double two_loop_tadpole[2];
      m.tadpole_hh_2loop(two_loop_tadpole);
      tadpole[0] -= two_loop_tadpole[0];
      tadpole[1] -= two_loop_tadpole[1];
   }

   // check that EWSB eqs. are fulfilled
   BOOST_CHECK_SMALL(std::fabs(tadpole[0]), 1.2e-7);
   BOOST_CHECK_SMALL(std::fabs(tadpole[1]), 1.5e-7);
}
