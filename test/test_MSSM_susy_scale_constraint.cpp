
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
   MSSM<Two_scale> m; MssmSoftsusy s;
   MSSM_input_parameters input;
   setup_MSSM(m, s, input);

   MSSM_susy_scale_constraint<Two_scale> constraint(input);
   constraint.set_model(&m);
   constraint.apply();

   // check that tree-level EWSB eqs. are fulfilled
   BOOST_CHECK_SMALL(std::fabs(m.get_ewsb_eq_vd() - m.tadpole_hh(0)), 1.2e-8);
   BOOST_CHECK_SMALL(std::fabs(m.get_ewsb_eq_vu() - m.tadpole_hh(1)), 1.2e-7);
}
