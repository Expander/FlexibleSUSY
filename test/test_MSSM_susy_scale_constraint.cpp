
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_susy_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test_MSSM.hpp"

#define private public

#include "MSSM_model.hpp"
#include "MSSM_susy_scale_constraint.hpp"
#include "softsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include <cmath>

BOOST_AUTO_TEST_CASE( test_susy_scale_constraint )
{
   MSSM m; MssmSoftsusy s;
   MSSM_input_parameters input;
   setup_MSSM(m, s, input);

   MSSM_susy_scale_constraint constraint(input);
   constraint.set_model(&m);
   constraint.apply();

   // check that tree-level EWSB eqs. are fulfilled
   BOOST_CHECK_SMALL(std::fabs(m.get_ewsb_eq_vd()), 4.0e-9);
   BOOST_CHECK_SMALL(std::fabs(m.get_ewsb_eq_vu()), 6.0e-10);
}
