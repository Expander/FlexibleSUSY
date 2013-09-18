
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSM_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "test_NMSSM.hpp"
#include "NMSSM_two_scale_model.hpp"
#include "NMSSM_two_scale_low_scale_constraint.hpp"
#include "nmssmsoftsusy.h"
#include "wrappers.hpp"
#include "ew_input.hpp"

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   NMSSM<Two_scale> m;
   NmssmSoftsusy s;
   NMSSM_input_parameters input;
   input.m0 = 250.; // avoids tree-level tachyons
   QedQcd oneset;
   setup_NMSSM(m, s, input);
   s.setData(oneset);

   m.calculate_DRbar_parameters();
   s.calcDrBarPars();

   NMSSM_low_scale_constraint<Two_scale> constraint(input, oneset);
   constraint.set_model(&m);

   const double alpha_em = oneset.displayAlpha(ALPHA);
   const double alpha_s  = oneset.displayAlpha(ALPHAS);
   const double scale = m.get_scale();

   const double delta_alpha_em_fs = constraint.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_fs  = constraint.calculate_delta_alpha_s(alpha_s);

   const double delta_alpha_em_ss = 1.0 - alpha_em / s.qedSusythresh(alpha_em, scale);
   const double delta_alpha_s_ss  = 1.0 - alpha_s  / s.qcdSusythresh(alpha_s , scale);

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_fs, delta_alpha_em_ss, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_fs , delta_alpha_s_ss , 1.0e-12);
}
