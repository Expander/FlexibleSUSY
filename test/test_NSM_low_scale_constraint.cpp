
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NSM_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "NSM_two_scale_model.hpp"
#include "NSM_two_scale_low_scale_constraint.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_model_properties )
{
   BOOST_CHECK(NSM_info::is_supersymmetric_model == false);
}

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   NSM<Two_scale> m;
   NSM_input_parameters input;
   input.LambdaInput1 = 0.1;
   input.LambdaInput2 = 0.1;
   input.LambdaInput3 = 0.1;
   input.LambdaInput4 = 0.0;
   input.LambdaInput5 = 0.0;
   QedQcd oneset;

   const double vev = 246.;

   m.set_scale(91.);
   m.set_vH(246.);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.calculate_DRbar_masses();

   NSM_low_scale_constraint<Two_scale> constraint(&m, input, oneset);

   const double alpha_em = oneset.displayAlpha(ALPHA);
   const double alpha_s  = oneset.displayAlpha(ALPHAS);
   const double scale = m.get_scale();

   const double delta_alpha_em_fs = constraint.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_fs  = constraint.calculate_delta_alpha_s(alpha_s);

   const double Mtop = m.get_MFu(2);

   // the extra singlet field does not couple electromagnetically and
   // does thus not contribute here (or below)
   const double delta_alpha_em =
      alpha_em / (2 * Pi) * (1./3. - 16./9. * log(Mtop/scale));

   // no MS-bar DR-bar conversion term appears here, because the NSM
   // is renormalized in MS-bar
   const double delta_alpha_s  =
      alpha_s / (2 * Pi) * (- 2./3. * log(Mtop/scale));

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_fs, delta_alpha_em, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_fs , delta_alpha_s , 1.0e-12);
}
