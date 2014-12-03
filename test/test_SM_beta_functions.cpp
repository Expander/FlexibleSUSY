
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "SM_two_scale_model.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_beta_functions )
{
   SM_input_parameters input;
   SM<Two_scale> m;

   const double gut_norm = Sqrt(0.6); // GUT normalization

   input.LambdaIN = 0.1;

   const double vev = 246.;
   const double g1 = 0.2;
   const double g2 = 0.3;
   const double g3 = 0.4;

   m.set_scale(91.);
   m.set_v(246.);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.set_loops(1);

   const SM_soft_parameters beta(m.calc_beta());

   const double beta_g1 = oneOver16PiSqr * Power(m.get_g1(),3) * Sqr(gut_norm) * (41./6.);
   const double beta_g2 = oneOver16PiSqr * Power(m.get_g2(),3) * (-19./6.);
   const double beta_g3 = oneOver16PiSqr * Power(m.get_g3(),3) * (-7.);

   BOOST_CHECK_CLOSE_FRACTION(beta.get_g1(), beta_g1, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(beta.get_g2(), beta_g2, 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(beta.get_g3(), beta_g3, 1.0e-10);
}
