
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_initial_guesser

#include <boost/test/unit_test.hpp>
#include "test_MSSM.hpp"

#include "MSSM_model.hpp"
#include "MSSM_initial_guesser.hpp"
#include "mssm_parameter_point.hpp"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"

BOOST_AUTO_TEST_CASE( test_initial_guess )
{
   MSSM m;
   MSSM_input_parameters input;

   MSSM_low_scale_constraint  low_constraint(input);
   MSSM_susy_scale_constraint susy_constraint(input);
   MSSM_high_scale_constraint hig_constraint(input);

   MSSM_initial_guesser guesser(&m, input, low_constraint,
                                susy_constraint, hig_constraint);

   guesser.guess();

   Mssm_parameter_point pp;
   pp.m0 = input.m0;
   pp.m12 = input.m12;
   pp.a0 = input.Azero;
   pp.mxGuess = low_constraint.get_scale();
   pp.signMu = input.SignMu;
   pp.tanBeta = input.TanBeta;
   QedQcd oneset;
   Mssm<Two_scale> smssm;
   Mssm_initial_guesser initial_guesser(&smssm, pp.mxGuess, pp.tanBeta, pp.signMu, pp.get_soft_pars(), false);
   initial_guesser.set_QedQcd(oneset);

   initial_guesser.guess();

   test_parameter_equality(smssm, m);
}
