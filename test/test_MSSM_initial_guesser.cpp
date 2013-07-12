
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_initial_guesser

#include <boost/test/unit_test.hpp>
#include "test_MSSM.hpp"

#include "MSSM_model.hpp"
#include "MSSM_initial_guesser.hpp"

BOOST_AUTO_TEST_CASE( test_initial_guess )
{
   MSSM m; MssmSoftsusy s;
   MSSM_input_parameters input;

   MSSM_low_scale_constraint  low_constraint(input);
   MSSM_susy_scale_constraint susy_constraint(input);
   MSSM_high_scale_constraint hig_constraint(input);

   MSSM_initial_guesser guesser(&m, input, low_constraint,
                                susy_constraint, hig_constraint);

   guesser.guess();

   const double mxGuess = low_constraint.get_scale();
   DoubleVector pars(3);
   pars(1) = input.m0;
   pars(2) = input.m12;
   pars(3) = input.Azero;
   const int signMu = input.SignMu;
   const double tanBeta = input.TanBeta;
   const bool uni = true;
   QedQcd oneset;

   softsusy::TOLERANCE = 10.; // ensures that itLowsoft will return immediately
   s.lowOrg(sugraBcs, mxGuess, pars, signMu, tanBeta, oneset, uni);

   // test_parameter_equality(s, m);
}
