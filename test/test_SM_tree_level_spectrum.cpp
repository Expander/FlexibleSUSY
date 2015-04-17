
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "wrappers.hpp"
#include "SM_two_scale_model.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_SM_tree_level_masses )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   // set mu2 which is not in agreement with EWSB
   const double mu2 = 100.;
   m.set_mu2(mu2);

   m.calculate_DRbar_masses();

   // check that mu2 was not changed
   BOOST_CHECK_CLOSE(m.get_mu2(), mu2, 1.0e-10);

   // neutral CP even Higgs
   const double hh(m.get_Mhh());
   const double lambda = m.get_Lambdax();
   const double v = m.get_v();
   const double hh_tree = Sqrt(lambda * Sqr(v));
   BOOST_CHECK_CLOSE(hh, hh_tree, 1.0e-12);

   // check Goldstone boson masses
   BOOST_CHECK_EQUAL(m.get_MHp(), m.get_MVWp());
   BOOST_CHECK_EQUAL(m.get_MAh(), m.get_MVZ());
}
