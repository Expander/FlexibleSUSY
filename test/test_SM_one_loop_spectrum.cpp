
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_one_loop_spectrum

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

   // check Higgs tree-level mass
   const double Mhh(m.get_Mhh());
   const double lambda = m.get_Lambdax();
   const double v = m.get_v();
   const double hh_tree = Sqrt(lambda * Sqr(v));

   BOOST_CHECK_CLOSE(Mhh, hh_tree, 1.0e-12);

   m.set_pole_mass_loop_order(1);
   m.do_calculate_sm_pole_masses(true);
   m.solve_ewsb_one_loop();
   m.calculate_pole_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   // check that mu2 was changed
   BOOST_CHECK(m.get_mu2() != mu2);

   // check Higgs pole mass
   const double Mhh_1l(m.get_physical().Mhh);
   const double hh_1l = Sqrt(-m.get_mu2() + 1.5*lambda*Sqr(v) - Re(m.self_energy_hh(Mhh_1l)));

   BOOST_CHECK_CLOSE(Mhh_1l, hh_1l, 2.0e-4);

   // check that tree-level Higgs mass has not changed
   BOOST_CHECK_CLOSE(m.get_Mhh(), hh_tree, 1.0e-12);
}
