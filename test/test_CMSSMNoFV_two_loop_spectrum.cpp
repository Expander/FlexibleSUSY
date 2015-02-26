
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMNoFV_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "test_CMSSMNoFV.hpp"
#include "wrappers.hpp"
#include "pv.hpp"
#include "CMSSMNoFV_two_scale_model.hpp"
#include "lowe.h"

using namespace flexiblesusy;
using namespace softsusy;
using namespace passarino_veltman;


BOOST_AUTO_TEST_CASE( test_CMSSMNoFV_two_loop_top_pole_mass )
{
   const QedQcd oneset;
   CMSSMNoFV_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.SignMu = 1;
   input.Azero = 0.;
   CMSSMNoFV<Two_scale> m;
   m.do_calculate_sm_pole_masses(true);
   setup_CMSSM_const(m, input);

   const double mt_pole_input = oneset.displayPoleMt();
   const double vu = m.get_vu();

   BOOST_MESSAGE("mt_pole(input) = " << mt_pole_input);

   // calculate DR-bar masses
   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.solve_ewsb();

   BOOST_MESSAGE("mt_drbar(guess) = " << m.get_MFt());

   unsigned iterations = 100;

   // calculate top DR-bar mass from top pole mass using two-loop
   // corrections
   do {
      Eigen::Matrix<double,3,3> mt_drbar_2loop(Eigen::Matrix<double,3,3>::Zero());
      mt_drbar_2loop(2,2) = m.calculate_MFt_DRbar(mt_pole_input, 2);
      m.set_Yu(((1.4142135623730951*mt_drbar_2loop)/vu).transpose());
      m.calculate_DRbar_masses();
      m.solve_ewsb();
   } while (--iterations);

   BOOST_MESSAGE("mt_drbar(2-loop) = " << m.get_MFt());

   m.set_pole_mass_loop_order(2);
   m.calculate_MFt_pole();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double mt_pole_2loop  = m.get_physical().MFt;

   BOOST_MESSAGE("mt_pole(2-loop) = " << mt_pole_2loop);

   BOOST_CHECK_CLOSE_FRACTION(mt_pole_input, mt_pole_2loop, 1.5e-3);
}
