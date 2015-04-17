
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_one_loop_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SM.hpp"
#include "wrappers.hpp"
#include "pv.hpp"
#include "SM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace passarino_veltman;

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

BOOST_AUTO_TEST_CASE( test_SM_wz_self_energies )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.calculate_DRbar_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double p = 100.;

   const double se_w_heavy = Re(m.self_energy_VWp_heavy(p));
   const double se_z_heavy = Re(m.self_energy_VZ_heavy(p));

   BOOST_CHECK_SMALL(se_w_heavy, 1.0e-10);
   BOOST_CHECK_SMALL(se_z_heavy, 1.0e-10);
}

BOOST_AUTO_TEST_CASE( test_SM_heavy_top_self_energy )
{
   SM_input_parameters input;
   input.LambdaIN = 0.25;
   SM<Two_scale> m;
   setup_SM_const(m, input);

   m.calculate_DRbar_masses();

   if (m.get_problems().have_problem()) {
      std::ostringstream ostr;
      m.get_problems().print_problems(ostr);
      BOOST_FAIL(ostr.str());
   }

   const double p = 100.;

   const Eigen::Array<double,3,1> MFu(m.get_MFu());
   const double g3 = m.get_g3();
   const double scale = m.get_scale();

   Eigen::Matrix<std::complex<double>,3,3> se_t;
   Eigen::Matrix<std::complex<double>,3,3> se_t_no_gluon;

   // top self-energies with and without gluon
   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         se_t(i,k) = m.self_energy_Fu_1(p,i,k);
         se_t_no_gluon(i,k) = m.self_energy_Fu_1_heavy_rotated(p,i,k);
      }
   }

   // adding gluon contrbution
   Eigen::Matrix<std::complex<double>,3,3> se_t_check(se_t_no_gluon);

   for (unsigned i = 0; i < 3; i++) {
      const double gluon_contrib =
         -5.333333333333333 *
         (-0.5 + ReB0(Sqr(p),Sqr(MFu(i)),0,Sqr(scale)))
            * Sqr(g3) * MFu(i);

      se_t_check(i,i) += gluon_contrib * oneOver16PiSqr;
   }

   const Eigen::Matrix<std::complex<double>,3,3> Uu(m.get_Uu());
   const Eigen::Matrix<std::complex<double>,3,3> Vu(m.get_Vu());

   se_t_check = Vu.transpose() * se_t_check * Uu;

   for (unsigned i = 0; i < 3; i++) {
      for (unsigned k = 0; k < 3; k++) {
         BOOST_CHECK_CLOSE(Re(se_t(i,k)), Re(se_t_check(i,k)), 1.0e-10);
         BOOST_CHECK_CLOSE(Im(se_t(i,k)), Im(se_t_check(i,k)), 1.0e-10);
      }
   }

}
