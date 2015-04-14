#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCPV_ewsb

#include <boost/test/unit_test.hpp>

#define private public

#include "test_CMSSMCPV.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "CMSSMCPV_two_scale_model.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_number_of_ewsb_eqs )
{
   const int neq = CMSSMCPV<Two_scale>::number_of_ewsb_equations;

   BOOST_CHECK_EQUAL(neq, 4);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_ewsb_tree_level )
{
   CMSSMCPV_input_parameters input;
   input.m0 = 125;
   input.m12 = 200;
   input.TanBeta = 10;
   input.Azero = 0;
   input.etaInput = 0.1;
   input.PhaseMu = std::complex<double>(1,0);

   CMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_CMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_Mu(input.m0);
   m.set_BMu(input.m0);

   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_1()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_2()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_3()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_4()) > precision);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_4(), precision);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_ewsb_one_loop )
{
   CMSSMCPV_input_parameters input;
   input.m0 = 125;
   input.m12 = 200;
   input.TanBeta = 10;
   input.Azero = 0;
   input.etaInput = 0.1;
   input.PhaseMu = std::complex<double>(1,0);

   CMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_CMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_Mu(input.m0);
   m.set_BMu(input.m0);

   const int neq = CMSSMCPV<Two_scale>::number_of_ewsb_equations;
   double tadpole[neq] = { 0. };

   m.tadpole_equations(tadpole);

   BOOST_CHECK(Abs(tadpole[0]) > 1.);
   BOOST_CHECK(Abs(tadpole[1]) > 1.);
   BOOST_CHECK(Abs(tadpole[2]) > 1.);
   BOOST_CHECK(Abs(tadpole[3]) > 1.);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_one_loop();

   BOOST_REQUIRE(error == 0);

   m.tadpole_equations(tadpole);

   BOOST_CHECK_SMALL(Abs(tadpole[0]), 0.003);
   BOOST_CHECK_SMALL(Abs(tadpole[1]), 0.0005);
   BOOST_CHECK_SMALL(Abs(tadpole[2]), precision);
   BOOST_CHECK_SMALL(Abs(tadpole[3]), precision);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_tree_level_tadpoles )
{
   CMSSMCPV_input_parameters input;
   input.m0 = 125;
   input.m12 = 200;
   input.TanBeta = 10;
   input.Azero = 0;
   input.etaInput = 0.1;
   input.PhaseMu = std::complex<double>(1,0);

   CMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_CMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_Mu(input.m0);
   m.set_BMu(input.m0);

   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_1()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_2()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_3()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_4()) > precision);

   // check relation between tree-level tadpoles 3 and 4
   const double vu = m.get_vu(), vd = m.get_vd();
   BOOST_CHECK_CLOSE_FRACTION(m.get_ewsb_eq_hh_3(), (vu/vd) * m.get_ewsb_eq_hh_4(), 1e-10);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_tree_level_tadpoles_real_limit )
{
   CMSSMCPV_input_parameters input;
   input.m0 = 125;
   input.m12 = 200;
   input.TanBeta = 10;
   input.Azero = 0;
   input.etaInput = 0.;
   input.PhaseMu = std::complex<double>(1,0);

   CMSSMCPV<Two_scale> m;
   const double precision = 1.0e-10;
   setup_CMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_Mu(input.m0);
   m.set_BMu(input.m0);

   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_1()) > 1.);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_2()) > 1.);
   BOOST_CHECK_LT(Abs(m.get_ewsb_eq_hh_3()), precision);
   BOOST_CHECK_LT(Abs(m.get_ewsb_eq_hh_4()), precision);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_one_loop_tadpoles_real_limit )
{
   CMSSMCPV_input_parameters input;
   input.m0 = 125;
   input.m12 = 200;
   input.TanBeta = 10;
   input.Azero = 0;
   input.etaInput = 0.;
   input.PhaseMu = std::complex<double>(1,0);

   CMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_CMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_Mu(input.m0);
   m.set_BMu(input.m0);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);

   m.calculate_DRbar_masses();

   BOOST_CHECK_GT(Abs(m.tadpole_hh(0)), precision);
   BOOST_CHECK_GT(Abs(m.tadpole_hh(1)), precision);
   BOOST_CHECK_SMALL(Abs(m.tadpole_hh(2)), precision);
   BOOST_CHECK_SMALL(Abs(m.tadpole_hh(3)), precision);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_one_loop_tadpoles )
{
   CMSSMCPV_input_parameters input;
   input.m0 = 125;
   input.m12 = 200;
   input.TanBeta = 10;
   input.Azero = 0;
   input.etaInput = 0.1;
   input.PhaseMu = std::complex<double>(1,0);

   CMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_CMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_Mu(input.m0);
   m.set_BMu(input.m0);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);

   m.calculate_DRbar_masses();

   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(0))), precision);
   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(1))), precision);
   BOOST_CHECK_GT(Abs(Im(m.tadpole_hh(2))), precision);
   BOOST_CHECK_GT(Abs(Im(m.tadpole_hh(3))), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(0)), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(1)), precision);
   BOOST_CHECK_SMALL(Re(m.tadpole_hh(2)), precision);
   BOOST_CHECK_SMALL(Re(m.tadpole_hh(3)), precision);

   // check relation between one-loop tadpoles 3 and 4
   const double vu = m.get_vu(), vd = m.get_vd();
   BOOST_CHECK_CLOSE_FRACTION(Im(m.tadpole_hh(2)), (vu/vd) * Im(m.tadpole_hh(3)), 1e-10);

   const int neq = CMSSMCPV<Two_scale>::number_of_ewsb_equations;
   double tadpole[neq] = { 0. };
   m.tadpole_equations(tadpole);

   // check relation between sum of tree-level and one-loop tadpoles 3 and 4
   BOOST_CHECK_CLOSE_FRACTION(tadpole[2], (vu/vd) * tadpole[3], 1e-10);
}
