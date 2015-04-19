#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_NMSSMCPV_ewsb

#include <boost/test/unit_test.hpp>

#define private public

#include "test_NMSSMCPV.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "NMSSMCPV_two_scale_model.hpp"

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_number_of_ewsb_eqs )
{
   const int neq = NMSSMCPV<Two_scale>::number_of_ewsb_equations;

   BOOST_CHECK_EQUAL(neq, 6);
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_ewsb_tree_level )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 250.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaTimesvSInput = 100;
   input.etaInput = 0.1;
   input.etaSInput = 0.1;

   NMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_NMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_1()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_2()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_3()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_4()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_5()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_6()) > precision);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_4(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_5(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_6(), precision);
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_ewsb_one_loop )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 500.;
   input.m12 = 500.;
   input.TanBeta = 5.;
   input.Azero = -300.;
   input.LambdaInput = -0.2;
   input.KappaInput = 0.3;
   input.LambdaTimesvSInput = 1000;
   input.etaInput = 0.1;
   input.etaSInput = 0.1;

   NMSSMCPV<Two_scale> m;
   const double precision = 1.0e-10;
   setup_NMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_ms2(-Sqr(input.m0));

   const int neq = NMSSMCPV<Two_scale>::number_of_ewsb_equations;
   double tadpole[neq] = { 0. };

   m.tadpole_equations(tadpole);

   BOOST_CHECK(Abs(tadpole[0]) > 100.);
   BOOST_CHECK(Abs(tadpole[1]) > 100.);
   BOOST_CHECK(Abs(tadpole[2]) > 100.);
   BOOST_CHECK(Abs(tadpole[3]) > 100.);
   BOOST_CHECK(Abs(tadpole[4]) > 100.);
   BOOST_CHECK(Abs(tadpole[5]) > 100.);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_one_loop();

   BOOST_REQUIRE(error == 0);

   m.tadpole_equations(tadpole);

   // tadpole[0,1,2,4,5] are fixed by the three EWSB eqs.  tadpole[3]
   // is not directly fixed by any EWSB eq.  However, the relation
   //
   //    tadpole[3]/tadpole[4] = vu/vd           (1)
   //
   // is valid, up to higher order effects.  After the iteration
   // tadpole[3] is not exactly zero, because during the iteration the
   // relation (1) is spoiled by higher order effects, see the
   // CMSSMCPV EWSB test for example.

   BOOST_CHECK_SMALL(Abs(tadpole[0]), 0.003);     // fixed by EWSB
   BOOST_CHECK_SMALL(Abs(tadpole[1]), 0.0005);    // fixed by EWSB
   BOOST_CHECK_SMALL(Abs(tadpole[2]), 0.008);     // fixed by EWSB
   BOOST_CHECK_SMALL(Abs(tadpole[3]), 3.);
   BOOST_CHECK_SMALL(Abs(tadpole[4]), 0.002);     // fixed by EWSB
   BOOST_CHECK_SMALL(Abs(tadpole[5]), 0.002);     // fixed by EWSB
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_tree_level_tadpoles )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 250.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaTimesvSInput = 100;
   input.etaInput = 0.1;
   input.etaSInput = 0.1;

   NMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_NMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_1()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_2()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_3()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_4()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_5()) > precision);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_6()) > precision);

   // check relation between tree-level tadpoles 3 and 4
   const double vu = m.get_vu(), vd = m.get_vd();
   BOOST_CHECK_CLOSE_FRACTION(m.get_ewsb_eq_hh_4(), (vu/vd) * m.get_ewsb_eq_hh_5(), 1e-10);
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_tree_level_tadpoles_real_limit )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 250.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaTimesvSInput = 100;
   input.etaInput = 0;
   input.etaSInput = 0;

   NMSSMCPV<Two_scale> m;
   const double precision = 1.0e-10;
   setup_NMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_1()) > 1.);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_2()) > 1.);
   BOOST_CHECK(Abs(m.get_ewsb_eq_hh_3()) > 1.);
   BOOST_CHECK_LT(Abs(m.get_ewsb_eq_hh_4()), precision);
   BOOST_CHECK_LT(Abs(m.get_ewsb_eq_hh_5()), precision);
   BOOST_CHECK_LT(Abs(m.get_ewsb_eq_hh_6()), precision);
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_one_loop_tadpoles_real_limit )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 250.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = -0.1;
   input.KappaInput = 0.1;
   input.LambdaTimesvSInput = 100;
   // eta and etaS must not be both equal to zero, because this will
   // lead to a division by zero in the tree-level EWSB eqs.
   input.etaInput = 1e-8;
   input.etaSInput = 1e-10;

   NMSSMCPV<Two_scale> m;
   setup_NMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));
   m.set_ms2(-Sqr(input.m0));

   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);

   m.calculate_DRbar_masses();

   BOOST_CHECK_GT(Abs(m.tadpole_hh(0)), 100);
   BOOST_CHECK_GT(Abs(m.tadpole_hh(1)), 100);
   BOOST_CHECK_GT(Abs(m.tadpole_hh(2)), 100);
   BOOST_CHECK_SMALL(Abs(m.tadpole_hh(3)), 0.2);
   BOOST_CHECK_SMALL(Abs(m.tadpole_hh(4)), 0.02);
   BOOST_CHECK_SMALL(Abs(m.tadpole_hh(5)), 0.1);
}

BOOST_AUTO_TEST_CASE( test_NMSSMCPV_one_loop_tadpoles )
{
   NMSSMCPV_input_parameters input;
   input.m0 = 250.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.LambdaTimesvSInput = 100;
   input.etaInput = 0.1;
   input.etaSInput = 0.1;

   NMSSMCPV<Two_scale> m;
   const double precision = 1.0e-5;
   setup_NMSSMCPV(m, input);

   // initial guess
   m.set_mHu2(-Sqr(input.m0));
   m.set_mHd2(Sqr(input.m0));

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);

   m.calculate_DRbar_masses();

   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(0))), precision);
   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(1))), precision);
   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(2))), precision);
   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(3))), precision);
   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(4))), precision);
   BOOST_CHECK_GT(Abs(Re(m.tadpole_hh(5))), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(0)), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(1)), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(2)), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(3)), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(4)), precision);
   BOOST_CHECK_SMALL(Im(m.tadpole_hh(5)), precision);

   const double vu = m.get_vu(), vd = m.get_vd();

   // check relation between one-loop tadpoles 3 and 4
   BOOST_CHECK_CLOSE_FRACTION(Re(m.tadpole_hh(3)), (vu/vd) * Re(m.tadpole_hh(4)), 1e-10);

   const int neq = NMSSMCPV<Two_scale>::number_of_ewsb_equations;
   double tadpole[neq] = { 0. };
   m.tadpole_equations(tadpole);

   // check relation between sum of tree-level and one-loop tadpoles 3 and 4
   BOOST_CHECK_CLOSE_FRACTION(tadpole[3], (vu/vd) * tadpole[4], 1e-10);
}
