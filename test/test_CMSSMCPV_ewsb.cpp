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

   // check relation between tree-level tadpoles 3 and 4
   const double vu = m.get_vu(), vd = m.get_vd();
   BOOST_CHECK_CLOSE_FRACTION(m.get_ewsb_eq_hh_3(), -(vu/vd) * m.get_ewsb_eq_hh_4(), 1e-10);

   m.set_ewsb_iteration_precision(precision);
   const int error = m.solve_ewsb_tree_level();

   BOOST_REQUIRE(error == 0);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3(), precision);
   BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_4(), precision);
}
