#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMSemiAnalytic_semi_analytic_solutions

#include <boost/test/unit_test.hpp>

#include "test.h"
#include "test_CMSSMSemiAnalytic.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_solutions.hpp"

using namespace flexiblesusy;

struct Boundary_values {
   double m12{};
   double Azero{};
   double m0Sq{};
   double BMu0{};
   double Mu{};
};

void setup_high_scale_CMSSMSemiAnalytic(
   CMSSMSemiAnalytic_mass_eigenstates& model, const Boundary_values& values)
{
   model.set_TYu(model.get_Yu() * values.Azero);
   model.set_TYd(model.get_Yd() * values.Azero);
   model.set_TYe(model.get_Ye() * values.Azero);

   model.set_mq2(values.m0Sq * UNITMATRIX(3));
   model.set_mu2(values.m0Sq * UNITMATRIX(3));
   model.set_md2(values.m0Sq * UNITMATRIX(3));
   model.set_ml2(values.m0Sq * UNITMATRIX(3));
   model.set_me2(values.m0Sq * UNITMATRIX(3));
   model.set_mHd2(values.m0Sq);
   model.set_mHu2(values.m0Sq);

   model.set_MassB(values.m12);
   model.set_MassWB(values.m12);
   model.set_MassG(values.m12);

   model.set_Mu(values.Mu);
   model.set_BMu(values.BMu0);
}

BOOST_AUTO_TEST_CASE( test_CMSSMSemiAnalytic_coefficients )
{
   CMSSMSemiAnalytic_input_parameters input;
   CMSSMSemiAnalytic_mass_eigenstates model(input);
   setup_CMSSMSemiAnalytic(model, input);

   const double high_scale = 2.e16;
   model.run_to(high_scale);

   Boundary_values values;
   values.m12 = 300.0;
   values.Azero = -200.;
   values.m0Sq = Sqr(250);
   values.BMu0 = -1.e4;
   values.Mu = model.get_Mu();

   setup_high_scale_CMSSMSemiAnalytic(model, values);

   CMSSMSemiAnalytic_semi_analytic_solutions solns;
   solns.set_input_scale(high_scale);
   solns.set_output_scale(Electroweak_constants::MZ);
   solns.set_AzeroBasis(values.Azero);
   solns.set_m12Basis(values.m12);
   solns.set_m0SqBasis(values.m0Sq);
   solns.set_BMu0Basis(values.BMu0);
   solns.set_MuBasis(values.Mu);

   solns.calculate_coefficients(model);

   model.run_to(Electroweak_constants::MZ);

   CMSSMSemiAnalytic_mass_eigenstates coeffs_model(model);
   solns.evaluate_solutions(coeffs_model);

   BOOST_CHECK_CLOSE_FRACTION(model.get_MassB(), coeffs_model.get_MassB(), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_MassWB(), coeffs_model.get_MassWB(), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_MassG(), coeffs_model.get_MassG(), 1.0e-3);

   BOOST_CHECK_CLOSE_FRACTION(model.get_BMu(), coeffs_model.get_BMu(), 1.0e-3);

   BOOST_CHECK_CLOSE_FRACTION(model.get_mHd2(), coeffs_model.get_mHd2(), 1.0e-3);
   BOOST_CHECK_CLOSE_FRACTION(model.get_mHu2(), coeffs_model.get_mHu2(), 1.0e-3);

   TEST_CLOSE_REL(model.get_TYu(), coeffs_model.get_TYu(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYd(), coeffs_model.get_TYd(), 1.0e-3);
   TEST_CLOSE_REL(model.get_TYe(), coeffs_model.get_TYe(), 1.0e-3);

   TEST_CLOSE_REL(model.get_mq2(), coeffs_model.get_mq2(), 1.0e-2);
   TEST_CLOSE_REL(model.get_mu2(), coeffs_model.get_mu2(), 1.0e-2);
   TEST_CLOSE_REL(model.get_md2(), coeffs_model.get_md2(), 1.0e-2);
   TEST_CLOSE_REL(model.get_ml2(), coeffs_model.get_ml2(), 1.0e-2);
   TEST_CLOSE_REL(model.get_me2(), coeffs_model.get_me2(), 1.0e-2);

   BOOST_CHECK_EQUAL(get_errors(), 0);
}
