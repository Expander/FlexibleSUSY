#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMSemiAnalytic_consistent_solutions

#include <boost/test/unit_test.hpp>

#define private public

#include "test_CMSSMSemiAnalytic.hpp"
#include "CMSSMSemiAnalytic_semi_analytic_spectrum_generator.hpp"

#include "CMSSM_input_parameters.hpp"
#include "CMSSM_slha_io.hpp"
#include "CMSSM_two_scale_ewsb_solver.hpp"
#include "CMSSM_two_scale_spectrum_generator.hpp"

using namespace flexiblesusy;

CMSSM_input_parameters initialize_two_scale_input(
   double m0Sq, const CMSSMSemiAnalytic_input_parameters& semi_analytic_input)
{
   CMSSM_input_parameters input;

   input.m0 = Sqrt(m0Sq);
   input.m12 = semi_analytic_input.m12;
   input.Azero = semi_analytic_input.Azero;
   input.TanBeta = semi_analytic_input.TanBeta;
   input.SignMu = Sign(semi_analytic_input.MuInput);

   return input;
}

CMSSM<Two_scale> initialize_two_scale_model(
   const CMSSMSemiAnalytic<Semi_analytic>& semi_analytic_model,
   const CMSSMSemiAnalytic_input_parameters& semi_analytic_input)
{
   CMSSM<Two_scale> two_scale_model;

   two_scale_model.set_loops(semi_analytic_model.get_loops());
   two_scale_model.set_scale(semi_analytic_model.get_scale());
   two_scale_model.set_ewsb_loop_order(semi_analytic_model.get_ewsb_loop_order());
   two_scale_model.set_thresholds(semi_analytic_model.get_thresholds());

   two_scale_model.set_input_parameters(
      initialize_two_scale_input(semi_analytic_model.get_m0Sq(),
                                 semi_analytic_input));

   two_scale_model.set_Yu(semi_analytic_model.get_Yu());
   two_scale_model.set_Yd(semi_analytic_model.get_Yd());
   two_scale_model.set_Ye(semi_analytic_model.get_Ye());

   two_scale_model.set_g1(semi_analytic_model.get_g1());
   two_scale_model.set_g2(semi_analytic_model.get_g2());
   two_scale_model.set_g3(semi_analytic_model.get_g3());

   two_scale_model.set_Mu(semi_analytic_model.get_Mu());

   two_scale_model.set_vd(semi_analytic_model.get_vd());
   two_scale_model.set_vu(semi_analytic_model.get_vu());

   two_scale_model.set_TYu(semi_analytic_model.get_TYu());
   two_scale_model.set_TYd(semi_analytic_model.get_TYd());
   two_scale_model.set_TYe(semi_analytic_model.get_TYe());

   two_scale_model.set_mHd2(semi_analytic_model.get_mHd2());
   two_scale_model.set_mHu2(semi_analytic_model.get_mHu2());
   two_scale_model.set_BMu(semi_analytic_model.get_BMu());

   two_scale_model.set_mq2(semi_analytic_model.get_mq2());
   two_scale_model.set_mu2(semi_analytic_model.get_mu2());
   two_scale_model.set_md2(semi_analytic_model.get_md2());
   two_scale_model.set_ml2(semi_analytic_model.get_ml2());
   two_scale_model.set_me2(semi_analytic_model.get_me2());

   two_scale_model.set_MassB(semi_analytic_model.get_MassB());
   two_scale_model.set_MassWB(semi_analytic_model.get_MassWB());
   two_scale_model.set_MassG(semi_analytic_model.get_MassG());

   two_scale_model.calculate_DRbar_masses();

   return two_scale_model;
}

CMSSM<Two_scale> run_single_two_scale_iteration(const CMSSM<Two_scale>& model, const CMSSM_scales& scales)
{
   CMSSM<Two_scale> next_model(model);

   CMSSM_ewsb_solver<Two_scale> ewsb_solver;
   next_model.set_ewsb_solver(&ewsb_solver);

   CMSSM_high_scale_constraint<Two_scale> high_scale_constraint;
   CMSSM_susy_scale_constraint<Two_scale> susy_scale_constraint;
   CMSSM_low_scale_constraint<Two_scale> low_scale_constraint;

   high_scale_constraint.set_model(&next_model);
   susy_scale_constraint.set_model(&next_model);
   low_scale_constraint.set_model(&next_model);

   softsusy::QedQcd qedqcd;
   qedqcd.to(qedqcd.displayPoleMZ());

   susy_scale_constraint.set_sm_parameters(qedqcd);
   low_scale_constraint.set_sm_parameters(qedqcd);

   high_scale_constraint.initialize();
   susy_scale_constraint.initialize();
   low_scale_constraint.initialize();

   high_scale_constraint.scale = scales.HighScale;
   susy_scale_constraint.scale = scales.SUSYScale;
   low_scale_constraint.scale = scales.LowScale;

   // apply constraints once
   next_model.run_to(scales.LowScale);
   low_scale_constraint.apply();

   next_model.run_to(scales.HighScale);
   high_scale_constraint.apply();

   next_model.run_to(scales.SUSYScale);
   susy_scale_constraint.apply();

   return next_model;
}

double max_mass_rel_diff(const CMSSM<Two_scale>& old_model, const CMSSM<Two_scale>& new_model)
{
   std::array<double, 34> diff{};

   diff[0] = MaxRelDiff(old_model.get_MGlu(), new_model.get_MGlu());
   for (int i = 0; i < 6; i++) {
      diff[i + 1] = MaxRelDiff(old_model.get_MSd()(i), new_model.get_MSd()(i));
   }
   for (int i = 0; i < 3; i++) {
      diff[i + 7] = MaxRelDiff(old_model.get_MSv()(i), new_model.get_MSv()(i));
   }
   for (int i = 0; i < 6; i++) {
      diff[i + 10] = MaxRelDiff(old_model.get_MSu()(i), new_model.get_MSu()(i));
   }
   for (int i = 0; i < 6; i++) {
      diff[i + 16] = MaxRelDiff(old_model.get_MSe()(i), new_model.get_MSe()(i));
   }
   for (int i = 0; i < 2; i++) {
      diff[i + 22] = MaxRelDiff(old_model.get_Mhh()(i), new_model.get_Mhh()(i));
   }
   for (int i = 1; i < 2; i++) {
      diff[i + 24] = MaxRelDiff(old_model.get_MAh()(i), new_model.get_MAh()(i));
   }
   for (int i = 1; i < 2; i++) {
      diff[i + 26] = MaxRelDiff(old_model.get_MHpm()(i), new_model.get_MHpm()(i));
   }
   for (int i = 0; i < 4; i++) {
      diff[i + 28] = MaxRelDiff(old_model.get_MChi()(i), new_model.get_MChi()(i));
   }
   for (int i = 0; i < 2; i++) {
      diff[i + 32] = MaxRelDiff(old_model.get_MCha()(i), new_model.get_MCha()(i));
   }

   return *std::max_element(diff.cbegin(), diff.cend());
}

BOOST_AUTO_TEST_CASE( test_semi_analytic_to_two_scale )
{
   const double precision = 1.0e-4;

   softsusy::QedQcd qedqcd;

   CMSSMSemiAnalytic_input_parameters semi_analytic_input;
   semi_analytic_input.m12 = 660.;
   semi_analytic_input.TanBeta = 40.;
   semi_analytic_input.Azero = 0.;
   semi_analytic_input.MuInput = 497.;

   CMSSMSemiAnalytic_spectrum_generator<Semi_analytic> spectrum_generator;
   spectrum_generator.set_precision_goal(precision);
   spectrum_generator.set_max_iterations(0);
   spectrum_generator.set_calculate_sm_masses(0);

   spectrum_generator.run(qedqcd, semi_analytic_input);

   const auto& semi_analytic_problems = spectrum_generator.get_problems();

   BOOST_CHECK_EQUAL(semi_analytic_problems.have_problem(), false);

   const double high_scale = spectrum_generator.get_high_scale();
   const double susy_scale = spectrum_generator.get_susy_scale();
   const double low_scale = spectrum_generator.get_low_scale();

   const CMSSMSemiAnalytic<Semi_analytic> semi_analytic_model = spectrum_generator.get_model();

   CMSSM_scales scales;
   scales.HighScale = high_scale;
   scales.SUSYScale = susy_scale;
   scales.LowScale = low_scale;

   const CMSSM<Two_scale> two_scale_model(
      initialize_two_scale_model(semi_analytic_model, semi_analytic_input));

   const CMSSM<Two_scale> single_iteration_model =
      run_single_two_scale_iteration(two_scale_model, scales);

   // check that the originally found solution is an (approximate)
   // fixed point of the two-scale iteration as well
   const double mass_rel_diff = max_mass_rel_diff(
      two_scale_model, single_iteration_model);

   BOOST_CHECK_LT(mass_rel_diff, precision);
}
