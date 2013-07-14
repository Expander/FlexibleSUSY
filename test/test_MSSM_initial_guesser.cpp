
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_initial_guesser

#include <boost/test/unit_test.hpp>
#include "test_MSSM.hpp"

#include "MSSM_model.hpp"
#include "MSSM_initial_guesser.hpp"
#include "mssm_parameter_point.hpp"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale_sugra_constraint.hpp"
#include "mssm_two_scale_msusy_constraint.hpp"
#include "mssm_two_scale_mz_constraint.hpp"

BOOST_AUTO_TEST_CASE( test_initial_guess )
{
   MSSM m;
   Mssm<Two_scale> smssm;

   // create MSSM initial guesser
   MSSM_input_parameters input;

   MSSM_low_scale_constraint  low_constraint(input);
   MSSM_susy_scale_constraint susy_constraint(input);
   MSSM_high_scale_constraint high_constraint(input);

   MSSM_initial_guesser guesser(&m, input, low_constraint,
                                susy_constraint, high_constraint);

   // create Mssm initial guesser
   Mssm_parameter_point pp;
   pp.m0 = input.m0;
   pp.m12 = input.m12;
   pp.a0 = input.Azero;
   pp.mxGuess = high_constraint.get_scale();
   pp.signMu = input.SignMu;
   pp.tanBeta = input.TanBeta;
   QedQcd oneset;
   Mssm_sugra_constraint mssm_sugra_constraint(pp);
   Mssm_mz_constraint mssm_mz_constraint(pp);
   Mssm_msusy_constraint mssm_msusy_constraint(pp);
   Mssm_initial_guesser initial_guesser(&smssm, pp, mssm_mz_constraint,
                                        mssm_msusy_constraint,
                                        mssm_sugra_constraint);
   initial_guesser.set_QedQcd(oneset);

   // guess both models
   guesser.guess();
   initial_guesser.guess();

   BOOST_CHECK_EQUAL(smssm.displayLoops()     , m.get_loops());
   BOOST_CHECK_EQUAL(smssm.displayMu()        , m.get_scale());
   BOOST_CHECK_EQUAL(smssm.displayThresholds(), m.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(smssm.displayGaugeCoupling(1), m.get_g1(), 1.0e-6);
   BOOST_CHECK_CLOSE_FRACTION(smssm.displayGaugeCoupling(2), m.get_g2(), 1.0e-5);
   BOOST_CHECK_CLOSE_FRACTION(smssm.displayGaugeCoupling(3), m.get_g3(), 1.0e-5);

   // test off-diagonal elements
   BOOST_MESSAGE("testing off-diagonal yukawa elements");
   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         if (i == k)
            continue;
         BOOST_MESSAGE("testing yukawa elements " << i << ", " << k);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(i-1,k-1), smssm.displayYukawaMatrix(YU)(i,k), 0.00001);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(i-1,k-1), smssm.displayYukawaMatrix(YD)(i,k), 0.00001);
         BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(i-1,k-1), smssm.displayYukawaMatrix(YE)(i,k), 0.00001);
      }
   }

   BOOST_MESSAGE("testing diagonal yukawa elements");
   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(0,0), smssm.displayYukawaMatrix(YU)(1,1), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(0,0), smssm.displayYukawaMatrix(YD)(1,1), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(0,0), smssm.displayYukawaMatrix(YE)(1,1), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(1,1), smssm.displayYukawaMatrix(YU)(2,2), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(1,1), smssm.displayYukawaMatrix(YD)(2,2), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(1,1), smssm.displayYukawaMatrix(YE)(2,2), 1.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(m.get_Yu()(2,2), smssm.displayYukawaMatrix(YU)(3,3), 1.0e-4);
   // BOOST_CHECK_CLOSE_FRACTION(m.get_Yd()(2,2), smssm.displayYukawaMatrix(YD)(3,3), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(m.get_Ye()(2,2), smssm.displayYukawaMatrix(YE)(3,3), 1.0e-5);

   const double vu = m.get_vu(), vd = m.get_vd();
   double tanBeta;
   if (is_zero(vu))
      tanBeta = 0.;
   else if (is_zero(vd))
      tanBeta = std::numeric_limits<double>::infinity();
   else
      tanBeta = vu/vd;
   const double vev = Sqrt(Sqr(m.get_vu()) + Sqr(m.get_vd()));
   BOOST_CHECK_CLOSE_FRACTION(smssm.displayTanb(), tanBeta, 3.0e-4);
   BOOST_CHECK_CLOSE_FRACTION(smssm.displayHvev(), vev    , 1.0e-5);

   // BOOST_CHECK_CLOSE_FRACTION(smssm.displayGaugino(1), m.get_MassB() , 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(smssm.displayGaugino(2), m.get_MassWB(), 1.0e-5);
   // BOOST_CHECK_CLOSE_FRACTION(smssm.displayGaugino(3), m.get_MassG() , 1.0e-5);
}
