
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_spectrum

#include <boost/test/unit_test.hpp>

#define private public

#include "nmssmsoftsusy.h"
#include "error.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "NMSSM_two_scale_model.hpp"
#include "NMSSM_input_parameters.hpp"
#include "NMSSM_two_scale_high_scale_constraint.hpp"
#include "NMSSM_two_scale_susy_scale_constraint.hpp"
#include "NMSSM_two_scale_low_scale_constraint.hpp"
#include "NMSSM_two_scale_convergence_tester.hpp"
#include "NMSSM_two_scale_initial_guesser.hpp"
#include "test_NMSSM.hpp"

class SoftSusy_error : public Error {
public:
   SoftSusy_error(const std::string& msg_)
      : msg(msg_) {}
   virtual ~SoftSusy_error() {}
   virtual std::string what() const { return msg; }
private:
   std::string msg;
};

class SoftSusy_NoConvergence_error : public SoftSusy_error {
public:
   SoftSusy_NoConvergence_error(const std::string& msg_)
      : SoftSusy_error(msg_) {}
   virtual ~SoftSusy_NoConvergence_error() {}
   virtual std::string what() const { return SoftSusy_error::what(); }
};

class SoftSusy_NonPerturbative_error : public SoftSusy_error {
public:
   SoftSusy_NonPerturbative_error(const std::string& msg_)
      : SoftSusy_error(msg_) {}
   virtual ~SoftSusy_NonPerturbative_error() {}
   virtual std::string what() const { return SoftSusy_error::what(); }
};

class SoftSusy_tester {
public:
   SoftSusy_tester()
      : mx(0.0), msusy(0.0), softSusy(), gaugeUnification(true) {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   sPhysical get_physical() const { return softSusy.displayPhys(); }
   NmssmSoftsusy get_model() const { return softSusy; }
   void test(const NMSSM_input_parameters& pp, double mxGuess, const QedQcd& oneset = QedQcd()) {
      // run softsusy
      softsusy::numRewsbLoops = 1;
      softsusy::numHiggsMassLoops = 1;
      softsusy::TOLERANCE = 1.0e-4;
      softsusy::Z3 = true;
      softsusy::GUTlambda = true;
#ifdef VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      DoubleVector pars(3);
      pars(1) = pp.m0;
      pars(2) = pp.m12;
      pars(3) = pp.Azero;
      DoubleVector nmpars(5);
      nmpars(1) = pp.LambdaInput;
      nmpars(2) = 0.1;   // initial guess at the low scale
      nmpars(3) = 1000.; // initial guess at the low scale
      nmpars(4) = 0.;
      nmpars(5) = 0.;

      softSusy.setAlternativeMs(true);
      softSusy.lowOrg(NmssmMsugraBcs, mxGuess, pars, nmpars, 1, pp.TanBeta,
                      oneset, gaugeUnification);
      mx = softSusy.displayMxBC();
      msusy = softSusy.displayMsusy();
      softsusy::PRINTOUT = 0;

      if (softSusy.displayProblem().test()) {
         std::stringstream ss;
         ss << "SoftSusy problem: " << softSusy.displayProblem();
         VERBOSE_MSG(ss.str());
         if (softSusy.displayProblem().noConvergence)
            throw SoftSusy_NoConvergence_error(ss.str());
         else if (softSusy.displayProblem().nonperturbative)
            throw SoftSusy_NonPerturbative_error(ss.str());
         else
            throw SoftSusy_error(ss.str());
      }
   }
private:
   double mx, msusy;
   NmssmSoftsusy softSusy;
   bool gaugeUnification;
};

class NMSSM_tester {
public:
   NMSSM_tester()
      : mx(0.0), msusy(0.0), mssm()
      , high_constraint(NULL), susy_constraint(NULL), low_constraint(NULL) {}
   ~NMSSM_tester() {
      delete high_constraint;
      delete susy_constraint;
      delete low_constraint;
   }
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   NMSSM_physical get_physical() const { return mssm.get_physical(); }
   NMSSM<Two_scale> get_model() const { return mssm; }
   void set_low_scale_constraint(NMSSM_low_scale_constraint<Two_scale>* c) { low_constraint = c; }
   void set_susy_scale_constraint(NMSSM_susy_scale_constraint<Two_scale>* c) { susy_constraint = c; }
   void set_high_scale_constraint(NMSSM_high_scale_constraint<Two_scale>* c) { high_constraint = c; }
   void setup_default_constaints() {
      if (!high_constraint)
         high_constraint = new NMSSM_high_scale_constraint<Two_scale>();
      if (!susy_constraint)
         susy_constraint = new NMSSM_susy_scale_constraint<Two_scale>();
      if (!low_constraint)
         low_constraint = new NMSSM_low_scale_constraint<Two_scale>();
   }
   void test(const NMSSM_input_parameters& pp, const QedQcd& oneset = QedQcd()) {
      setup_default_constaints();
      high_constraint->set_input_parameters(pp);
      low_constraint->set_input_parameters(pp);
      low_constraint->set_sm_parameters(oneset);
      susy_constraint->set_input_parameters(pp);

      NMSSM_convergence_tester<Two_scale> convergence_tester(&mssm, 1.0e-4);
      NMSSM_initial_guesser<Two_scale> initial_guesser(&mssm, pp, oneset,
                                                      *low_constraint,
                                                      *susy_constraint,
                                                      *high_constraint);
      Two_scale_increasing_precision precision(10.0, 1.0e-6);

      mssm.set_input(pp);
      mssm.set_precision(1.0e-4); // == softsusy::TOLERANCE

      std::vector<Constraint<Two_scale>*> upward_constraints;
      upward_constraints.push_back(low_constraint);
      upward_constraints.push_back(high_constraint);

      std::vector<Constraint<Two_scale>*> downward_constraints;
      downward_constraints.push_back(high_constraint);
      downward_constraints.push_back(susy_constraint);
      downward_constraints.push_back(low_constraint);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&convergence_tester);
      solver.set_running_precision(&precision);
      solver.set_initial_guesser(&initial_guesser);
      solver.add_model(&mssm, upward_constraints, downward_constraints);
      solver.solve();
      mssm.run_to(susy_constraint->get_scale());
      mssm.calculate_spectrum();
      mssm.run_to(Electroweak_constants::MZ);

      mx = high_constraint->get_scale();
      msusy = susy_constraint->get_scale();
   }
private:
   double mx, msusy;
   NMSSM<Two_scale> mssm;
   NMSSM_high_scale_constraint<Two_scale>* high_constraint;
   NMSSM_susy_scale_constraint<Two_scale>* susy_constraint;
   NMSSM_low_scale_constraint<Two_scale>*  low_constraint;
};

BOOST_AUTO_TEST_CASE( test_NMSSM_spectrum )
{
   NMSSM_input_parameters pp;
   const NMSSM_high_scale_constraint<Two_scale> high_constraint(pp);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   NMSSM_tester nmssm_tester;
   BOOST_REQUIRE_NO_THROW(nmssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(nmssm_tester.get_mx(), softSusy_tester.get_mx(), 0.18);
   BOOST_CHECK_CLOSE_FRACTION(nmssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.006);

   // compare model parameters
   const NmssmSoftsusy ss(softSusy_tester.get_model());
   const NMSSM<Two_scale> fs(nmssm_tester.get_model());

   BOOST_CHECK_EQUAL(ss.displayLoops()     , fs.get_loops());
   BOOST_CHECK_EQUAL(ss.displayMu()        , fs.get_scale());
   BOOST_CHECK_EQUAL(ss.displayThresholds(), fs.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00076);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0015);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.00013);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Lambdax(), ss.displayLambda(), 0.0009);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Kappa()  , ss.displayKappa() , 0.007);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(0,0), ss.displayYukawaMatrix(YU)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(1,1), ss.displayYukawaMatrix(YU)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(2,2), ss.displayYukawaMatrix(YU)(3,3), 0.0012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(0,0), ss.displayYukawaMatrix(YD)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(1,1), ss.displayYukawaMatrix(YD)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(2,2), ss.displayYukawaMatrix(YD)(3,3), 0.0075);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(0,0), ss.displayYukawaMatrix(YE)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(1,1), ss.displayYukawaMatrix(YE)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(2,2), ss.displayYukawaMatrix(YE)(3,3), 0.0096);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassB() , ss.displayGaugino(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassWB(), ss.displayGaugino(2), 0.0046);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassG() , ss.displayGaugino(3), 0.0051);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.015);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.009);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ms2() , ss.displayMsSquared() , 0.015);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(0,0), ss.displaySoftMassSquared(mQl)(1,1), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(1,1), ss.displaySoftMassSquared(mQl)(2,2), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(2,2), ss.displaySoftMassSquared(mQl)(3,3), 0.012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(0,0), ss.displaySoftMassSquared(mUr)(1,1), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(1,1), ss.displaySoftMassSquared(mUr)(2,2), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(2,2), ss.displaySoftMassSquared(mUr)(3,3), 0.012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(0,0), ss.displaySoftMassSquared(mDr)(1,1), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(1,1), ss.displaySoftMassSquared(mDr)(2,2), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(2,2), ss.displaySoftMassSquared(mDr)(3,3), 0.012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(0,0), ss.displaySoftMassSquared(mLl)(1,1), 0.01);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(1,1), ss.displaySoftMassSquared(mLl)(2,2), 0.01);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(2,2), ss.displaySoftMassSquared(mLl)(3,3), 0.01);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(0,0), ss.displaySoftMassSquared(mEr)(1,1), 0.0025);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(1,1), ss.displaySoftMassSquared(mEr)(2,2), 0.0025);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(2,2), ss.displaySoftMassSquared(mEr)(3,3), 0.0013);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(0,0), ss.displayTrilinear(UA)(1,1), 0.018);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(1,1), ss.displayTrilinear(UA)(2,2), 0.018);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(2,2), ss.displayTrilinear(UA)(3,3), 0.008);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(0,0), ss.displayTrilinear(DA)(1,1), 0.019);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(1,1), ss.displayTrilinear(DA)(2,2), 0.019);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(2,2), ss.displayTrilinear(DA)(3,3), 0.016);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(0,0), ss.displayTrilinear(EA)(1,1), 0.021);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(1,1), ss.displayTrilinear(EA)(2,2), 0.021);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(2,2), ss.displayTrilinear(EA)(3,3), 0.021);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double vs = fs.get_vS();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 3.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(vs     , ss.displaySvev(), 0.004);
}
