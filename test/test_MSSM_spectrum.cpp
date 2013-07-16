
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_spectrum

#include <boost/test/unit_test.hpp>

#include "softsusy.h"
#include "error.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "MSSM_model.hpp"
#include "MSSM_input_parameters.hpp"
#include "MSSM_high_scale_constraint.hpp"
#include "MSSM_susy_scale_constraint.hpp"
#include "MSSM_low_scale_constraint.hpp"
#include "MSSM_convergence_tester.hpp"
#include "MSSM_initial_guesser.hpp"

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
      : mx(0.0), msusy(0.0), softSusy() {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   sPhysical get_physical() const { return softSusy.displayPhys(); }
   MssmSoftsusy get_model() const { return softSusy; }
   void test(const MSSM_input_parameters& pp, double mxGuess, const QedQcd& oneset = QedQcd()) {
      // run softsusy
      softsusy::TOLERANCE = 1.0e-4;
#ifdef VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      DoubleVector pars(3);
      pars(1) = pp.m0;
      pars(2) = pp.m12;
      pars(3) = pp.Azero;
      mx = softSusy.lowOrg(sugraBcs, mxGuess, pars, pp.SignMu, pp.TanBeta, oneset, true);
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
   MssmSoftsusy softSusy;
};

class MSSM_tester {
public:
   MSSM_tester()
      : mx(0.0), msusy(0.0), mssm() {}
   ~MSSM_tester() {}
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   MSSM_physical get_physical() const { return mssm.get_physical(); }
   MSSM get_model() const { return mssm; }
   void test(const MSSM_input_parameters& pp) {
      MSSM_high_scale_constraint sugra_constraint(pp);
      MSSM_low_scale_constraint  low_constraint(pp);
      MSSM_susy_scale_constraint susy_constraint(pp);
      MSSM_convergence_tester    convergence_tester(&mssm, 1.0e-4);
      MSSM_initial_guesser initial_guesser(&mssm, pp, low_constraint,
                                           susy_constraint,
                                           sugra_constraint);
      Two_scale_increasing_precision precision(10.0, 1.0e-5);

      mssm.set_input(pp);

      std::vector<Constraint<Two_scale>*> upward_constraints;
      upward_constraints.push_back(&low_constraint);
      upward_constraints.push_back(&sugra_constraint);

      std::vector<Constraint<Two_scale>*> downward_constraints;
      downward_constraints.push_back(&sugra_constraint);
      downward_constraints.push_back(&susy_constraint);
      downward_constraints.push_back(&low_constraint);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&convergence_tester);
      solver.set_running_precision(&precision);
      solver.set_initial_guesser(&initial_guesser);
      solver.add_model(&mssm, upward_constraints, downward_constraints);
      solver.solve();
      mssm.run_to(susy_constraint.get_scale());
      mssm.calculate_spectrum();
      mssm.run_to(Electroweak_constants::MZ);

      mx = sugra_constraint.get_scale();
      msusy = susy_constraint.get_scale();
   }
private:
   double mx, msusy;
   MSSM mssm;
};

BOOST_AUTO_TEST_CASE( test_MSSM_GUT_scale )
{
   MSSM_input_parameters pp;
   const MSSM_high_scale_constraint high_constraint(pp);
   const double mxGuess = high_constraint.get_scale();

   MSSM_tester mssm_tester;
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 0.14);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.015);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const MSSM fs(mssm_tester.get_model());

   BOOST_CHECK_EQUAL(ss.displayLoops()     , fs.get_loops());
   BOOST_CHECK_EQUAL(ss.displayMu()        , fs.get_scale());
   BOOST_CHECK_EQUAL(ss.displayThresholds(), fs.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00027);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0013);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 8.0e-5);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassB() , ss.displayGaugino(1), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassWB(), ss.displayGaugino(2), 0.002);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassG() , ss.displayGaugino(3), 0.0022);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.009);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0012);
}
