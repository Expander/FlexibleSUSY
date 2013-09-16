#include "mssm_solver.h"
#include "mssm_parameter_point.hpp"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale_sugra_constraint.hpp"
#include "mssm_two_scale_susy_scale_constraint.hpp"
#include "mssm_two_scale_low_scale_constraint.hpp"
#include "mssm_two_scale_convergence_tester.hpp"
#include "two_scale_running_precision.hpp"
#include "softsusy.h"
#include "two_scale_solver.hpp"
#include "logger.hpp"
#include "stopwatch.hpp"
#include "error.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_scale_mssm_solver

#include <boost/test/unit_test.hpp>

using namespace flexiblesusy;

BOOST_AUTO_TEST_CASE( test_softsusy_mssm_solver )
{
   DoubleVector highScaleSoftPars(3);
   const double m12 = 500., a0 = 0., m0 = 125.;
   highScaleSoftPars(1) = m0;
   highScaleSoftPars(2) = m12;
   highScaleSoftPars(3) = a0;

   const int signMu = 1;
   const double tanBeta = 10.0;
   const bool uni = true;
   const double mxGuess = 1.0e16;

   QedQcd oneset;
   const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
   oneset.setAlpha(ALPHAS, alphasMZ);
   oneset.setPoleMt(mtop);
   oneset.setMass(mBottom, mbmb);
   oneset.toMz();

   RGFlow<SoftSusy_t> mssmSolver;
   mssmSolver.setSoftHighScalePars(highScaleSoftPars);
   mssmSolver.setSignMu(signMu);
   mssmSolver.setTanBeta(tanBeta);
   mssmSolver.setLowScaleBoundaryConditions(oneset);
   mssmSolver.setGaugeUnification(uni);
   mssmSolver.setMxGuess(mxGuess);
   mssmSolver.setHighScaleBoundaryCondition(sugraBcs);
   mssmSolver.solve();

   MssmSoftsusy softSusy;
   softSusy.lowOrg(sugraBcs, mxGuess, highScaleSoftPars, signMu, tanBeta, oneset, uni);

   BOOST_CHECK_EQUAL(mssmSolver.displayPhys(), softSusy.displayPhys());
}

/**
 * Test equality of physical MSSM parameters
 *
 * @param a first parameter set
 * @param b second parameter set
 * @param tolerance max. allowed deviation
 */
void test_equality(const sPhysical& a, const sPhysical& b, double tolerance)
{
   BOOST_CHECK_CLOSE(a.mh0(1)      , b.mh0(1)      , tolerance);
   BOOST_CHECK_CLOSE(a.mA0(1)      , b.mA0(1)      , tolerance);
   BOOST_CHECK_CLOSE(a.mh0(2)      , b.mh0(2)      , tolerance);
   BOOST_CHECK_CLOSE(a.mHpm        , b.mHpm        , tolerance);
   BOOST_CHECK_CLOSE(a.mGluino     , b.mGluino     , tolerance);
   BOOST_CHECK_CLOSE(a.thetaL      , b.thetaL      , tolerance);
   BOOST_CHECK_CLOSE(a.thetaR      , b.thetaR      , tolerance);
   BOOST_CHECK_CLOSE(a.thetat      , b.thetat      , tolerance);
   BOOST_CHECK_CLOSE(a.thetab      , b.thetab      , tolerance);
   BOOST_CHECK_CLOSE(a.thetatau    , b.thetatau    , tolerance);
   BOOST_CHECK_CLOSE(a.thetaH      , b.thetaH      , tolerance);

   // sneutrino masses
   for (int i = a.msnu.displayStart();
        i < a.msnu.displayEnd(); ++i)
      BOOST_CHECK_CLOSE(a.msnu(i), b.msnu(i), tolerance);

   // chargino masses
   for (int i = a.mch.displayStart();
        i < a.mch.displayEnd(); ++i)
      BOOST_CHECK_CLOSE(a.mch(i), b.mch(i), tolerance);

   // neutralino masses
   for (int i = a.mneut.displayStart();
        i < a.mneut.displayEnd(); ++i)
      BOOST_CHECK_CLOSE(a.mneut(i), b.mneut(i), tolerance);

   // neuralino mixing matrix
   for (int i = 1; i < a.mixNeut.displayRows(); ++i)
      for (int k = 1; k < a.mixNeut.displayCols(); ++k)
         BOOST_CHECK_CLOSE(a.mixNeut(i,k), b.mixNeut(i,k), tolerance);

   // up squarks
   for (int i = 1; i < a.mu.displayRows(); ++i)
      for (int k = 1; k < a.mu.displayCols(); ++k)
         BOOST_CHECK_CLOSE(a.mu(i,k), b.mu(i,k), tolerance);

   // down squarks
   for (int i = 1; i < a.md.displayRows(); ++i)
      for (int k = 1; k < a.md.displayCols(); ++k)
         BOOST_CHECK_CLOSE(a.md(i,k), b.md(i,k), tolerance);

   // down sleptons
   for (int i = 1; i < a.me.displayRows(); ++i)
      for (int k = 1; k < a.me.displayCols(); ++k)
         BOOST_CHECK_CLOSE(a.me(i,k), b.me(i,k), tolerance);
}

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
      : mx(0.0), softSusy() {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   sPhysical get_physical() const { return softSusy.displayPhys(); }
   void test(const Mssm_parameter_point& pp, const QedQcd& oneset) {
      Stopwatch stopwatch;
      stopwatch.start(); // record time for SoftSusy to solve the MSSM

      // run softsusy
      softsusy::TOLERANCE = 1.0e-4;
#ifdef VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      softSusy.lowOrg(sugraBcs, pp.mxGuess, pp.get_soft_pars(), pp.signMu, pp.tanBeta, oneset, true);
      mx = softSusy.displayMxBC();
      softsusy::PRINTOUT = 0;

      stopwatch.stop();
      VERBOSE_MSG("MssmSoftsusy solved in " << stopwatch.get_time_in_seconds()
                  << " seconds");

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
   double mx;
   MssmSoftsusy softSusy;
};

class Two_scale_tester {
public:
   Two_scale_tester()
      : mx(0.0), mssm() {}
   ~Two_scale_tester() {}
   double get_mx() const { return mx; }
   sPhysical get_physical() const { return mssm.displayPhys(); }
   void test(const Mssm_parameter_point& pp, const QedQcd& oneset) {
      Stopwatch stopwatch;
      stopwatch.start(); // record time for the two scale method to solve the MSSM

      softsusy::TOLERANCE = 1.0e-4;
      // setup the MSSM with the two scale method
      Mssm_sugra_constraint mssm_sugra_constraint(pp);
      Mssm_low_scale_constraint mssm_mz_constraint(pp);
      Mssm_susy_scale_constraint mssm_msusy_constraint(pp);
      Mssm_convergence_tester mssm_convergence_tester(&mssm, 1.0e-4);
      Mssm_initial_guesser initial_guesser(&mssm, pp, mssm_mz_constraint,
                                           mssm_msusy_constraint,
                                           mssm_sugra_constraint);
      initial_guesser.set_QedQcd(oneset);
      Two_scale_increasing_precision two_scale_increasing_precision(10.0, 1.0e-5);

      std::vector<Constraint<Two_scale>*> mssm_upward_constraints;
      mssm_upward_constraints.push_back(&mssm_mz_constraint);
      mssm_upward_constraints.push_back(&mssm_sugra_constraint);

      std::vector<Constraint<Two_scale>*> mssm_downward_constraints;
      mssm_downward_constraints.push_back(&mssm_sugra_constraint);
      mssm_downward_constraints.push_back(&mssm_msusy_constraint);
      mssm_downward_constraints.push_back(&mssm_mz_constraint);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&mssm_convergence_tester);
      solver.set_running_precision(&two_scale_increasing_precision);
      solver.set_initial_guesser(&initial_guesser);
      solver.add_model(&mssm, mssm_upward_constraints, mssm_downward_constraints);
      solver.solve();
      mssm.calculate_spectrum();

      stopwatch.stop();
      VERBOSE_MSG("Mssm<Two_scale> solved in " << stopwatch.get_time_in_seconds()
                  << " seconds");

      mx = mssm_sugra_constraint.get_scale();
   }
private:
   double mx;
   Mssm<Two_scale> mssm;
};

/**
 * Tests if our two scale algorithm calculates the same spectrum as
 * SoftSusy
 *
 * @param pp CMSSM parameter point
 */
void test_point(const Mssm_parameter_point& pp)
{
   QedQcd oneset;

   Two_scale_tester two_scale_tester;
   BOOST_REQUIRE_NO_THROW(two_scale_tester.test(pp, oneset));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, oneset));

   // check equality of physical parameters
   test_equality(softSusy_tester.get_physical(), two_scale_tester.get_physical(), 0.1);
   BOOST_CHECK_CLOSE(softSusy_tester.get_mx(), two_scale_tester.get_mx(), 0.1);
}

BOOST_AUTO_TEST_CASE( test_default_cmssm_parameter_point )
{
   Mssm_parameter_point pp;
   BOOST_MESSAGE("testing " << pp);
   test_point(pp);
}

BOOST_AUTO_TEST_CASE( test_cmssm_tanb_scan )
{
   // do small parameter scan of tan(beta)
   Mssm_parameter_point pp;
   for (double tanb = 3.0; tanb <= 44.; tanb += 3.0) {
      pp.tanBeta = tanb;
      BOOST_MESSAGE("testing " << pp);
      BOOST_CHECK_NO_THROW(test_point(pp));
   }
}

BOOST_AUTO_TEST_CASE( test_slow_convergence_point )
{
   // slowly convergent point taken from arXiv:1211.3231 Fig. 7
   Mssm_parameter_point pp;
   pp.tanBeta = 10.0;
   pp.a0 = 0.0;
   pp.m12 = 337.5;
   pp.m0 = 3400.0;
   QedQcd oneset;
   const double alphasMZ = 0.1187, mtop = 173.5, mbmb = 4.2;
   oneset.setAlpha(ALPHAS, alphasMZ);
   oneset.setPoleMt(mtop);
   oneset.setMass(mBottom, mbmb);
   oneset.toMz();

   BOOST_MESSAGE("testing slow convergent " << pp);
   Two_scale_tester two_scale_tester;
   BOOST_CHECK_THROW(two_scale_tester.test(pp, oneset), NoConvergenceError);
   SoftSusy_tester softSusy_tester;
   // BOOST_CHECK_THROW(softSusy_tester.test(pp, oneset), SoftSusy_NoConvergence_error);
   try {
      softSusy_tester.test(pp, oneset);
   } catch (SoftSusy_NoConvergence_error) {
      INFO("SoftSusy_NoConvergence_error thrown");
   } catch (std::string& s) {
      INFO("string: " << s);
   } catch (const char* s) {
      INFO("char: " << s);
   } catch (...) {
      BOOST_FAIL("Softsusy has thrown an unknown exception");
   }
}

BOOST_AUTO_TEST_CASE( test_non_perturbative_point )
{
   // non-perturbative point
   Mssm_parameter_point pp;
   pp.tanBeta = 100.0;
   pp.a0 = 0.0;
   pp.m12 = 337.5;
   pp.m0 = 3400.0;
   QedQcd oneset;
   oneset.setPoleMt(173.5);

   BOOST_MESSAGE("testing non-perturbative " << pp);
   Two_scale_tester two_scale_tester;
   BOOST_CHECK_THROW(two_scale_tester.test(pp, oneset), NonPerturbativeRunningError);
   SoftSusy_tester softSusy_tester;
   // BOOST_CHECK_THROW(softSusy_tester.test(pp, oneset), SoftSusy_NonPerturbative_error);
   BOOST_CHECK_THROW(softSusy_tester.test(pp, oneset), SoftSusy_error);
}
