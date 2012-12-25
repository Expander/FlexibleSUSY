#include "mssm_solver.h"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale_sugra_constraint.hpp"
#include "mssm_two_scale_msusy_constraint.hpp"
#include "mssm_two_scale_mz_constraint.hpp"
#include "mssm_two_scale_convergence_tester.hpp"
#include "two_scale_running_precision.hpp"
#include "softsusy.h"
#include "two_scale_solver.hpp"
#include "logger.hpp"
#include "stopwatch.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_mssm_solver

#include <boost/test/unit_test.hpp>

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

struct Mssm_parameter_point {
   Mssm_parameter_point()
      : m0(125.0)
      , m12(500.0)
      , a0(0.0)
      , mxGuess(1.0e16)
      , signMu(1)
      , tanBeta(10.0)
      , oneset()
   {
      const double alphasMZ = 0.1187, mtop = 173.4, mbmb = 4.2;
      oneset.setAlpha(ALPHAS, alphasMZ);
      oneset.setPoleMt(mtop);
      oneset.setMass(mBottom, mbmb);
      oneset.toMz();
   }
   DoubleVector get_soft_pars() const {
      DoubleVector highScaleSoftPars(3);
      highScaleSoftPars(1) = m0;
      highScaleSoftPars(2) = m12;
      highScaleSoftPars(3) = a0;
      return highScaleSoftPars;
   }
   friend std::ostream& operator<<(std::ostream& os, const Mssm_parameter_point& pp) {
      os << "CMSSM parameter point:"
         << " m0=" << pp.m0
         << ", m12=" << pp.m12
         << ", a0=" << pp.a0
         << ", mxGuess=" << pp.mxGuess
         << ", signMu=" << pp.signMu
         << ", tanBeta=" << pp.tanBeta
         << ", QedQcd=" << pp.oneset
         << '\n';
      return os;
   }

   double m0, m12, a0, mxGuess, signMu, tanBeta;
   QedQcd oneset;
};

/**
 * Test equality of physical MSSM parameters
 *
 * @param a first parameter set
 * @param b second parameter set
 * @param tolerance max. allowed deviation
 */
void test_equality(const sPhysical& a, const sPhysical& b, double tolerance)
{
   BOOST_CHECK_CLOSE(a.mh0         , b.mh0         , tolerance);
   BOOST_CHECK_CLOSE(a.mA0         , b.mA0         , tolerance);
   BOOST_CHECK_CLOSE(a.mH0         , b.mH0         , tolerance);
   BOOST_CHECK_CLOSE(a.mHpm        , b.mHpm        , tolerance);
   BOOST_CHECK_CLOSE(a.mGluino     , b.mGluino     , tolerance);
   BOOST_CHECK_CLOSE(a.thetaL      , b.thetaL      , tolerance);
   BOOST_CHECK_CLOSE(a.thetaR      , b.thetaR      , tolerance);
   BOOST_CHECK_CLOSE(a.thetat      , b.thetat      , tolerance);
   BOOST_CHECK_CLOSE(a.thetab      , b.thetab      , tolerance);
   BOOST_CHECK_CLOSE(a.thetatau    , b.thetatau    , tolerance);
   BOOST_CHECK_CLOSE(a.thetaH      , b.thetaH      , tolerance);
   BOOST_CHECK_CLOSE(a.t1OV1Ms     , b.t1OV1Ms     , tolerance);
   BOOST_CHECK_CLOSE(a.t2OV2Ms     , b.t2OV2Ms     , tolerance);
   BOOST_CHECK_CLOSE(a.t1OV1Ms1loop, b.t1OV1Ms1loop, tolerance);
   BOOST_CHECK_CLOSE(a.t2OV2Ms1loop, b.t2OV2Ms1loop, tolerance);

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

class SoftSusy_tester {
public:
   SoftSusy_tester()
      : mx(0.0), softSusy() {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   sPhysical get_physical() { return softSusy.displayPhys(); }
   void test(const Mssm_parameter_point& pp) {
      Stopwatch stopwatch;
      stopwatch.start(); // record time for SoftSusy to solve the MSSM

      // run softsusy
      softsusy::TOLERANCE = 1.0e-4;
#ifdef VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      mx = softSusy.lowOrg(sugraBcs, pp.mxGuess, pp.get_soft_pars(), pp.signMu, pp.tanBeta, pp.oneset, true);
      softsusy::PRINTOUT = 0;

      stopwatch.stop();
      VERBOSE_MSG("MssmSoftsusy solved in " << stopwatch.get_time_in_seconds()
                  << " seconds (" << stopwatch.get_clicks() << " clicks)");

      if (softSusy.displayProblem().test()) {
         std::stringstream ss;
         ss << "SoftSusy problem: " << softSusy.displayProblem();
         throw SoftSusy_error(ss.str());
      }
   }
private:
   double mx;
   MssmSoftsusy softSusy;
};

/**
 * Tests if our two scale algorithm calculates the same spectrum as
 * SoftSusy
 *
 * @param pp CMSSM parameter point
 */
void test_point(const Mssm_parameter_point& pp)
{
   // record time for the two scale method to solve the MSSM
   Stopwatch stopwatch;
   stopwatch.start();

   // setup the MSSM with the two scale method
   Mssm<Two_scale> mssm;
   Mssm_sugra_constraint mssm_sugra_constraint(&mssm, pp.mxGuess, pp.m0, pp.m12, pp.a0, pp.signMu);
   Mssm_mz_constraint mssm_mz_constraint(&mssm, pp.tanBeta);
   Mssm_msusy_constraint mssm_msusy_constraint(&mssm, pp.get_soft_pars(), 1000.0, pp.signMu);
   Mssm_convergence_tester mssm_convergence_tester(&mssm, 1.0e-4);
   Mssm_initial_guesser initial_guesser(&mssm, pp.oneset, pp.mxGuess, pp.tanBeta, pp.signMu, pp.get_soft_pars(), false);
   Two_scale_increasing_precision two_scale_increasing_precision(10.0, 1.0e-5);

   std::vector<Constraint<Two_scale>*> mssm_upward_constraints;
   mssm_upward_constraints.push_back(&mssm_mz_constraint);
   mssm_upward_constraints.push_back(&mssm_sugra_constraint);

   std::vector<Constraint<Two_scale>*> mssm_downward_constraints;
   mssm_downward_constraints.push_back(&mssm_sugra_constraint);
   mssm_downward_constraints.push_back(&mssm_msusy_constraint);
   mssm_downward_constraints.push_back(&mssm_mz_constraint);

   RGFlow<Two_scale> solver;
   solver.set_max_iterations(10);
   solver.set_convergence_tester(&mssm_convergence_tester);
   solver.set_running_precision(&two_scale_increasing_precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_model(&mssm, mssm_upward_constraints, mssm_downward_constraints);
   solver.solve();
   mssm.calculate_spectrum();

   stopwatch.stop();
   VERBOSE_MSG("Mssm<Two_scale> solved in " << stopwatch.get_time_in_seconds()
               << " seconds (" << stopwatch.get_clicks() << " clicks)");

   SoftSusy_tester softSusy_tester;
   softSusy_tester.test(pp);

   // check equality of physical parameters
   test_equality(softSusy_tester.get_physical(), mssm.displayPhys(), 0.1);
   BOOST_CHECK_CLOSE(softSusy_tester.get_mx(), mssm_sugra_constraint.get_scale(), 0.1);
}

BOOST_AUTO_TEST_CASE( test_default_cmssm_parameter_point )
{
   Mssm_parameter_point pp;
   BOOST_MESSAGE("testing parameter point " << pp);
   test_point(pp);
}

BOOST_AUTO_TEST_CASE( test_slow_convergence_point )
{
   // slowly convergent point taken from arXiv:1211.3231 Fig. 7
   Mssm_parameter_point pp;
   pp.oneset.setPoleMt(173.5);
   pp.tanBeta = 10.0;
   pp.a0 = 0.0;
   pp.m12 = 337.5;
   pp.m0 = 3400.0;

   BOOST_MESSAGE("testing slow convergent parameter point " << pp);
   BOOST_CHECK_THROW(test_point(pp), RGFlow<Two_scale>::NoConvergenceError);
}

BOOST_AUTO_TEST_CASE( test_non_perturbative_point )
{
   // non-perturbative point
   Mssm_parameter_point pp;
   pp.oneset.setPoleMt(173.5);
   pp.tanBeta = 100.0;
   pp.a0 = 0.0;
   pp.m12 = 337.5;
   pp.m0 = 3400.0;

   BOOST_MESSAGE("testing non-perturbative parameter point " << pp);
   BOOST_CHECK_THROW(test_point(pp), RGFlow<Two_scale>::NonPerturbativeRunningError);
}
