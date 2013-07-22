
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
      softsusy::numRewsbLoops = 1;
      softsusy::TOLERANCE = 1.0e-4;
#ifdef VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      DoubleVector pars(3);
      pars(1) = pp.m0;
      pars(2) = pp.m12;
      pars(3) = pp.Azero;
      softSusy.setAlternativeMs(true);
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
      Two_scale_increasing_precision precision(10.0, 1.0e-6);

      mssm.set_input(pp);
      mssm.set_precision(1.0e-4); // == softsusy::TOLERANCE

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
   const double mxGuess = high_constraint.get_initial_scale_guess();

   MSSM_tester mssm_tester;
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 0.14);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.003);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const MSSM fs(mssm_tester.get_model());

   BOOST_CHECK_EQUAL(ss.displayLoops()     , fs.get_loops());
   BOOST_CHECK_EQUAL(ss.displayMu()        , fs.get_scale());
   BOOST_CHECK_EQUAL(ss.displayThresholds(), fs.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00075);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0028);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.00011);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(0,0), ss.displayYukawaMatrix(YU)(1,1), 0.0012);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(1,1), ss.displayYukawaMatrix(YU)(2,2), 0.0012);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(2,2), ss.displayYukawaMatrix(YU)(3,3), 0.0012);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(0,0), ss.displayYukawaMatrix(YD)(1,1), 0.0012);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(1,1), ss.displayYukawaMatrix(YD)(2,2), 0.0012);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(2,2), ss.displayYukawaMatrix(YD)(3,3), 0.0075);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(0,0), ss.displayYukawaMatrix(YE)(1,1), 0.0012);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(1,1), ss.displayYukawaMatrix(YE)(2,2), 0.0012);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(2,2), ss.displayYukawaMatrix(YE)(3,3), 0.0016);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_MassB() , ss.displayGaugino(1), 0.005);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_MassWB(), ss.displayGaugino(2), 0.002);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_MassG() , ss.displayGaugino(3), 0.0022);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.009);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.009);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(0,0), ss.displaySoftMassSquared(mQl)(1,1), 0.0051);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(1,1), ss.displaySoftMassSquared(mQl)(2,2), 0.0051);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_mq2()(2,2), ss.displaySoftMassSquared(mQl)(3,3), 0.0044);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(0,0), ss.displaySoftMassSquared(mUr)(1,1), 0.0051);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(1,1), ss.displaySoftMassSquared(mUr)(2,2), 0.0051);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_mu2()(2,2), ss.displaySoftMassSquared(mUr)(3,3), 0.0037);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(0,0), ss.displaySoftMassSquared(mDr)(1,1), 0.0051);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(1,1), ss.displaySoftMassSquared(mDr)(2,2), 0.0051);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_md2()(2,2), ss.displaySoftMassSquared(mDr)(3,3), 0.0047);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(0,0), ss.displaySoftMassSquared(mLl)(1,1), 0.0043);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(1,1), ss.displaySoftMassSquared(mLl)(2,2), 0.0043);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_ml2()(2,2), ss.displaySoftMassSquared(mLl)(3,3), 0.0042);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(0,0), ss.displaySoftMassSquared(mEr)(1,1), 0.0011);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(1,1), ss.displaySoftMassSquared(mEr)(2,2), 0.0011);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_me2()(2,2), ss.displaySoftMassSquared(mEr)(3,3), 0.00079);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(0,0), ss.displayTrilinear(UA)(1,1), 0.0043);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(1,1), ss.displayTrilinear(UA)(2,2), 0.0043);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYu()(2,2), ss.displayTrilinear(UA)(3,3), 0.0027);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(0,0), ss.displayTrilinear(DA)(1,1), 0.005);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(1,1), ss.displayTrilinear(DA)(2,2), 0.005);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYd()(2,2), ss.displayTrilinear(DA)(3,3), 0.011);

   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(0,0), ss.displayTrilinear(EA)(1,1), 0.0057);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(1,1), ss.displayTrilinear(EA)(2,2), 0.0057);
   // BOOST_CHECK_CLOSE_FRACTION(fs.get_TYe()(2,2), ss.displayTrilinear(EA)(3,3), 0.006);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0071);

   // comparing tree-level masses

   const DoubleVector MCha(fs.get_MCha()), MChi(fs.get_MChi()),
      MHpm(fs.get_MHpm()), MAh(fs.get_MAh()), Mhh(fs.get_Mhh());
   const DoubleVector mch(ss.displayDrBarPars().mchBpmz),
      mn(ss.displayDrBarPars().mnBpmz);
   const double MwRun = ss.displayMwRun();
   const double MzRun = ss.displayMzRun();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const double mA0 = ss.displayDrBarPars().mA0;
   const double mh0 = ss.displayDrBarPars().mh0;
   const double mH0 = ss.displayDrBarPars().mH0;

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha(1), mch(1), 0.004);
   BOOST_CHECK_CLOSE_FRACTION(MCha(2), mch(2), 0.02);

   BOOST_CHECK_CLOSE_FRACTION(MChi(1), mn(1), 0.0044);
   BOOST_CHECK_CLOSE_FRACTION(MChi(2), mn(2), 0.0041);
   BOOST_CHECK_CLOSE_FRACTION(MChi(3), mn(3), 0.021);
   BOOST_CHECK_CLOSE_FRACTION(MChi(4), mn(4), 0.020);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 0.018); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm , 0.04);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 0.014); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0  , 0.041);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0, 0.0004);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mH0, 0.041);

   // down-type squarks
   const DoubleVector Sd(fs.get_MSd());
   const DoubleVector md(ss.displayDrBarPars().md.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Sd(1), md(1), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Sd(2), md(2), 0.0032);
   BOOST_CHECK_CLOSE_FRACTION(Sd(3), md(3), 0.0034);
   BOOST_CHECK_CLOSE_FRACTION(Sd(4), md(4), 0.0034);
   BOOST_CHECK_CLOSE_FRACTION(Sd(5), md(5), 0.0034);
   BOOST_CHECK_CLOSE_FRACTION(Sd(6), md(6), 0.0034);

   // up-type squarks
   const DoubleVector Su(fs.get_MSu());
   const DoubleVector mu(ss.displayDrBarPars().mu.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Su(1), mu(1), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Su(2), mu(2), 0.0034);
   BOOST_CHECK_CLOSE_FRACTION(Su(3), mu(3), 0.0034);
   BOOST_CHECK_CLOSE_FRACTION(Su(4), mu(4), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Su(5), mu(5), 0.0034);
   BOOST_CHECK_CLOSE_FRACTION(Su(6), mu(6), 0.0034);

   // sleptons
   const DoubleVector Se(fs.get_MSe());
   const DoubleVector me(ss.displayDrBarPars().me.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Se(1), me(1), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(2), me(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(3), me(3), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(4), me(4), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(5), me(5), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(6), me(6), 0.003);

   // sneutrinos
   const DoubleVector msnu(ss.displayDrBarPars().msnu);
   const DoubleVector Snu(fs.get_MSv());
   BOOST_CHECK_CLOSE_FRACTION(Snu(1), msnu(1), 0.0055);
   BOOST_CHECK_CLOSE_FRACTION(Snu(2), msnu(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Snu(3), msnu(3), 0.003);

   BOOST_CHECK_EQUAL(fs.get_MVP(), 0.0);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MVZ() , MzRun, 0.02);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MVWm(), MwRun, 0.022);

   BOOST_CHECK_EQUAL(fs.get_MVG(), 0.0);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MGlu(), ss.displayDrBarPars().mGluino, 0.005);

   BOOST_CHECK_EQUAL(fs.get_MFv()(1), 0.0);
   BOOST_CHECK_EQUAL(fs.get_MFv()(2), 0.0);
   BOOST_CHECK_EQUAL(fs.get_MFv()(3), 0.0);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFe()(3), ss.displayDrBarPars().mtau, 0.00044);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFu()(3), ss.displayDrBarPars().mt  , 0.00047);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFd()(3), ss.displayDrBarPars().mb  , 0.0067);
}
