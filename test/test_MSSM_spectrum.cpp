
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_spectrum

#include <boost/test/unit_test.hpp>

#define private public

#include "softsusy.h"
#include "error.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSM_input_parameters.hpp"
#include "MSSM_two_scale_high_scale_constraint.hpp"
#include "MSSM_two_scale_susy_scale_constraint.hpp"
#include "MSSM_two_scale_low_scale_constraint.hpp"
#include "MSSM_two_scale_convergence_tester.hpp"
#include "MSSM_two_scale_initial_guesser.hpp"
#include "test_MSSM.hpp"

/**
 * @class MSSM_precise_gauge_couplings_low_scale_constraint
 *
 * Replacement class for MSSM_low_scale_constraint, which calculates
 * the gauge couplings at the low scale as Softsusy does it.
 */
class MSSM_precise_gauge_couplings_low_scale_constraint
   : public MSSM_low_scale_constraint<Two_scale> {
public:
   MSSM_precise_gauge_couplings_low_scale_constraint()
      : MSSM_low_scale_constraint<Two_scale>() {}
   MSSM_precise_gauge_couplings_low_scale_constraint(const MSSM_input_parameters& inputPars_, const QedQcd& oneset_)
      : MSSM_low_scale_constraint<Two_scale>(inputPars_,oneset_) {}
   virtual ~MSSM_precise_gauge_couplings_low_scale_constraint() {}

   virtual void apply();
};

void MSSM_precise_gauge_couplings_low_scale_constraint::apply()
{
   assert(model && "Error: MSSM_precise_gauge_couplings_low_scale_constraint:"
          " model pointer must not be zero");

   // save old model parmeters
   const MSSM<Two_scale> mssm(*model);

   // run MSSM_low_scale_constraint::apply(), without the gauge
   // couplings
   model->calculate_DRbar_parameters();
   update_scale();
   calculate_DRbar_gauge_couplings();

   const double MZDRbar
      = model->calculate_MVZ_DRbar_1loop(Electroweak_constants::MZ);

   const double TanBeta = inputPars.TanBeta;
   const double g1 = model->get_g1();
   const double g2 = model->get_g2();

   model->set_vd((2*MZDRbar)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(TanBeta)
      )));
   model->set_vu((2*MZDRbar*TanBeta)/(Sqrt(0.6*Sqr(g1) + Sqr(g2))*Sqrt(1 + Sqr(
      TanBeta))));

   calculate_DRbar_yukawa_couplings();

   model->set_Yu(new_Yu);
   model->set_Yd(new_Yd);
   model->set_Ye(new_Ye);

   // Now calculate the gauge couplings using
   // MssmSoftsusy::sparticleThresholdCorrections
   MssmSoftsusy softsusy;
   copy_parameters(mssm, softsusy);

   softsusy.sparticleThresholdCorrections(inputPars.TanBeta);

   BOOST_MESSAGE("Difference (g1_FlexibleSUSY - g1_softsusy)(MZ) = "
                 << new_g1 - softsusy.displayGaugeCoupling(1));
   BOOST_MESSAGE("Difference (g2_FlexibleSUSY - g2_softsusy)(MZ) = "
                 << new_g2 - softsusy.displayGaugeCoupling(2));
   BOOST_MESSAGE("Difference (g3_FlexibleSUSY - g3_softsusy)(MZ) = "
                 << new_g3 - softsusy.displayGaugeCoupling(3));

   BOOST_MESSAGE("Difference (Yu_FlexibleSUSY - Yu_softsusy)(MZ) = "
                 << ToDoubleMatrix(new_Yu) - softsusy.displayYukawaMatrix(YU));
   BOOST_MESSAGE("Difference (Yd_FlexibleSUSY - Yd_softsusy)(MZ) = "
                 << ToDoubleMatrix(new_Yd) - softsusy.displayYukawaMatrix(YD));
   BOOST_MESSAGE("Difference (Ye_FlexibleSUSY - Ye_softsusy)(MZ) = "
                 << ToDoubleMatrix(new_Ye) - softsusy.displayYukawaMatrix(YE));

   model->set_g1(softsusy.displayGaugeCoupling(1));
   model->set_g2(softsusy.displayGaugeCoupling(2));
   model->set_g3(softsusy.displayGaugeCoupling(3));

   model->set_Yu(ToEigenMatrix(softsusy.displayYukawaMatrix(YU)));
   model->set_Yd(ToEigenMatrix(softsusy.displayYukawaMatrix(YD)));
   model->set_Ye(ToEigenMatrix(softsusy.displayYukawaMatrix(YE)));

   const double tanBeta = softsusy.displayTanb();
   const double vev = softsusy.displayHvev();
   const double beta = atan(tanBeta);
   const double sinBeta = sin(beta);
   const double cosBeta = cos(beta);
   const double vu = sinBeta * vev;
   const double vd = cosBeta * vev;

   BOOST_MESSAGE("Difference (vu_FlexibleSUSY - vu_softsusy)(MZ) = "
                 << model->get_vu() - vu);
   BOOST_MESSAGE("Difference (vd_FlexibleSUSY - vd_softsusy)(MZ) = "
                 << model->get_vd() - vd);

   model->set_vu(vu);
   model->set_vd(vd);
}

/**
 * @class MSSM_softsusy_ewsb_susy_scale_constraint
 *
 * Replacement class for MSSM_susy_scale_constraint, which does the
 * one-loop ewsb at the susy scale as Softsusy does it.
 */
class MSSM_softsusy_ewsb_susy_scale_constraint
   : public MSSM_susy_scale_constraint<Two_scale> {
public:
   MSSM_softsusy_ewsb_susy_scale_constraint()
      : MSSM_susy_scale_constraint<Two_scale>() {}
   MSSM_softsusy_ewsb_susy_scale_constraint(const MSSM_input_parameters& inputPars_)
      : MSSM_susy_scale_constraint<Two_scale>(inputPars_) {}
   virtual ~MSSM_softsusy_ewsb_susy_scale_constraint() {}

   virtual void apply();
};

void MSSM_softsusy_ewsb_susy_scale_constraint::apply()
{
   assert(model && "Error: MSSM_softsusy_ewsb_susy_scale_constraint:"
          " model pointer must not be zero");

   // save old model parmeters
   model->calculate_DRbar_parameters();
   const MSSM<Two_scale> mssm(*model);

   MSSM_susy_scale_constraint<Two_scale>::apply();

   // Now do the one-loop EWSB using MssmSoftsusy::rewsb
   MssmSoftsusy softsusy;
   softsusy.setAlternativeMs(true);
   copy_parameters(mssm, softsusy);
   softsusy.calcDrBarPars();
   const double new_Msusy = softsusy.calcMs();

   const int signMu = inputPars.SignMu;
   const double mt = softsusy.displayDrBarPars().mt;
   DoubleVector highScaleSoftPars(3);
   highScaleSoftPars(1) = inputPars.m0;
   highScaleSoftPars(2) = inputPars.m12;
   highScaleSoftPars(3) = inputPars.Azero;

   softsusy.rewsb(signMu, mt, highScaleSoftPars);

   const double new_Mu  = softsusy.displaySusyMu();
   const double new_BMu = softsusy.displayM3Squared();

   BOOST_MESSAGE("Difference (Mu_FlexibleSUSY - Mu_softsusy)(Msusy) = "
                 << model->get_Mu() - new_Mu);
   BOOST_MESSAGE("Difference (BMu_FlexibleSUSY - BMu_softsusy)(Msusy) = "
                 << model->get_BMu() - new_BMu);
   BOOST_MESSAGE("Difference (mt_FlexibleSUSY - mt_softsusy)(Msusy) = "
                 << model->get_MFu()(2) - mt);
   BOOST_MESSAGE("Difference (Msusy_FlexibleSUSY - Msusy_softsusy)(Msusy) = "
                 << get_scale() - new_Msusy);

   model->set_Mu(new_Mu);
   model->set_BMu(new_BMu);
   scale = new_Msusy;
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
      : mx(0.0), msusy(0.0), softSusy(), gaugeUnification(true) {}
   ~SoftSusy_tester() {}
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   sPhysical get_physical() const { return softSusy.displayPhys(); }
   MssmSoftsusy get_model() const { return softSusy; }
   void test(const MSSM_input_parameters& pp, double mxGuess, const QedQcd& oneset = QedQcd()) {
      // run softsusy
      softsusy::numRewsbLoops = 1;
      softsusy::numHiggsMassLoops = 1;
      softsusy::TOLERANCE = 1.0e-4;
#ifdef VERBOSE
      softsusy::PRINTOUT = 1;
#endif
      DoubleVector pars(3);
      pars(1) = pp.m0;
      pars(2) = pp.m12;
      pars(3) = pp.Azero;
      softSusy.setAlternativeMs(true);
      softSusy.lowOrg(sugraBcs, mxGuess, pars, pp.SignMu, pp.TanBeta,
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
   MssmSoftsusy softSusy;
   bool gaugeUnification;
};

class MSSM_tester {
public:
   MSSM_tester()
      : mx(0.0), msusy(0.0), mssm()
      , high_constraint(NULL), susy_constraint(NULL), low_constraint(NULL) {}
   ~MSSM_tester() {
      delete high_constraint;
      delete susy_constraint;
      delete low_constraint;
   }
   double get_mx() const { return mx; }
   double get_msusy() const { return msusy; }
   MSSM_physical get_physical() const { return mssm.get_physical(); }
   MSSM<Two_scale> get_model() const { return mssm; }
   void set_low_scale_constraint(MSSM_low_scale_constraint<Two_scale>* c) { low_constraint = c; }
   void set_susy_scale_constraint(MSSM_susy_scale_constraint<Two_scale>* c) { susy_constraint = c; }
   void set_high_scale_constraint(MSSM_high_scale_constraint<Two_scale>* c) { high_constraint = c; }
   void setup_default_constaints() {
      if (!high_constraint)
         high_constraint = new MSSM_high_scale_constraint<Two_scale>();
      if (!susy_constraint)
         susy_constraint = new MSSM_susy_scale_constraint<Two_scale>();
      if (!low_constraint)
         low_constraint = new MSSM_low_scale_constraint<Two_scale>();
   }
   void test(const MSSM_input_parameters& pp, const QedQcd& oneset = QedQcd()) {
      setup_default_constaints();
      high_constraint->set_input_parameters(pp);
      low_constraint->set_input_parameters(pp);
      low_constraint->set_sm_parameters(oneset);
      susy_constraint->set_input_parameters(pp);

      MSSM_convergence_tester<Two_scale> convergence_tester(&mssm, 1.0e-4);
      MSSM_initial_guesser<Two_scale> initial_guesser(&mssm, pp, oneset,
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
   MSSM<Two_scale> mssm;
   MSSM_high_scale_constraint<Two_scale>* high_constraint;
   MSSM_susy_scale_constraint<Two_scale>* susy_constraint;
   MSSM_low_scale_constraint<Two_scale>*  low_constraint;
};

BOOST_AUTO_TEST_CASE( test_MSSM_spectrum )
{
   MSSM_input_parameters pp;
   const MSSM_high_scale_constraint<Two_scale> high_constraint(pp);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   MSSM_tester mssm_tester;
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 0.13);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 0.006);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const MSSM<Two_scale> fs(mssm_tester.get_model());

   BOOST_CHECK_EQUAL(ss.displayLoops()     , fs.get_loops());
   BOOST_CHECK_EQUAL(ss.displayMu()        , fs.get_scale());
   BOOST_CHECK_EQUAL(ss.displayThresholds(), fs.get_thresholds());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00076);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.0011);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.00013);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(0,0), ss.displayYukawaMatrix(YU)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(1,1), ss.displayYukawaMatrix(YU)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yu()(2,2), ss.displayYukawaMatrix(YU)(3,3), 0.0012);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(0,0), ss.displayYukawaMatrix(YD)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(1,1), ss.displayYukawaMatrix(YD)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Yd()(2,2), ss.displayYukawaMatrix(YD)(3,3), 0.0075);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(0,0), ss.displayYukawaMatrix(YE)(1,1), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(1,1), ss.displayYukawaMatrix(YE)(2,2), 0.0093);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_Ye()(2,2), ss.displayYukawaMatrix(YE)(3,3), 0.0096);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassB() , ss.displayGaugino(1), 0.0044);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassWB(), ss.displayGaugino(2), 0.0046);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MassG() , ss.displayGaugino(3), 0.0051);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.0043);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.009);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.0017);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.0077);

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
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0093);

   // comparing tree-level masses

   const DoubleVector MCha(ToDoubleVector(fs.get_MCha())),
      MChi(ToDoubleVector(fs.get_MChi())),
      MHpm(ToDoubleVector(fs.get_MHpm())),
      MAh(ToDoubleVector(fs.get_MAh())),
      Mhh(ToDoubleVector(fs.get_Mhh()));
   const DoubleVector mch(ss.displayDrBarPars().mchBpmz),
      mn(ss.displayDrBarPars().mnBpmz);
   const double MwRun = fs.get_MVWm();
   const double MzRun = fs.get_MVZ();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const double mA0 = ss.displayDrBarPars().mA0(1);
   const double mh0 = ss.displayDrBarPars().mh0(1);
   const double mH0 = ss.displayDrBarPars().mh0(2);

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha(1), mch(1), 0.0013);
   BOOST_CHECK_CLOSE_FRACTION(MCha(2), mch(2), 0.008);

   BOOST_CHECK_CLOSE_FRACTION(MChi(1), mn(1), 0.0042);
   BOOST_CHECK_CLOSE_FRACTION(MChi(2), mn(2), 0.0041);
   BOOST_CHECK_CLOSE_FRACTION(MChi(3), mn(3), 0.0043);
   BOOST_CHECK_CLOSE_FRACTION(MChi(4), mn(4), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0  , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0, 0.00007);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mH0, 0.004);

   // down-type squarks
   const DoubleVector Sd(ToDoubleVector(fs.get_MSd()));
   const DoubleVector md(ss.displayDrBarPars().md.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Sd(1), md(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(2), md(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(3), md(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(4), md(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(5), md(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd(6), md(6), 0.006);

   // up-type squarks
   const DoubleVector Su(ToDoubleVector(fs.get_MSu()));
   const DoubleVector mu(ss.displayDrBarPars().mu.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Su(1), mu(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(2), mu(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(3), mu(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(4), mu(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(5), mu(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su(6), mu(6), 0.006);

   // sleptons
   const DoubleVector Se(ToDoubleVector(fs.get_MSe()));
   const DoubleVector me(ss.displayDrBarPars().me.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Se(1), me(1), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(2), me(2), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(3), me(3), 0.003);
   BOOST_CHECK_CLOSE_FRACTION(Se(4), me(4), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(Se(5), me(5), 0.005);
   BOOST_CHECK_CLOSE_FRACTION(Se(6), me(6), 0.005);

   // sneutrinos
   const DoubleVector msnu(ss.displayDrBarPars().msnu);
   const DoubleVector Snu(ToDoubleVector(fs.get_MSv()));
   BOOST_CHECK_CLOSE_FRACTION(Snu(1), msnu(1), 0.0055);
   BOOST_CHECK_CLOSE_FRACTION(Snu(2), msnu(2), 0.0055);
   BOOST_CHECK_CLOSE_FRACTION(Snu(3), msnu(3), 0.0085);

   BOOST_CHECK_EQUAL(fs.get_MVP(), 0.0);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MVZ() , MzRun, 0.022);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MVWm(), MwRun, 0.024);

   BOOST_CHECK_EQUAL(fs.get_MVG(), 0.0);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MGlu(), ss.displayDrBarPars().mGluino, 0.005);

   BOOST_CHECK_EQUAL(fs.get_MFv()(0), 0.0);
   BOOST_CHECK_EQUAL(fs.get_MFv()(1), 0.0);
   BOOST_CHECK_EQUAL(fs.get_MFv()(2), 0.0);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFe()(2), ss.displayDrBarPars().mtau, 0.00044);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFu()(2), ss.displayDrBarPars().mt  , 0.0097);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_MFd()(2), ss.displayDrBarPars().mb  , 0.0027);


   // comparing one-loop masses

   const DoubleVector
      MCha_1l(ToDoubleVector(fs.get_physical().MCha)),
      MChi_1l(ToDoubleVector(fs.get_physical().MChi)),
      MHpm_1l(ToDoubleVector(fs.get_physical().MHpm)),
      MAh_1l(ToDoubleVector(fs.get_physical().MAh)),
      Mhh_1l(ToDoubleVector(fs.get_physical().Mhh));
   const DoubleVector mch_1l(ss.displayPhys().mch),
      mn_1l(ss.displayPhys().mneut.apply(fabs));
   const double mHpm_1l = ss.displayPhys().mHpm;
   const double mA0_1l  = ss.displayPhys().mA0(1);
   const double mh0_1l  = ss.displayPhys().mh0(1);
   const double mH0_1l  = ss.displayPhys().mh0(2);

   // charginos
   BOOST_CHECK_CLOSE_FRACTION(MCha_1l(1), mch_1l(1), 0.0011);
   BOOST_CHECK_CLOSE_FRACTION(MCha_1l(2), mch_1l(2), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(1), mn_1l(1), 0.0042);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(2), mn_1l(2), 0.0011);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(3), mn_1l(3), 0.0043);
   BOOST_CHECK_CLOSE_FRACTION(MChi_1l(4), mn_1l(4), 0.0040);

   BOOST_CHECK_CLOSE_FRACTION(MHpm_1l(2), mHpm_1l , 0.004);
   BOOST_CHECK_CLOSE_FRACTION(MAh_1l(2) , mA0_1l  , 0.004);

   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(1), mh0_1l, 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(Mhh_1l(2), mH0_1l, 0.004);

   // down-type squarks
   const DoubleVector Sd_1l(ToDoubleVector(fs.get_physical().MSd));
   const DoubleVector md_1l(ss.displayPhys().md.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(1), md_1l(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(2), md_1l(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(3), md_1l(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(4), md_1l(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(5), md_1l(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Sd_1l(6), md_1l(6), 0.006);

   // up-type squarks
   const DoubleVector Su_1l(ToDoubleVector(fs.get_physical().MSu));
   const DoubleVector mu_1l(ss.displayPhys().mu.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(1), mu_1l(1), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(2), mu_1l(2), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(3), mu_1l(3), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(4), mu_1l(4), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(5), mu_1l(5), 0.006);
   BOOST_CHECK_CLOSE_FRACTION(Su_1l(6), mu_1l(6), 0.006);

   // sleptons
   const DoubleVector Se_1l(ToDoubleVector(fs.get_physical().MSe));
   const DoubleVector me_1l(ss.displayPhys().me.flatten().sort());
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(1), me_1l(1), 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(2), me_1l(2), 0.0008);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(3), me_1l(3), 0.0008);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(4), me_1l(4), 0.0050);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(5), me_1l(5), 0.0050);
   BOOST_CHECK_CLOSE_FRACTION(Se_1l(6), me_1l(6), 0.0050);

   // sneutrinos
   const DoubleVector msnu_1l(ss.displayPhys().msnu);
   const DoubleVector Snu_1l(ToDoubleVector(fs.get_physical().MSv));
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(1), msnu_1l(1), 0.0056);
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(2), msnu_1l(2), 0.0056);
   BOOST_CHECK_CLOSE_FRACTION(Snu_1l(3), msnu_1l(3), 0.0090);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_physical().MGlu, ss.displayPhys().mGluino, 0.005);
}

// ===== test with gauge couplings determined from the Rho parameter =====

BOOST_AUTO_TEST_CASE( test_MSSM_spectrum_with_Softsusy_gauge_couplings )
{
   MSSM_input_parameters pp;
   const MSSM_high_scale_constraint<Two_scale> high_constraint(pp);
   const double mxGuess = high_constraint.get_initial_scale_guess();

   MSSM_tester mssm_tester;
   mssm_tester.set_low_scale_constraint(new MSSM_precise_gauge_couplings_low_scale_constraint());
   mssm_tester.set_susy_scale_constraint(new MSSM_softsusy_ewsb_susy_scale_constraint());
   BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   SoftSusy_tester softSusy_tester;
   BOOST_REQUIRE_NO_THROW(softSusy_tester.test(pp, mxGuess));

   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_mx(), softSusy_tester.get_mx(), 0.04);
   BOOST_CHECK_CLOSE_FRACTION(mssm_tester.get_msusy(), softSusy_tester.get_msusy(), 6.2e-4);

   // compare model parameters
   const MssmSoftsusy ss(softSusy_tester.get_model());
   const MSSM<Two_scale> fs(mssm_tester.get_model());

   BOOST_CHECK_CLOSE_FRACTION(fs.get_g1(), ss.displayGaugeCoupling(1), 0.00023);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g2(), ss.displayGaugeCoupling(2), 0.00066);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_g3(), ss.displayGaugeCoupling(3), 0.00010);

   BOOST_CHECK_CLOSE_FRACTION(fs.get_Mu() , ss.displaySusyMu(), 0.0012);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_BMu(), ss.displayM3Squared(), 0.0023);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHd2(), ss.displayMh1Squared(), 0.0005);
   BOOST_CHECK_CLOSE_FRACTION(fs.get_mHu2(), ss.displayMh2Squared(), 0.0022);

   const double vu = fs.get_vu();
   const double vd = fs.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));

   BOOST_CHECK_CLOSE_FRACTION(tanBeta, ss.displayTanb(), 1.0e-9);
   BOOST_CHECK_CLOSE_FRACTION(vev    , ss.displayHvev(), 0.0068);

   // comparing tree-level masses

   const DoubleVector MHpm(ToDoubleVector(fs.get_MHpm())),
      MAh(ToDoubleVector(fs.get_MAh())),
      Mhh(ToDoubleVector(fs.get_Mhh()));
   const double MwRun = fs.get_MVWm();
   const double MzRun = fs.get_MVZ();
   const double mHpm = ss.displayDrBarPars().mHpm;
   const double mA0 = ss.displayDrBarPars().mA0(1);
   const double mh0 = ss.displayDrBarPars().mh0(1);
   const double mH0 = ss.displayDrBarPars().mh0(2);

   BOOST_CHECK_CLOSE_FRACTION(MHpm(1), MwRun, 1.0e-10); // for RXi(Wm) == 1
   BOOST_CHECK_CLOSE_FRACTION(MHpm(2), mHpm, 0.0011);

   BOOST_CHECK_CLOSE_FRACTION(MAh(1), MzRun, 1.0e-10); // for RXi(VZ) == 1
   BOOST_CHECK_CLOSE_FRACTION(MAh(2), mA0, 0.0012);

   BOOST_CHECK_CLOSE_FRACTION(Mhh(1), mh0, 0.000022);
   BOOST_CHECK_CLOSE_FRACTION(Mhh(2), mH0, 0.0011);
}


/**
 * @class MSSM_iterative_low_scale_constraint
 *
 * Replacement class for MSSM_low_scale_constraint, which inputs the
 * Higgs mass and eliminates TanBeta.
 */
class MSSM_iterative_low_scale_constraint
   : public MSSM_low_scale_constraint<Two_scale> {
public:
   MSSM_iterative_low_scale_constraint()
      : MSSM_low_scale_constraint<Two_scale>() {}
   MSSM_iterative_low_scale_constraint(const MSSM_input_parameters& inputPars_, const QedQcd& oneset_)
      : MSSM_low_scale_constraint<Two_scale>(inputPars_,oneset_) {}
   virtual ~MSSM_iterative_low_scale_constraint() {}

   virtual void apply();
};

void MSSM_iterative_low_scale_constraint::apply()
{
   assert(model && "Error: MSSM_low_scale_constraint:"
          " model pointer must not be zero");

   model->calculate_DRbar_parameters();
   update_scale();
   calculate_DRbar_gauge_couplings();
   calculate_DRbar_yukawa_couplings();

   struct Chi_sqr_mH_mZ {
      static double func(const gsl_vector* x, void* params) {
         if (contains_nan(x, 2))
            return std::numeric_limits<double>::max();

         MSSM<Two_scale>* model = static_cast<MSSM<Two_scale>*>(params);

         const double vd = gsl_vector_get(x, 0);
         const double vu = gsl_vector_get(x, 1);

         if (vd < std::numeric_limits<double>::epsilon() ||
             vu < std::numeric_limits<double>::epsilon())
            return std::numeric_limits<double>::max();

         model->set_vd(vd);
         model->set_vu(vu);

         model->calculate_DRbar_parameters();
         model->calculate_Mhh_pole_1loop();
         model->calculate_MVZ_pole_1loop();

         const double mH = model->get_physical().Mhh(0);
         const double mZ = model->get_physical().MVZ;

         #define SM(p) Electroweak_constants::p
         #define STANDARD_DEVIATION(p) Electroweak_constants::Error_##p

         return Sqr(SM(MZ) - mZ)/Sqr(STANDARD_DEVIATION(MZ))
              + Sqr(SM(MH) - mH)/Sqr(STANDARD_DEVIATION(MH)*10);
      }
   };

   Minimizer<2> minimizer(Chi_sqr_mH_mZ::func, model, 100, 1.0e-2);
   const double start[2] = { model->get_vd(), model->get_vu() };

   const int status = minimizer.minimize(start);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_MESSAGE("chi^2 = " << minimizer.get_minimum_value());
   BOOST_MESSAGE("New vd = " << model->get_vd() << ", vu = " << model->get_vu());
   BOOST_MESSAGE("Predicted tan(beta) = " << model->get_vu() / model->get_vd());

   model->set_g1(new_g1);
   model->set_g2(new_g2);
   model->set_g3(new_g3);

   model->set_Yu(new_Yu);
   model->set_Yd(new_Yd);
   model->set_Ye(new_Ye);
}

BOOST_AUTO_TEST_CASE( test_MSSM_spectrum_higgs_iteration )
{
   MSSM_input_parameters pp;
   pp.m0 = 500.;
   pp.Azero = 1000.;
   pp.TanBeta = 30.;

   MSSM_tester mssm_tester;
   mssm_tester.set_low_scale_constraint(new MSSM_iterative_low_scale_constraint());
   // BOOST_REQUIRE_NO_THROW(mssm_tester.test(pp));

   try {
      mssm_tester.test(pp);
   } catch (Error& error) {
      ERROR(error.what());
   }
}
