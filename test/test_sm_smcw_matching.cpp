#include "two_scale_solver.hpp"
#include "sm_two_scale.hpp"
#include "smcw_two_scale.hpp"
#include "linalg.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_smcw_matching

#include <boost/test/unit_test.hpp>

class Trivial_SM_SMCW_matching_condition: public Two_scale_matching {
public:
   virtual ~Trivial_SM_SMCW_matching_condition() {}
   virtual DoubleVector calcHighFromLowScaleParameters(const DoubleVector& v) const {
      DoubleVector pars(v);
      assert(pars.displayStart() == 1 &&
             pars.displayEnd() == StandardModel<Two_scale>::numStandardModelPars);
      pars.setEnd(StandardModelCW<Two_scale>::numStandardModelCWPars);
      return pars;
   }
   virtual DoubleVector calcLowFromHighScaleParameters(const DoubleVector& v) const {
      DoubleVector pars(v);
      assert(pars.displayStart() == 1 &&
             pars.displayEnd() == StandardModelCW<Two_scale>::numStandardModelCWPars);
      pars.setEnd(StandardModel<Two_scale>::numStandardModelPars);
      return pars;
   }
};

BOOST_AUTO_TEST_CASE( test_trival_matching )
{
   const double vev = 246;
   const double root2 = sqrt(2.0);
   const double mtoprun = 165;
   const double mbrun = 2.9;
   const double mtau = 1.77699;
   const double yt = mtoprun * root2 / vev;
   const double yb = mbrun * root2 / vev;
   const double ytau = mtau * root2 / vev;

   const double MZ = 91.1876;
   const double aem = 1.0 / 127.918; // at MZ
   const double sinthWsq = 0.23122;
   const double alpha1 = 5.0 * aem / (3.0 * (1.0 - sinthWsq));
   const double alpha2 = aem / sinthWsq;
   const double alpha3 = 0.1187; // at MZ

   StandardModel<Two_scale>* sm = new StandardModel<Two_scale>();
   sm->setMu(MZ);
   sm->setYukawaElement(StandardModel<Two_scale>::YU, 3, 3, yt);
   sm->setYukawaElement(StandardModel<Two_scale>::YD, 3, 3, yb);
   sm->setYukawaElement(StandardModel<Two_scale>::YE, 3, 3, ytau);
   sm->setGaugeCoupling(1, sqrt(4 * PI * alpha1));
   sm->setGaugeCoupling(2, sqrt(4 * PI * alpha2));
   sm->setGaugeCoupling(3, sqrt(4 * PI * alpha3));

   StandardModelCW<Two_scale>* smcw = new StandardModelCW<Two_scale>();

   // this trivial matching condition simply forwards the parameters
   // of one model to the other
   Trivial_SM_SMCW_matching_condition mc;

   Two_scale_solver solver;
   solver.add_model(sm);
   solver.add_model(smcw);
   solver.add_matching_condition(&mc);

   try {
      solver.solve();
   } catch (Two_scale_solver::Error& e) {
      BOOST_ERROR(e.what());
   }

   DoubleVector smPars(sm->getParameters());
   DoubleVector smcwPars(smcw->getParameters());

   // check that g4 and lambda are 0.0
   BOOST_CHECK_EQUAL(smcwPars(StandardModelCW<Two_scale>::numStandardModelCWPars - 1), 0.0);
   BOOST_CHECK_EQUAL(smcwPars(StandardModelCW<Two_scale>::numStandardModelCWPars), 0.0);

   // check that the SM parameters are the same in both models
   smcwPars.setEnd(StandardModel<Two_scale>::numStandardModelPars);
   BOOST_CHECK_EQUAL(smcwPars, smPars);

   delete smcw;
   delete sm;
}
