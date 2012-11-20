#include "two_scale_solver.hpp"
#include "sm_two_scale.hpp"
#include "sm_two_scale_experimental_constraint.hpp"
#include "smcw_two_scale.hpp"
#include "smcw_two_scale_gut_constraint.hpp"
#include "linalg.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_sm_smcw_matching

#include <boost/test/unit_test.hpp>

#define YU StandardModel<Two_scale>::YU
#define YD StandardModel<Two_scale>::YD
#define YE StandardModel<Two_scale>::YE

class Trivial_SM_SMCW_matching_condition: public Matching<Two_scale> {
public:
   Trivial_SM_SMCW_matching_condition(StandardModel<Two_scale>* sm_,
                                      StandardModelCW<Two_scale>* smcw_)
      : sm(sm_)
      , smcw(smcw_)
      {
         BOOST_REQUIRE(sm != NULL);
         BOOST_REQUIRE(smcw != NULL);
      }
   virtual ~Trivial_SM_SMCW_matching_condition() {}
   virtual void matchLowToHighScaleModel() const {
      smcw->setYukawaMatrix(YU, sm->displayYukawaMatrix(YU));
      smcw->setYukawaMatrix(YD, sm->displayYukawaMatrix(YD));
      smcw->setYukawaMatrix(YE, sm->displayYukawaMatrix(YE));
      for (int i = 1; i <= 3; ++i)
         smcw->setGaugeCoupling(i, sm->displayGaugeCoupling(i));
      smcw->setScale(sm->getScale());
   }
   virtual void matchHighToLowScaleModel() const {
      sm->setYukawaMatrix(YU, smcw->displayYukawaMatrix(YU));
      sm->setYukawaMatrix(YD, smcw->displayYukawaMatrix(YD));
      sm->setYukawaMatrix(YE, smcw->displayYukawaMatrix(YE));
      for (int i = 1; i <= 3; ++i)
         sm->setGaugeCoupling(i, smcw->displayGaugeCoupling(i));
      sm->setScale(smcw->getScale());
   }
   virtual double get_scale() const {
      return 3000;
   }
private:
   StandardModel<Two_scale>* sm;
   StandardModelCW<Two_scale>* smcw;
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
   sm->setScale(MZ);
   sm->setYukawaElement(YU, 3, 3, yt);
   sm->setYukawaElement(YD, 3, 3, yb);
   sm->setYukawaElement(YE, 3, 3, ytau);
   sm->setGaugeCoupling(1, sqrt(4 * PI * alpha1));
   sm->setGaugeCoupling(2, sqrt(4 * PI * alpha2));
   sm->setGaugeCoupling(3, sqrt(4 * PI * alpha3));

   StandardModelCW<Two_scale>* smcw = new StandardModelCW<Two_scale>();

   // this trivial matching condition simply forwards the parameters
   // of one model to the other
   Trivial_SM_SMCW_matching_condition mc(sm, smcw);

   RGFlow<Two_scale> solver;
   solver.add_model(sm, &mc);
   solver.add_model(smcw, NULL);

   try {
      solver.solve();
   } catch (RGFlow<Two_scale>::Error& e) {
      BOOST_ERROR(e.what());
   }

   // check that g4 and lambda are 0.0
   BOOST_CHECK_EQUAL(smcw->displayGaugeCoupling(4), 0.0);
   BOOST_CHECK_EQUAL(smcw->displayLambda(), 0.0);

   // check that the SM parameters are the same in both models
   for (int i = 1; i <= 3; ++i)
      BOOST_CHECK_EQUAL(smcw->displayGaugeCoupling(i),
                        sm->displayGaugeCoupling(i));

   BOOST_CHECK_EQUAL(smcw->displayYukawaMatrix(YU), sm->displayYukawaMatrix(YU));
   BOOST_CHECK_EQUAL(smcw->displayYukawaMatrix(YD), sm->displayYukawaMatrix(YD));
   BOOST_CHECK_EQUAL(smcw->displayYukawaMatrix(YE), sm->displayYukawaMatrix(YE));

   delete smcw;
   delete sm;
}

BOOST_AUTO_TEST_CASE( test_sm_smcw_running )
{
   StandardModel<Two_scale> sm;
   sm.setScale(ewConstants::MZ);
   StandardModelExpConstraint sm_ew_constraint(&sm);
   const std::vector<Constraint<Two_scale>*> sm_constraints(1, &sm_ew_constraint);

   StandardModelCW<Two_scale> smcw;
   StandardModelCWGUTConstraint smcw_gut_constraint(&smcw, 1.0e12);
   const std::vector<Constraint<Two_scale>*> smcw_constraints(1, &smcw_gut_constraint);

   Trivial_SM_SMCW_matching_condition mc(&sm, &smcw);

   RGFlow<Two_scale> solver;
   solver.add_model(&sm, &mc, sm_constraints);
   solver.add_model(&smcw, NULL, smcw_constraints);

   try {
      solver.solve();
   } catch (RGFlow<Two_scale>::Error& e) {
      BOOST_ERROR(e.what());
   }

   // check that the SM parameters are the same in both models
   for (int i = 1; i <= 3; ++i)
      BOOST_CHECK_EQUAL(smcw.displayGaugeCoupling(i),
                        sm.displayGaugeCoupling(i));

   BOOST_CHECK_EQUAL(smcw.displayYukawaMatrix(YU), sm.displayYukawaMatrix(YU));
   BOOST_CHECK_EQUAL(smcw.displayYukawaMatrix(YD), sm.displayYukawaMatrix(YD));
   BOOST_CHECK_EQUAL(smcw.displayYukawaMatrix(YE), sm.displayYukawaMatrix(YE));

   // get the GUT scale value
   const double gut_scale = smcw_gut_constraint.estimate_scale();

   // run smcw to the GUT scale and test equality of g1 and g2
   smcw.run_to(gut_scale);
   const double g1_at_mgut = smcw.displayGaugeCoupling(1);
   const double g2_at_mgut = smcw.displayGaugeCoupling(2);
   BOOST_CHECK_CLOSE(g1_at_mgut, g2_at_mgut, 1.0e-8);
}
