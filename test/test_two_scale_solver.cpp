#include "two_scale_solver.hpp"
#include "linalg.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_scale_solver

#include <boost/test/unit_test.hpp>

class Static_model: public Two_scale_model {
public:
   Static_model() : parameters(1) {}
   Static_model(const DoubleVector& pars) : parameters(pars) {}
   virtual ~Static_model() {}
   virtual void run_up() {}
   virtual void run_down() {}
   virtual void setParameters(const DoubleVector& v) { parameters = v; }
   virtual DoubleVector getParameters() const { return parameters; }
private:
   DoubleVector parameters;
};

class Trivial_matching_condition: public Two_scale_matching {
public:
   Trivial_matching_condition(Static_model* mLow_, Static_model* mHigh_)
      : mLow(mLow_)
      , mHigh(mHigh_)
      {
         BOOST_REQUIRE(mLow != NULL);
         BOOST_REQUIRE(mHigh != NULL);
      }
   virtual ~Trivial_matching_condition() {}
   virtual void matchLowToHighScaleModel() const {
      mHigh->setParameters(mLow->getParameters());
   }
   virtual void matchHighToLowScaleModel() const {
      mLow->setParameters(mHigh->getParameters());
   }
private:
   Static_model *mLow, *mHigh;
};

BOOST_AUTO_TEST_CASE( test_unchanged_parameters )
{
   const DoubleVector parameters(10);
   Static_model model(parameters);
   Two_scale_solver solver;
   solver.add_model(&model);

   try {
      solver.solve();
   } catch (Two_scale_solver::Error& e) {
      BOOST_ERROR(e.what());
   }

   BOOST_CHECK_EQUAL(model.getParameters(), parameters);
}

BOOST_AUTO_TEST_CASE( test_trival_matching )
{
   DoubleVector parameters(10);
   const DoubleVector zeros(10);
   for (int i = 1; i <= 10; ++i)
      parameters(i) = i;

   // init two models:
   // * low-energy model with parameters
   // * high-energy model with zeros
   Static_model model1(parameters), model2(zeros);

   // this trivial matching condition simply forwards the parameters
   // of one model to the other
   Trivial_matching_condition mc(&model1, &model2);

   Two_scale_solver solver;
   solver.add_model(&model1);
   solver.add_model(&model2);
   solver.add_matching_condition(&mc);

   try {
      solver.solve();
   } catch (Two_scale_solver::Error& e) {
      BOOST_ERROR(e.what());
   }

   // the high scale parameters should be the same as the low scale
   // parameters
   BOOST_CHECK_EQUAL(model1.getParameters(), parameters);
   BOOST_CHECK_EQUAL(model2.getParameters(), parameters);
}
