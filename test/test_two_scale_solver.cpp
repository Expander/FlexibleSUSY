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
   virtual int run_to(double) { return 0; }
   virtual void set_parameters(const DoubleVector& v) { parameters = v; }
   virtual DoubleVector get_parameters() const { return parameters; }
private:
   DoubleVector parameters;
};

class Trivial_matching_condition: public Matching<Two_scale> {
public:
   Trivial_matching_condition(Static_model* mLow_, Static_model* mHigh_)
      : mLow(mLow_)
      , mHigh(mHigh_)
      {
         BOOST_REQUIRE(mLow != NULL);
         BOOST_REQUIRE(mHigh != NULL);
      }
   virtual ~Trivial_matching_condition() {}
   virtual void match_low_to_high_scale_model() {
      mHigh->set_parameters(mLow->get_parameters());
   }
   virtual void match_high_to_low_scale_model() {
      mLow->set_parameters(mHigh->get_parameters());
   }
   virtual double get_scale() const {
      return 100.0;
   }
   virtual void update_scale() {}
private:
   Static_model *mLow, *mHigh;
};

class Counting_model: public Two_scale_model {
public:
   Counting_model() : number_of_runs(0) {}
   virtual ~Counting_model() {}
   virtual int run_to(double) { ++number_of_runs; return 0; }
   unsigned get_number_of_runs() const {
      return number_of_runs;
   }
private:
   unsigned number_of_runs;
};

class Counting_matching_condition: public Matching<Two_scale> {
public:
   Counting_matching_condition(double scale_)
      : scale(scale_)
      , number_of_low_to_high_matches(0)
      , number_of_high_to_low_matches(0)
      , number_of_get_scale(0)
      , number_of_update_scale(0)
      {}
   virtual ~Counting_matching_condition() {}
   virtual void match_low_to_high_scale_model() {
      ++number_of_low_to_high_matches;
   }
   virtual void match_high_to_low_scale_model() {
      ++number_of_high_to_low_matches;
   }
   virtual double get_scale() const {
      ++number_of_get_scale;
      return scale;
   }
   virtual void update_scale() {
      ++number_of_update_scale;
   }
   unsigned get_number_of_low_to_high_matches() const {
      return number_of_low_to_high_matches;
   }
   unsigned get_number_of_high_to_low_matches() const {
      return number_of_high_to_low_matches;
   }
   unsigned get_number_of_get_scale() const {
      return number_of_get_scale;
   }
   unsigned get_number_of_update_scale() const {
      return number_of_update_scale;
   }
private:
   double   scale;
   unsigned number_of_low_to_high_matches;
   unsigned number_of_high_to_low_matches;
   mutable unsigned number_of_get_scale;
   unsigned number_of_update_scale;
};

BOOST_AUTO_TEST_CASE( test_unchanged_parameters )
{
   const DoubleVector parameters(10);
   Static_model model(parameters);
   RGFlow<Two_scale> solver;
   solver.add_model(&model);

   try {
      solver.solve();
   } catch (RGFlow<Two_scale>::Error& e) {
      BOOST_ERROR(e.what());
   }

   BOOST_CHECK_EQUAL(model.get_parameters(), parameters);
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

   RGFlow<Two_scale> solver;
   solver.add_model(&model1, &mc);
   solver.add_model(&model2);

   try {
      solver.solve();
   } catch (RGFlow<Two_scale>::Error& e) {
      BOOST_ERROR(e.what());
   }

   // the high scale parameters should be the same as the low scale
   // parameters
   BOOST_CHECK_EQUAL(model1.get_parameters(), parameters);
   BOOST_CHECK_EQUAL(model2.get_parameters(), parameters);
}

BOOST_AUTO_TEST_CASE( test_count_method_calls )
{
   Counting_model model1, model2;
   Counting_matching_condition mc(3000);
   const unsigned number_of_iterations = 0;
   const unsigned number_of_constraints = 0;

   RGFlow<Two_scale> solver;
   solver.set_max_iterations(number_of_iterations);
   solver.add_model(&model1, &mc);
   solver.add_model(&model2);

   try {
      solver.solve();
   } catch (RGFlow<Two_scale>::Error& e) {
      BOOST_ERROR(e.what());
   }

   BOOST_CHECK_EQUAL(model1.get_number_of_runs(),
                     number_of_constraints + 1);
   BOOST_CHECK_EQUAL(model2.get_number_of_runs(),
                     number_of_constraints + 1);
   BOOST_CHECK_EQUAL(mc.get_number_of_low_to_high_matches(),
                     number_of_iterations + 1);
   BOOST_CHECK_EQUAL(mc.get_number_of_high_to_low_matches(),
                     number_of_iterations + 1);
}
