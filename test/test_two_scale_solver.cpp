#include "two_scale_solver.hpp"
#include "two_scale_matching.hpp"
#include "two_scale_model.hpp"
#include "two_scale_constraint.hpp"
#include "two_scale_convergence_tester.hpp"
#include "linalg.h"
#include "error.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_scale_solver

#include <boost/test/unit_test.hpp>

using namespace flexiblesusy;

class Static_model: public Two_scale_model {
public:
   Static_model() : parameters(1) {}
   Static_model(const DoubleVector& pars) : parameters(pars) {}
   virtual ~Static_model() {}
   virtual void calculate_spectrum() {}
   virtual std::string name() const { return "Static_model"; }
   virtual void run_to(double, double) {}
   virtual void set_parameters(const DoubleVector& v) { parameters = v; }
   virtual DoubleVector get_parameters() const { return parameters; }
   virtual void set_precision(double) {}
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
private:
   Static_model *mLow, *mHigh;
};

class Counting_model: public Two_scale_model {
public:
   Counting_model() : number_of_runs(0) {}
   virtual ~Counting_model() {}
   virtual void calculate_spectrum() {}
   virtual void run_to(double, double) { ++number_of_runs; }
   virtual void set_precision(double) {}
   unsigned get_number_of_runs() const {
      return number_of_runs;
   }
private:
   unsigned number_of_runs;
};

class Counting_constraint : public Constraint<Two_scale> {
public:
   Counting_constraint(double scale_)
      : scale(scale_)
      , number_of_apply_calls(0) {}
   virtual ~Counting_constraint() {}
   virtual void apply() { ++number_of_apply_calls; }
   virtual double get_scale() const { return scale; }
   virtual void set_model(Two_scale_model*) {}

   unsigned get_number_of_apply_calls() const {
      return number_of_apply_calls;
   }

private:
   double   scale;
   unsigned number_of_apply_calls;
};

class Counting_matching_condition: public Matching<Two_scale> {
public:
   Counting_matching_condition(double scale_)
      : scale(scale_)
      , number_of_low_to_high_matches(0)
      , number_of_high_to_low_matches(0)
      , number_of_get_scale(0)
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
   unsigned get_number_of_low_to_high_matches() const {
      return number_of_low_to_high_matches;
   }
   unsigned get_number_of_high_to_low_matches() const {
      return number_of_high_to_low_matches;
   }
   unsigned get_number_of_get_scale() const {
      return number_of_get_scale;
   }
private:
   double   scale;
   unsigned number_of_low_to_high_matches;
   unsigned number_of_high_to_low_matches;
   mutable unsigned number_of_get_scale;
};

class Counting_convergence_tester : public Convergence_tester<Two_scale> {
public:
   Counting_convergence_tester(unsigned int max_iterations_)
      : iteration(0), maximum_iterations(max_iterations_) {}
   virtual ~Counting_convergence_tester() {}
   virtual bool accuracy_goal_reached() {
      return false;
   }
   virtual unsigned int max_iterations() const {
      return maximum_iterations;
   }
private:
   unsigned int iteration, maximum_iterations;
};

BOOST_AUTO_TEST_CASE( test_unchanged_parameters )
{
   const DoubleVector parameters(10);
   Static_model model(parameters);
   Counting_convergence_tester ccc(10);

   RGFlow<Two_scale> solver;
   solver.add_model(&model);
   solver.set_convergence_tester(&ccc);

   BOOST_CHECK_THROW(solver.solve(), NoConvergenceError);
   BOOST_CHECK_EQUAL(model.get_parameters(), parameters);
   BOOST_CHECK_EQUAL(solver.get_model(), &model);
}

BOOST_AUTO_TEST_CASE( test_trival_matching )
{
   DoubleVector parameters(10);
   const DoubleVector zeros(10);
   for (int i = 1; i <= 10; ++i)
      parameters(i) = i;

   Counting_convergence_tester ccc(10);

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
   solver.set_convergence_tester(&ccc);

   BOOST_CHECK_THROW(solver.solve(), NoConvergenceError);

   // the high scale parameters should be the same as the low scale
   // parameters
   BOOST_CHECK_EQUAL(model1.get_parameters(), parameters);
   BOOST_CHECK_EQUAL(model2.get_parameters(), parameters);
}

BOOST_AUTO_TEST_CASE( test_count_method_calls )
{
   for (unsigned number_of_iterations = 0;
        number_of_iterations < 10; ++number_of_iterations ) {
      Counting_model model1, model2;
      Counting_constraint model1_c1(1000), model1_c2(2000);
      Counting_constraint model2_c1(4000), model2_c2(5000);

      std::vector<Constraint<Two_scale>*> model1_constraints;
      model1_constraints.push_back(&model1_c1);
      model1_constraints.push_back(&model1_c2);

      std::vector<Constraint<Two_scale>*> model2_constraints;
      model2_constraints.push_back(&model2_c1);
      model2_constraints.push_back(&model2_c2);

      Counting_matching_condition mc(3000);
      Counting_convergence_tester ccc(number_of_iterations);

      RGFlow<Two_scale> solver;
      solver.set_convergence_tester(&ccc);
      solver.add_model(&model1, &mc, model1_constraints);
      solver.add_model(&model2, model2_constraints);

      if (number_of_iterations == 0) {
         try {
            solver.solve();
         } catch (Error& e) {
            BOOST_ERROR(e.what());
         }
      } else {
         // expect NoConvergenceError because the accuracy_goal_reached()
         // function of the convergence tester returns always false
         BOOST_CHECK_THROW(solver.solve(), NoConvergenceError);
      }

      // check that all iterations were done
      BOOST_CHECK_EQUAL(solver.number_of_iterations_done(),
                        number_of_iterations);
      // check how often the matching is appied
      BOOST_CHECK_EQUAL(mc.get_number_of_low_to_high_matches(),
                        number_of_iterations);
      BOOST_CHECK_EQUAL(mc.get_number_of_high_to_low_matches(),
                        number_of_iterations);
      // lowest constraint applied first time and once in each iteration
      BOOST_CHECK_EQUAL(model1_c1.get_number_of_apply_calls(),
                        (number_of_iterations == 0 ? 0 : number_of_iterations + 1));
      // highest constraint applied once in each iteration
      BOOST_CHECK_EQUAL(model2_c2.get_number_of_apply_calls(),
                        number_of_iterations);
      // intermediate constraints applied twice in each iteration
      BOOST_CHECK_EQUAL(model1_c2.get_number_of_apply_calls(),
                        2 * number_of_iterations);
      BOOST_CHECK_EQUAL(model2_c1.get_number_of_apply_calls(),
                        2 * number_of_iterations);
   }
}

BOOST_AUTO_TEST_CASE( test_run_to_with_zero_models )
{
   RGFlow<Two_scale> solver;

   int status = 0;
   try { solver.run_to(1000.); } catch (Error&) { status = 1; }

   BOOST_CHECK_EQUAL(status, 1);
   BOOST_CHECK_EQUAL(solver.get_model(), (void*)NULL);
}

BOOST_AUTO_TEST_CASE( test_run_to_with_one_model )
{
   Static_model model(DoubleVector(10));
   RGFlow<Two_scale> solver;
   solver.add_model(&model);

   int status = 0;
   try { solver.run_to(1000.); } catch (Error&) { status = 1; }

   BOOST_CHECK_EQUAL(status, 0);
   BOOST_CHECK_EQUAL(solver.get_model(), &model);
}

BOOST_AUTO_TEST_CASE( test_run_to_with_two_models )
{
   Static_model model1(DoubleVector(10)), model2(DoubleVector(10));
   Trivial_matching_condition mc(&model1, &model2);
   const double mc_scale = mc.get_scale();

   RGFlow<Two_scale> solver;
   solver.add_model(&model1, &mc);
   solver.add_model(&model2);

   solver.run_to(mc_scale);
   BOOST_CHECK_EQUAL(solver.get_model(), &model1);

   solver.run_to(mc_scale / 2.);
   BOOST_CHECK_EQUAL(solver.get_model(), &model1);

   solver.run_to(mc_scale * 2.);
   BOOST_CHECK_EQUAL(solver.get_model(), &model2);
}
