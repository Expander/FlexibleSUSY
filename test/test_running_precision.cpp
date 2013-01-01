
#include "two_scale_solver.hpp"
#include "two_scale_model.hpp"
#include "two_scale_constraint.hpp"
#include "two_scale_convergence_tester.hpp"
#include "two_scale_running_precision.hpp"

#include <cmath>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_running_precision

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( test_constant_running_precision )
{
   for (int i = 0; i < 5; ++i) {
      const double default_precison = static_cast<double>(i);
      Two_scale_constant_precision const_prec(default_precison);

      // test that all iterations yield the same precison
      for (int iteration = 0; iteration < 10; ++iteration) {
         double precision = const_prec.get_precision(iteration);
         BOOST_CHECK_EQUAL(precision, default_precison);
      }
   }
}

BOOST_AUTO_TEST_CASE( test_increasing_running_precision )
{
   for (int i = 1; i < 10; ++i) {
      const double increasing_factor = static_cast<double>(i);
      const double min_precision = 1.0e-12;
      Two_scale_increasing_precision incr_prec(increasing_factor, min_precision);

      // test that all iterations yield expected precision
      for (int iteration = 0; iteration < 20; ++iteration) {
         double precision = incr_prec.get_precision(iteration);
         double expected_precision = std::pow(1./increasing_factor, iteration + 1);
         if (expected_precision < min_precision)
            expected_precision = min_precision;
         BOOST_CHECK_CLOSE(precision, expected_precision, 1.0e-5);
      }
   }
}

class Static_convergence_tester : public Convergence_tester<Two_scale> {
public:
   Static_convergence_tester(unsigned int max_iterations_)
      : iteration(0), maximum_iterations(max_iterations_) {}
   virtual ~Static_convergence_tester() {}
   virtual bool accuracy_goal_reached() {
      return false;
   }
   virtual unsigned int max_iterations() const {
      return maximum_iterations;
   }
private:
   unsigned int iteration, maximum_iterations;
};

class Static_model: public Two_scale_model {
public:
   Static_model() {}
   virtual ~Static_model() {}
   virtual void calculate_spectrum() {}
   virtual std::string name() const { return "Static_model"; }
   virtual int run_to(double, double) { return 0; }
};

class Static_constraint : public Constraint<Two_scale> {
public:
   Static_constraint(double scale_)
      : scale(scale_) {}
   virtual ~Static_constraint() {}
   virtual void apply() {}
   virtual double get_scale() const { return scale; }
private:
   double scale;
};

class Test_increasing_precision : public Two_scale_running_precision {
public:
   Test_increasing_precision() : call(0) {}
   virtual ~Test_increasing_precision() {}
   virtual double get_precision(unsigned it) {
      // Check that this function is called 3 times with the same
      // argument, before the argument is increased, i.e. it is called
      // with 0 0 0  1 1 1  2 2 2  3 3 3 ...
      BOOST_CHECK_EQUAL(call / 3, it);
      ++call;
      return 1.0e-3;
   }
private:
   unsigned call; ///< stores number of calls of get_precision()
};

BOOST_AUTO_TEST_CASE( test_increasing_iteration_number )
{
   // This checks that during each iteration (within
   // RGFlow<Two_scale>) the iteration number is increased and passed
   // to Two_scale_running_precision.
   //
   // We use a model with 2 constraints.  This means that
   // Two_scale_model::run_to() is called 3 times (and so is
   // Two_scale_running_precision::get_precision()): 2 times when we
   // run up, and 1 time to run down.
   //
   // The test in Test_increasing_precision ensures that
   // get_precision() is called exactly 3 times with the same
   // argument.

   Static_model model;

   Static_constraint c1(100), c2(200);
   std::vector<Constraint<Two_scale>*> constraints;
   constraints.push_back(&c1);
   constraints.push_back(&c2);

   Static_convergence_tester ccc(10);
   Test_increasing_precision incr_prec;

   RGFlow<Two_scale> solver;
   solver.add_model(&model, constraints);
   solver.set_convergence_tester(&ccc);
   solver.set_running_precision(&incr_prec);

   BOOST_CHECK_THROW(solver.solve(), RGFlow<Two_scale>::NoConvergenceError);
}
