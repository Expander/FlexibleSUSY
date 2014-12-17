
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_ewsb_solver

#include <boost/test/unit_test.hpp>

#include "ewsb_solver.hpp"
#include "root_finder.hpp"
#include "minimizer.hpp"
#include "fixed_point_iterator.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

class Parabola {
public:
   static void reset() { number_of_calls = 0; }
   static unsigned get_number_of_calls() { return number_of_calls; }

   /**
    * Function f(x,y) = ((x-5)^2, (y-1)^2) ,
    *
    * @param x touple (x,y)
    * @return f(x,y)
    */
   static int func(const gsl_vector* x, void*, gsl_vector* f) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      gsl_vector_set(f, 0, (y - 5.)*(y - 5.));
      gsl_vector_set(f, 1, (z - 1.)*(z - 1.));
      number_of_calls++;
      return GSL_SUCCESS;
   }

   /**
    * Function f(x,y) = ((x-5)^2, (y-1)^2) ,
    *
    * @param x touple (x,y)
    * @return f_1(x,y) + f_2(x,y)
    */
   static double chi2(const gsl_vector* x, void*) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      const double chi2 = (y - 5.)*(y - 5.) + (z - 1.)*(z - 1.);
      number_of_calls++;
      return chi2;
   }

   /**
    * Update steps of f(x,y) = ((x-5)^2, (y-1)^2) :
    *
    * (x,y) = (-25/(x-10), -1/(y-2))
    *
    * @param x touple (x,y)
    *
    * @return fixed point iteration update steps
    */
   static int update(const double x[2], void*, double f[2]) {
      const double y = x[0];
      const double z = x[1];
      f[0] = -25./(y - 10.);
      f[1] = -1./(z - 2.);
      number_of_calls++;
      return GSL_SUCCESS;
   }

private:
   static unsigned number_of_calls;
};

unsigned Parabola::number_of_calls = 0;

BOOST_AUTO_TEST_CASE( test_parabola_2dim )
{
   const int dimension = 2;
   const double precision = 1.0e-4;
   const int max_iterations = 1000;

   double start[dimension] = { 9, 9 };

   EWSB_solver* solvers[] = {
      new Root_finder<dimension>(Parabola::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_hybrid),
      new Root_finder<dimension>(Parabola::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_hybrids),
      new Root_finder<dimension>(Parabola::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_broyden),
      new Root_finder<dimension>(Parabola::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_dnewton),
      new Minimizer<dimension>(Parabola::chi2, NULL, max_iterations, precision, gsl_multimin_fminimizer_nmsimplex),
      new Minimizer<dimension>(Parabola::chi2, NULL, max_iterations, precision, gsl_multimin_fminimizer_nmsimplex2),
      new Minimizer<dimension>(Parabola::chi2, NULL, max_iterations, precision, gsl_multimin_fminimizer_nmsimplex2rand),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_relative>(Parabola::update, NULL, max_iterations, precision),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_absolute>(Parabola::update, NULL, max_iterations, precision)
   };

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      Parabola::reset();

      const int status = solvers[i]->solve(start);

      BOOST_CHECK(status == EWSB_solver::SUCCESS);

      BOOST_CHECK_CLOSE_FRACTION(solvers[i]->get_solution(0), 5., 0.01);
      BOOST_CHECK_CLOSE_FRACTION(solvers[i]->get_solution(1), 1., 0.01);

      BOOST_MESSAGE("solver " << i << ": "
                    << (status == EWSB_solver::SUCCESS ?
                        "success" : "fail")
                    << " ("
                    << Parabola::get_number_of_calls() << " function calls)");
   }

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      delete solvers[i];
   }
}

/**
 * More realistic example closer to EWSB eqs (including tadpoles)
 */
class Perturbation {
public:
   static void reset() { number_of_calls = 0; }
   static unsigned get_number_of_calls() { return number_of_calls; }

   /**
    * @param x touple (x,y)
    * @return f(x,y)
    */
   static int func(const gsl_vector* x, void*, gsl_vector* f) {
      const double y = gsl_vector_get(x, 0);
      const double z = gsl_vector_get(x, 1);
      gsl_vector_set(f, 0,  y + z + small  *(y*y + z*z) + 1);
      gsl_vector_set(f, 1, -y + z + small*2*(y*y + z*z) + 2);
      number_of_calls++;
      return GSL_SUCCESS;
   }

   /**
    * @param x touple (x,y)
    * @return (f_1(x,y))^2 + (f_2(x,y))^2
    */
   static double chi2(const gsl_vector* x, void*) {
      gsl_vector* f = gsl_vector_alloc(2);
      func(x, NULL, f);
      const double chi2 = Sqr(gsl_vector_get(f,0)) + Sqr(gsl_vector_get(f,1));
      gsl_vector_free(f);
      number_of_calls++;
      return chi2;
   }

   /**
    * Update steps of f(x,y)
    * @param x touple (x,y)
    * @return fixed point iteration update steps
    */
   static int update(const double x[2], void*, double f[2]) {
      const double y = x[0];
      const double z = x[1];
      f[0] =  0.5 + small*0.5*(y*y + z*z);
      f[1] = -1.5 - small*1.5*(y*y + z*z);
      number_of_calls++;
      return GSL_SUCCESS;
   }

private:
   static unsigned number_of_calls;
   static const double small;
};

unsigned Perturbation::number_of_calls = 0;
const double Perturbation::small = 1/(16 * Pi*Pi);

BOOST_AUTO_TEST_CASE( test_perturbation )
{
   const int dimension = 2;
   const double precision = 1.0e-4;
   const int max_iterations = 1000;

   double start[dimension] = { 10, 10 };

   EWSB_solver* solvers[] = {
      new Root_finder<dimension>(Perturbation::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_hybrid),
      new Root_finder<dimension>(Perturbation::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_hybrids),
      new Root_finder<dimension>(Perturbation::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_broyden),
      new Root_finder<dimension>(Perturbation::func, NULL, max_iterations, precision, gsl_multiroot_fsolver_dnewton),
      new Minimizer<dimension>(Perturbation::chi2, NULL, max_iterations, precision, gsl_multimin_fminimizer_nmsimplex),
      new Minimizer<dimension>(Perturbation::chi2, NULL, max_iterations, precision, gsl_multimin_fminimizer_nmsimplex2),
      new Minimizer<dimension>(Perturbation::chi2, NULL, max_iterations, precision, gsl_multimin_fminimizer_nmsimplex2rand),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_relative>(Perturbation::update, NULL, max_iterations, precision),
      new Fixed_point_iterator<dimension,fixed_point_iterator::Convergence_tester_absolute>(Perturbation::update, NULL, max_iterations, precision)
   };

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      Perturbation::reset();

      const int status = solvers[i]->solve(start);

      BOOST_CHECK(status == EWSB_solver::SUCCESS);

      const double solution_1 = solvers[i]->get_solution(0);
      const double solution_2 = solvers[i]->get_solution(1);

      // true solution calculated with Mathematica 7
      BOOST_CHECK_CLOSE_FRACTION(solution_1, 0.508177, 0.01);
      BOOST_CHECK_CLOSE_FRACTION(solution_2, -1.52453, 0.01);

      BOOST_MESSAGE("solver " << i << ": "
                    << (status == EWSB_solver::SUCCESS ?
                        "success" : "fail")
                    << ", solution = (" << solution_1 << ", " << solution_2
                    << ") ("
                    << Perturbation::get_number_of_calls() << " function calls)");
   }

   for (int i = 0; i < sizeof(solvers)/sizeof(solvers[0]); i++) {
      delete solvers[i];
   }
}
