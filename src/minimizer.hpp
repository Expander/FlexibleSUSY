// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef MINIMIZER_H
#define MINIMIZER_H

#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>

#include "logger.hpp"

namespace flexiblesusy {

/**
 * @class Minimizer
 * @brief Function minimizer
 *
 * The user has to provide the function to be minimized of the type
 * Function_t.  This function gets as arguments a GSL vector of lenght
 * `dimension' and a pointer to the model (of type void*).
 */
template <class Model_t, std::size_t dimension>
class Minimizer {
public:
   /// pointer to function to minimize
   typedef double (*Function_t)(const gsl_vector*, void*);

   Minimizer(Model_t*, Function_t, std::size_t, double);
   Minimizer(const Minimizer&);
   ~Minimizer();

   double get_minimum_value() const { return minimum_value; }
   int minimize(const double[dimension]);

private:
   std::size_t max_iterations; ///< maximum number of iterations
   double precision;           ///< precision goal
   double initial_step_size;   ///< initial step size
   double minimum_value;       ///< minimum function value found
   gsl_vector* starting_point; ///< GSL vector of starting point
   gsl_vector* step_size;      ///< GSL vector of initial step size
   Model_t* model;             ///< pointer to model
   Function_t function;        ///< function to minimize

   void print_state(gsl_multimin_fminimizer*, std::size_t) const;
};

/**
 * Constructor
 *
 * @param model_ pointer to the model
 * @param function_ pointer to the function to minimize
 * @param max_iterations_ maximum number of iterations
 * @param precision_ precision goal
 */
template <class Model_t, std::size_t dimension>
Minimizer<Model_t,dimension>::Minimizer(Model_t* model_, Function_t function_,
                                        std::size_t max_iterations_,
                                        double precision_)
   : max_iterations(max_iterations_)
   , precision(precision_)
   , initial_step_size(1.0)
   , minimum_value(0.0)
   , model(model_)
   , function(function_)
{
   starting_point = gsl_vector_alloc(dimension);
   step_size = gsl_vector_alloc(dimension);
}

template <class Model_t, std::size_t dimension>
Minimizer<Model_t,dimension>::Minimizer(const Minimizer& other)
   : max_iterations(other.max_iterations)
   , precision(other.precision)
   , initial_step_size(other.initial_step_size)
   , minimum_value(other.minimum_value)
   , model(other.model)
   , function(other.function)
{
   starting_point = gsl_vector_alloc(dimension);
   step_size = gsl_vector_alloc(dimension);
   // copy vectors
   gsl_vector_memcpy(starting_point, other.starting_point);
   gsl_vector_memcpy(step_size, other.step_size);
}

template <class Model_t, std::size_t dimension>
Minimizer<Model_t,dimension>::~Minimizer()
{
   gsl_vector_free(starting_point);
   gsl_vector_free(step_size);
}

/**
 * Start the minimization
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if minimum found)
 */
template <class Model_t, std::size_t dimension>
int Minimizer<Model_t,dimension>::minimize(const double start[dimension])
{
   assert(model && "Minimizer<dimension>::minimize: model pointer"
          " must not be zero!");
   assert(function && "Minimizer<dimension>::minimize: function pointer"
          " must not be zero!");

   const gsl_multimin_fminimizer_type *type =
      gsl_multimin_fminimizer_nmsimplex2;
   gsl_multimin_fminimizer *minimizer;
   gsl_multimin_function minex_func;

   // Set starting point
   for (std::size_t i = 0; i < dimension; i++)
      gsl_vector_set(starting_point, i, start[i]);

   // Set initial step sizes
   gsl_vector_set_all(step_size, initial_step_size);

   // Initialize method and iterate
   minex_func.n = dimension;
   minex_func.f = function;
   minex_func.params = model;

   minimizer = gsl_multimin_fminimizer_alloc(type, dimension);
   gsl_multimin_fminimizer_set(minimizer, &minex_func, starting_point, step_size);

   size_t iter = 0;
   int status;

   do {
      iter++;
      status = gsl_multimin_fminimizer_iterate(minimizer);

      if (status)
         break;

      const double size = gsl_multimin_fminimizer_size(minimizer);
      status = gsl_multimin_test_size(size, precision);

#ifdef VERBOSE
      print_state(minimizer, iter);
#endif
   } while (status == GSL_CONTINUE && iter < max_iterations);

#ifdef VERBOSE
   printf("\tMinimization status = %s\n", gsl_strerror(status));
#endif

   minimum_value = minimizer->fval;
   gsl_multimin_fminimizer_free(minimizer);

   return status;
}

/**
 * Print state of the minimizer
 *
 * @param minimizer minimizer
 * @param iteration iteration number
 */
template <class Model_t, std::size_t dimension>
void Minimizer<Model_t,dimension>::print_state(gsl_multimin_fminimizer* minimizer,
                                               std::size_t iteration) const
{
   std::cout << "\tIteration " << iteration << ": x =";
   for (std::size_t i = 0; i < dimension; ++i)
      std::cout << " " << gsl_vector_get(minimizer->x, i);
   std::cout << ", f(x) = " << minimizer->fval << '\n';
}

} // namespace flexiblesusy

#endif
