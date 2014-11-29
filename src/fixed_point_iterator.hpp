// ====================================================================
// Class to do fixed point iteration. Uses std::vector instead of
// gsl_vector at the moment, as the syntax is slightly nicer.
//
// TODO:
//   - implement check for no progress towards solution
// ====================================================================

#ifndef FIXED_POINT_ITERATOR_H
#define FIXED_POINT_ITERATOR_H

#include <iostream>
#include <cassert>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>

#include "wrappers.hpp"
#include "error.hpp"

namespace flexiblesusy {

/**
 * @class Fixed_point_iterator
 * @brief Does fixed point iteration
 * @author Dylan Harries, Alexander Voigt
 * @tparam dimension dimension of function
 * @tparam Compare_relative function for relative comparison
 *    of subsequent iteration steps
 *
 * The user has to provide the function (of which a fixed point
 * should be found) of the type \a Function_t. This function gets as
 * arguments a GSL vector of length \a dimension, a pointer to the
 * parameters (of type \a void*) and a GSL vector where the next
 * point must be stored.
 */
template <std::size_t dimension, double (*Compare_relative)(double,double) = MaxRelDiff>
class Fixed_point_iterator {
public:
   typedef int (*Function_t)(const gsl_vector*, void*, gsl_vector*);

   Fixed_point_iterator();
   Fixed_point_iterator(Function_t, void*, std::size_t, double, bool = false);
   Fixed_point_iterator(const Fixed_point_iterator&);
   ~Fixed_point_iterator();

   double get_fixed_point(std::size_t) const;
   void set_function(Function_t f) { function = f; }
   void set_parameters(void* m) { parameters = m; }
   void set_precision(double p) { precision = p; }
   void set_max_iterations(std::size_t n) { max_iterations = n; }
   void set_test_absolute_errors(bool ae) { test_on_absolute = ae; }
   int find_fixed_point(const double[dimension]);

private:
   std::size_t max_iterations;       ///< maximum number of iterations
   double precision;                 ///< precision goal
   bool test_on_absolute;            ///< use absolute convergence criterion
   gsl_vector* xn;                   ///< current iteration point
   gsl_vector* fixed_point;          ///< vector of fixed point estimate
   void* parameters;                 ///< pointer to parameters
   Function_t function;              ///< function defining fixed point

   int fixed_point_iterator_iterate();
   int fixed_point_iterator_test_relative() const;
   int fixed_point_iterator_test_absolute() const;
   void print_state(std::size_t) const;
};

/**
 * Default constructor
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
Fixed_point_iterator<dimension,Compare_relative>::Fixed_point_iterator()
   : max_iterations(100)
   , precision(1.0e-2)
   , test_on_absolute(false)
   , parameters(NULL)
   , function(NULL)
{
   xn = gsl_vector_alloc(dimension);
   fixed_point = gsl_vector_alloc(dimension);

   if (!xn || !fixed_point)
      throw OutOfMemoryError("GSL vector allocation failed in"
                             " Fixed_point_iterator()");
}

/**
 * Constructor
 *
 * @param function_ pointer to the function to find fixed point for
 * @param parameters_ pointer to the parameters (for example the model)
 * @param max_iterations_ maximum number of iterations
 * @param precision_ precision goal
 * @param absolute_ use absolute convergence test
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
Fixed_point_iterator<dimension,Compare_relative>::Fixed_point_iterator(
   Function_t function_,
   void* parameters_,
   std::size_t max_iterations_,
   double precision_, bool absolute_
)
   : max_iterations(max_iterations_)
   , precision(precision_)
   , test_on_absolute(absolute_)
   , parameters(parameters_)
   , function(function_)
{
   xn = gsl_vector_alloc(dimension);
   fixed_point = gsl_vector_alloc(dimension);

   if (!xn || !fixed_point)
      throw OutOfMemoryError("GSL vector allocation failed in"
                             " Fixed_point_iterator(Function_t, void*,"
                             " size_t, double, bool)");
}

template <std::size_t dimension, double (*Compare_relative)(double,double)>
Fixed_point_iterator<dimension,Compare_relative>::Fixed_point_iterator(
   const Fixed_point_iterator& other
)
   : max_iterations(other.max_iterations)
   , precision(other.precision)
   , test_on_absolute(other.test_on_absolute)
   , parameters(other.parameters)
   , function(other.function)
{
   xn = gsl_vector_alloc(dimension);
   gsl_vector_memcpy(xn, other.xn);

   fixed_point = gsl_vector_alloc(dimension);
   gsl_vector_memcpy(fixed_point, other.fixed_point);
}

template <std::size_t dimension, double (*Compare_relative)(double,double)>
Fixed_point_iterator<dimension,Compare_relative>::~Fixed_point_iterator()
{
   gsl_vector_free(xn);
   gsl_vector_free(fixed_point);
}

/**
 * Start the iteration
 *
 * @param start starting point
 *
 * @return GSL error code (GSL_SUCCESS if fixed point found)
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
int Fixed_point_iterator<dimension,Compare_relative>::find_fixed_point(
   const double start[dimension]
)
{
   assert(function && "Fixed_point_iterator<dimension,Compare_relative>"
          "::find_fixed_point: function pointer must not be zero!");

   int status;
   std::size_t iter = 0;

#ifndef ENABLE_DEBUG
   gsl_set_error_handler_off();
#endif

   for (std::size_t i = 0; i < dimension; ++i) {
      gsl_vector_set(xn, i, start[i]);
      gsl_vector_set(fixed_point, i, start[i]);
   }

#ifdef ENABLE_VERBOSE
   print_state(iter);
#endif

   do {
      iter++;
      status = fixed_point_iterator_iterate();

#ifdef ENABLE_VERBOSE
      print_state(iter);
#endif

      if (status)   // check if iterator has problems
         break;

      status = (test_on_absolute ?
                fixed_point_iterator_test_absolute()
                : fixed_point_iterator_test_relative());

   } while (status == GSL_CONTINUE && iter < max_iterations);

#ifdef ENABLE_VERBOSE
   std::cout << "\tFixed_point_iterator status = "
             << gsl_strerror(status) << '\n';
#endif

   return status;
}

/**
 * Perform a single step of the fixed point iteration
 *
 * @return GSL error code
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
int Fixed_point_iterator<dimension,Compare_relative>::fixed_point_iterator_iterate()
{
   gsl_vector_memcpy(xn, fixed_point);

   int status = function(xn, parameters, fixed_point);

   if (status != GSL_SUCCESS)
      return GSL_EBADFUNC;

   // For safety, include a check for nans or infs here (which
   // should be sufficient for now)
   for (std::size_t i = 0; i < dimension; ++i) {
      if (!gsl_finite(gsl_vector_get(fixed_point, i)))
         GSL_ERROR("update point is not finite", GSL_EBADFUNC);
   }

   return GSL_SUCCESS;
}

/**
 * Test whether the relative difference is less than the set
 * precision. The relative difference test used here is carried out by
 * applying \a Compare_relative to each element of the vector.
 *
 * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
int Fixed_point_iterator<dimension,Compare_relative>::fixed_point_iterator_test_relative() const
{
   double rel_diff = 0.;

   if (precision < 0.)
      GSL_ERROR("relative tolerance is negative", GSL_EBADTOL);

   for (std::size_t i = 0; i < dimension; ++i) {
      rel_diff = Compare_relative(gsl_vector_get(xn, i),
                                  gsl_vector_get(fixed_point, i));

      if (rel_diff > precision)
         return GSL_CONTINUE;
   }

   return GSL_SUCCESS;
}

/**
 * Test whether the absolute value of the residual, defined by
 * |x_{n+1}-x_n| = \f$\sqrt{\sum_i (x_{n+1}(i) - x_n(i))^2}\f$,
 * is less than the set precision.
 *
 * @return GSL error code (GSL_SUCCESS or GSL_CONTINUE)
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
int Fixed_point_iterator<dimension,Compare_relative>::fixed_point_iterator_test_absolute() const
{
   double residual = 0;

   if (precision < 0.)
      GSL_ERROR("absolute tolerance is negative", GSL_EBADTOL);

   for (std::size_t i = 0; i < dimension; ++i) {
      residual += Sqr(gsl_vector_get(fixed_point, i) -
                      gsl_vector_get(xn, i));
   }

   residual = Sqrt(residual);

   return (residual < precision ? GSL_SUCCESS : GSL_CONTINUE);
}

/**
 * Print state of the fixed point iterator
 *
 * @param iteration iteration number
 */
template <std::size_t dimension, double (*Compare_relative)(double,double)>
void Fixed_point_iterator<dimension,Compare_relative>::print_state(std::size_t iteration) const
{
   std::cout << "\tIteration n = " << iteration << ": x_{n} =";
   for (std::size_t i = 0; i < dimension; ++i) {
      std::cout << " " << gsl_vector_get(xn, i);
   }
   std::cout << ", x_{n+1} =";
   for (std::size_t i = 0; i < dimension; ++i) {
      std::cout << " " << gsl_vector_get(fixed_point, i);
   }
   std::cout << '\n';
}

template <std::size_t dimension, double (*Compare_relative)(double,double)>
double Fixed_point_iterator<dimension,Compare_relative>::get_fixed_point(std::size_t i) const
{
   assert(i < dimension && "Fixed_point_iterator<>::get_fixed_point: index out"
          " of bounds");
   return gsl_vector_get(fixed_point, i);
}

} // namespace flexiblesusy

#endif
