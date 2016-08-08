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

#include "two_scale_solver.hpp"
#include "two_scale_constraint.hpp"
#include "two_scale_convergence_tester.hpp"
#include "two_scale_initial_guesser.hpp"
#include "two_scale_matching.hpp"
#include "two_scale_model.hpp"
#include "two_scale_running_precision.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "functors.hpp"

#include <cmath>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cassert>
#include <sstream>

/**
 * @file two_scale_solver.cpp
 * @brief contains the implementation of the RGFlow<Two_scale> class members
 */

namespace flexiblesusy {

/**
 * Create empty two scale solver.
 * The RG running precision is set to the default value 0.001.
 */
RGFlow<Two_scale>::RGFlow()
   : sliders()
   , iteration(0)
   , convergence_tester(NULL)
   , initial_guesser(NULL)
   , running_precision_calculator(NULL)
   , running_precision(1.0e-3)
   , model_at_this_scale(NULL)
{
}

RGFlow<Two_scale>::~RGFlow()
{
}

/**
 * Adding a model constraint
 *
 * @param c constraint
 * @param m model
 */
void RGFlow<Two_scale>::add(Constraint<Two_scale>* c, Two_scale_model* m)
{
   if (!c) throw SetupError("constraint pointer is NULL");
   if (!m) throw SetupError("model pointer is NULL");
   sliders.push_back(std::make_shared<Constraint_slider>(m, c));
}

/**
 * Adds an upwards matching condition.  This matching condition
 * matches the low-scale model to the high-scale model by calling
 * match_low_to_high_scale_model().
 *
 * @param mc matching condition
 * @param m1 low-scale model
 * @param m2 high-scale model
 */
void RGFlow<Two_scale>::add_upwards(Matching<Two_scale>* mc, Two_scale_model* m1, Two_scale_model* m2)
{
   if (!mc) throw SetupError("matching condition pointer is NULL");
   if (!m1) throw SetupError("model pointer 1 is NULL");
   if (!m2) throw SetupError("model pointer 2 is NULL");
   mc->set_models(m1, m2);
   sliders.push_back(std::make_shared<Matching_up_slider>(m1, m2, mc));
}

/**
 * Adds a downwards matching condition.  This matching condition
 * matches the high-scale model to the low-scale model by calling
 * match_high_to_low_scale_model().
 *
 * @param mc matching condition
 * @param m1 high-scale model
 * @param m2 low-scale model
 */
void RGFlow<Two_scale>::add_downwards(Matching<Two_scale>* mc, Two_scale_model* m1, Two_scale_model* m2)
{
   if (!mc) throw SetupError("matching condition pointer is NULL");
   if (!m1) throw SetupError("model pointer 1 is NULL");
   if (!m2) throw SetupError("model pointer 2 is NULL");
   mc->set_models(m2, m1);
   sliders.push_back(std::make_shared<Matching_down_slider>(m1, m2, mc));
}

/**
 * @brief Solves the boundary value problem.
 *
 * At first the initial_guess() is called.  Afterwards, the function
 * iteratively runs the tower up and down and imposes the boundary
 * conditions.  The iteration stops if either the maximum number of
 * iterations is reached or the precision goal is achieved (defined by
 * the convergence_tester).
 */
void RGFlow<Two_scale>::solve()
{
   check_setup();

   const unsigned max_iterations = get_max_iterations();
   if (sliders.empty() || max_iterations == 0)
      return;

   initial_guess();

   iteration = 0;
   bool accuracy_reached = false;

   while (iteration < max_iterations && !accuracy_reached) {
      update_running_precision();
      clear_problems();
      run_sliders();
      accuracy_reached = accuracy_goal_reached();
      ++iteration;
   }

   if (!accuracy_reached)
      throw NoConvergenceError(max_iterations);

   VERBOSE_MSG("convergence reached after " << iteration << " iterations");
}

/**
 * Sanity checks the models and boundary condtitions.
 */
void RGFlow<Two_scale>::check_setup() const
{
   if (!convergence_tester) {
      throw SetupError("RGFlow<Two_scale>::Error: convergence tester must "
                       "not be NULL");
   }
}

void RGFlow<Two_scale>::clear_problems()
{
   VERBOSE_MSG("> clearing problems ...");

   for (auto& s: sliders)
      s->clear_problems();
}

/**
 * Does the initial guess by calling the guess() method of the initial
 * guesser (if given).
 */
void RGFlow<Two_scale>::initial_guess()
{
   if (initial_guesser)
      initial_guesser->guess();
}

void RGFlow<Two_scale>::run_sliders()
{
   VERBOSE_MSG("> running all models (iteration " << iteration << ") ...");

   for (auto& s: sliders) {
      s->set_precision(get_precision());
      s->slide();
   }

   VERBOSE_MSG("> running slides finished");
}

/**
 * Returns the precision of the RG running.
 *
 * @return RG running precision
 */
double RGFlow<Two_scale>::get_precision()
{
   return running_precision;
}

/**
 * Recalculates the precision of the RG running using the user defined
 * Two_scale_running_precision_calculator class.
 */
void RGFlow<Two_scale>::update_running_precision()
{
   if (running_precision_calculator)
      running_precision = running_precision_calculator->get_precision(iteration);
}

/**
 * Returns the value returned by the accuracy_goal_reached() method of
 * the convergence tester.
 */
bool RGFlow<Two_scale>::accuracy_goal_reached() const
{
   return convergence_tester->accuracy_goal_reached();
}

/**
 * Set the convergence tester to be used during the iteration.
 *
 * @param convergence_tester_ the convergence tester to be used
 */
void RGFlow<Two_scale>::set_convergence_tester(Convergence_tester<Two_scale>* convergence_tester_)
{
   convergence_tester = convergence_tester_;
}

void RGFlow<Two_scale>::set_initial_guesser(Initial_guesser<Two_scale>* ig)
{
   initial_guesser = ig;
}

/**
 * Set RG running precision calculator.
 *
 * @param rp running precision calculator
 */
void RGFlow<Two_scale>::set_running_precision(Two_scale_running_precision* rp)
{
   running_precision_calculator = rp;
}

/**
 * Returns the number of performed iterations
 * @return number of performed iterations
 */
unsigned int RGFlow<Two_scale>::number_of_iterations_done() const
{
   return iteration;
}

/**
 * Returns the maximum number of iterations set in the convergence
 * tester.
 */
unsigned int RGFlow<Two_scale>::get_max_iterations() const
{
   return convergence_tester->max_iterations();
}

/**
 * @brief resets the solver to the initial condition
 *
 * The pointers to the models, matching conditions, convergence
 * tester, initial guesser, and running precision calculator are set
 * to zero.  The runnin precision is set to the default value 0.001.
 */
void RGFlow<Two_scale>::reset()
{
   sliders.clear();

   iteration = 0;
   convergence_tester = NULL;
   initial_guesser = NULL;
   running_precision_calculator = NULL;
   running_precision = 1.0e-3;
   model_at_this_scale = NULL;
}

/* Implementation of sliders */

void RGFlow<Two_scale>::Constraint_slider::clear_problems() {
   model->clear_problems();
}

void RGFlow<Two_scale>::Constraint_slider::slide() {
   VERBOSE_MSG("> \trunning " << model->name() << " to scale " << constraint->get_scale() << " GeV");
   model->run_to(constraint->get_scale());
   VERBOSE_MSG("> \tapplying constraint");
   constraint->apply();
}

void RGFlow<Two_scale>::Constraint_slider::set_precision(double p) {
   model->set_precision(p);
}

void RGFlow<Two_scale>::Matching_up_slider::clear_problems() {
   model1->clear_problems();
   model2->clear_problems();
}

void RGFlow<Two_scale>::Matching_up_slider::slide() {
   VERBOSE_MSG("> \trunning " << model1->name() << " to scale " << matching->get_scale() << " GeV");
   model1->run_to(matching->get_scale());
   VERBOSE_MSG("> \trunning " << model2->name() << " to scale " << matching->get_scale() << " GeV");
   model2->run_to(matching->get_scale());
   VERBOSE_MSG("> \tmatching " << model1->name() << " -> " << model2->name());
   matching->match_low_to_high_scale_model();
}

void RGFlow<Two_scale>::Matching_up_slider::set_precision(double p) {
   model1->set_precision(p);
   model2->set_precision(p);
}

void RGFlow<Two_scale>::Matching_down_slider::clear_problems() {
   model1->clear_problems();
   model2->clear_problems();
}

void RGFlow<Two_scale>::Matching_down_slider::slide() {
   VERBOSE_MSG("> \trunning " << model1->name() << " to scale " << matching->get_scale() << " GeV");
   model1->run_to(matching->get_scale());
   VERBOSE_MSG("> \trunning " << model2->name() << " to scale " << matching->get_scale() << " GeV");
   model2->run_to(matching->get_scale());
   VERBOSE_MSG("> \tmatching " << model1->name() << " -> " << model2->name());
   matching->match_high_to_low_scale_model();
}

void RGFlow<Two_scale>::Matching_down_slider::set_precision(double p) {
   model1->set_precision(p);
   model2->set_precision(p);
}

} // namespace flexiblesusy
