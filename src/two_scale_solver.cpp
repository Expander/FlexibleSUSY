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
#include "logger.hpp"

#include <cmath>
#include <algorithm>
#include <iterator>

/**
 * Create empty two scale solver.  Sets maximum number of iterations
 * to 10 (default).
 */
RGFlow<Two_scale>::RGFlow()
   : models()
   , max_iterations(10)
   , iteration(0)
   , convergence_tester(NULL)
   , initial_guesser(NULL)
{
}

RGFlow<Two_scale>::~RGFlow()
{
   for (size_t m = 0; m < models.size(); ++m)
      delete models[m];
}

void RGFlow<Two_scale>::solve()
{
   if (models.empty() || max_iterations == 0)
      return;

   check_setup();
   initial_guess();

   unsigned int iter = max_iterations;
   while (iter) {
      iteration = max_iterations - iter;
      run_up();
      run_down();
      if (accuracy_goal_reached())
         break;
      --iter;
   }

   apply_lowest_constaint();

   // save number of iterations that were done
   iteration = max_iterations - iter;

   if (iter == 0 && convergence_tester)
      throw NoConvergenceError(max_iterations);

   VERBOSE_MSG("convergence reached after " << iteration << " iterations");
}

void RGFlow<Two_scale>::check_setup() const
{
   for (size_t m = 0; m < models.size(); ++m) {
      TModel* model = models[m];
      if (!model->model) {
         std::stringstream message;
         message << "RGFlow<Two_scale>::Error: model pointer ["
                 << m << "] is NULL";
         throw SetupError(message.str());
      }

      // check wether last model has a non-zero matching condition
      if (m + 1 == models.size()) {
         if (model->matching_condition)
            WARNING("the matching condition of the " << model->model->name()
                    << " is non-zero but will not be used");
      } else {
         if (model->matching_condition == NULL) {
            std::stringstream message;
            message << "RGFlow<Two_scale>::Error: matching condition "
                    << "of the " << model->model->name() << " to the "
                    << models[m + 1]->model->name() << " is NULL";
            throw SetupError(message.str());
         }
      }
   }
}

void RGFlow<Two_scale>::initial_guess()
{
   if (initial_guesser)
      initial_guesser->guess();
}

void RGFlow<Two_scale>::run_up()
{
   VERBOSE_MSG("> running tower up ...");
   for (size_t m = 0; m < models.size(); ++m) {
      TModel* model = models[m];
      VERBOSE_MSG("> \tselecting model " << model->model->name());
      // apply all constraints
      for (size_t c = 0; c < model->upwards_constraints.size(); ++c) {
         Constraint<Two_scale>* constraint = model->upwards_constraints[c];
         const double scale = constraint->get_scale();
         VERBOSE_MSG("> \t\tselecting constraint " << c << " at scale " << scale);
         VERBOSE_MSG("> \t\t\trunning model to scale " << scale);
         model->model->run_to(scale, get_precision());
         VERBOSE_MSG("> \t\t\tapplying constraint");
         constraint->apply();
      }
      // apply matching condition if this is not the last model
      if (m != models.size() - 1) {
         VERBOSE_MSG("> \tmatching to model " << models[m + 1]->model->name());
         Matching<Two_scale>* mc = model->matching_condition;
         mc->match_low_to_high_scale_model();
      }
   }
   VERBOSE_MSG("> running up finished");
}

void RGFlow<Two_scale>::run_down()
{
   VERBOSE_MSG("< running tower down ...");
   for (long m = models.size() - 1; m >= 0; --m) {
      TModel* model = models[m];
      VERBOSE_MSG("< \tselecting model " << model->model->name());
      // apply all constraints:
      // If m is the last model, do not apply the highest constraint,
      // because it was already appied when we ran up.
      const long c_begin = (m == static_cast<long>(models.size() - 1) ? 1 : 0);
      for (long c = c_begin; c < model->downwards_constraints.size(); ++c) {
         Constraint<Two_scale>* constraint = model->downwards_constraints[c];
         const double scale = constraint->get_scale();
         VERBOSE_MSG("< \t\tselecting constraint " << c << " at scale " << scale);
         VERBOSE_MSG("< \t\t\trunning model to scale " << scale);
         model->model->run_to(scale, get_precision());
         // If m is the first model, do not apply the lowest
         // constraint, because it will be appied when we run up next
         // time.
         if (m != 0 || c != model->downwards_constraints.size() - 1) {
            VERBOSE_MSG("< \t\t\tapplying constraint");
            constraint->apply();
         }
      }
      // apply matching condition if this is not the first model
      if (m > 0) {
         Matching<Two_scale>* mc = models[m - 1]->matching_condition;
         VERBOSE_MSG("< \tmatching to model " << models[m - 1]->model->name());
         mc->match_high_to_low_scale_model();
      }
   }
   VERBOSE_MSG("< running down finished");
}

void RGFlow<Two_scale>::apply_lowest_constaint()
{
   if (models.empty())
      return;

   TModel* model = models[0];

   if (model->downwards_constraints.empty())
      return;

   Constraint<Two_scale>* constraint = model->downwards_constraints.back();
   const double scale = constraint->get_scale();
   VERBOSE_MSG("| selecting constraint 0 at scale " << scale);
   VERBOSE_MSG("| \trunning model " << model->model->name() << " to scale " << scale);
   model->model->run_to(scale, get_precision());
   VERBOSE_MSG("| \tapplying constraint");
   constraint->apply();
}

/**
 * Returns the precision of the RG running.  If
 * use_increasing_precision == false the function returns the default
 * precision \f$1.0e-3\f$.  Otherwise the returned precision is
 * \f$1/10^n\f$ where \f$n\f$ is the iteration number (starting at 1).
 *
 * @return RG running precision
 */
double RGFlow<Two_scale>::get_precision()
{
   static const double default_precision = 1.0e-3;
   static const double minimum_precision = default_precision * 0.01;

   if (use_increasing_precision) {
      return std::max(exp(- static_cast<double>(iteration + 1) * log(10.0)),
                      minimum_precision);
   }

   return default_precision;
}

/**
 * Add a model and the corresponding model constraints.  Note that the
 * order of the model registration is important: Models that are added
 * later are assumed to be valid at a higher scale.
 *
 * @param model model
 * @param constraints vector of model constraints
 */
void RGFlow<Two_scale>::add_model(Two_scale_model* model,
                                  const std::vector<Constraint<Two_scale>*>& constraints)
{
   add_model(model, NULL, constraints);
}

void RGFlow<Two_scale>::add_model(Two_scale_model* model,
                                  const std::vector<Constraint<Two_scale>*>& upwards_constraints,
                                  const std::vector<Constraint<Two_scale>*>& downwards_constraints)
{
   add_model(model, NULL, upwards_constraints, downwards_constraints);
}

/**
 * Add a model, the corresponding model constraints and the matching
 * condition to the next model.  Note that the order of the model
 * registration is important: Models that are added later are assumed
 * to be valid at a higher scale.
 *
 * @param model model
 * @param mc matching condition to the next higher model
 * @param constraints vector of model constraints
 */
void RGFlow<Two_scale>::add_model(Two_scale_model* model,
                                  Matching<Two_scale>* mc,
                                  const std::vector<Constraint<Two_scale>*>& constraints)
{
   // create vector of downward constraints
   std::vector<Constraint<Two_scale>*> downward_constraints;
   std::reverse_copy(constraints.begin(), constraints.end(),
                     std::back_inserter(downward_constraints));

   add_model(model, mc, constraints, downward_constraints);
}

void RGFlow<Two_scale>::add_model(Two_scale_model* model,
                                  Matching<Two_scale>* mc,
                                  const std::vector<Constraint<Two_scale>*>& upwards_constraints,
                                  const std::vector<Constraint<Two_scale>*>& downwards_constraints)
{
   models.push_back(new TModel(model, upwards_constraints, downwards_constraints, mc));
}

bool RGFlow<Two_scale>::accuracy_goal_reached() const
{
   if (convergence_tester)
      return convergence_tester->accuracy_goal_reached();
   return false;
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

/**
 * Enable/ Disable the increasing RG running precision in each
 * iteration.
 *
 * @param use_increasing_precision_ enable (true), disable (false)
 * increasing precision
 */
void RGFlow<Two_scale>::set_increasing_running_precision(bool use_increasing_precision_)
{
   use_increasing_precision = use_increasing_precision_;
}

void RGFlow<Two_scale>::set_initial_guesser(Initial_guesser<Two_scale>* ig)
{
   initial_guesser = ig;
}

/**
 * Set the maximum number of iterations.  If 0 is given, there will be
 * no iteration at all.
 *
 * @param max_it maximum number of iterations
 */
void RGFlow<Two_scale>::set_max_iterations(unsigned int max_it)
{
   max_iterations = max_it;
}

unsigned int RGFlow<Two_scale>::number_of_iterations_done() const
{
   return iteration;
}
