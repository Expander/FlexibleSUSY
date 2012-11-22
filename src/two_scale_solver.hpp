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

#ifndef TWO_SCALE_SOLVER_H
#define TWO_SCALE_SOLVER_H

#include "rg_flow.hpp"
#include "two_scale_model.hpp"
#include "two_scale_constraint.hpp"
#include "two_scale_matching.hpp"
#include "two_scale_convergence_tester.hpp"
#include "logger.hpp"

#include <vector>
#include <string>
#include <sstream>

class Two_scale;

template<>
class RGFlow<Two_scale> {
public:
   class Error {
   public:
      virtual ~Error() {}
      virtual std::string what() const = 0;
   };

   class SetupError : public Error {
   public:
      SetupError(const std::string& message_) : message(message_) {}
      virtual ~SetupError() {}
      virtual std::string what() const { return message; }
   private:
      std::string message;
   };

   class NoConvergenceError : public Error {
   public:
      NoConvergenceError(unsigned number_of_iterations_)
         : number_of_iterations(number_of_iterations_) {}
      virtual ~NoConvergenceError() {}
      virtual std::string what() const {
         std::stringstream message;
         message << "RGFlow<Two_scale>::NoConvergenceError: no convergence"
                 << " after " << number_of_iterations << " iterations";
         return message.str();
      }
   private:
      unsigned number_of_iterations;
   };


   RGFlow();
   ~RGFlow();

   /// add models and constraints
   void add_model(Two_scale_model*,
                  const std::vector<Constraint<Two_scale>*>&);
   /// add models, constraints and matching condition
   void add_model(Two_scale_model*,
                  Matching<Two_scale>* m = NULL,
                  const std::vector<Constraint<Two_scale>*>& constraints = std::vector<Constraint<Two_scale>*>());
   /// set convergence tester
   void set_convergence_tester(Convergence_tester<Two_scale>*);
   /// set maximum number of iterations
   void set_max_iterations(unsigned int);
   /// solve all models
   void solve();

private:
   /**
    * @class TModel
    * @brief contains model, constraints and matching condition
    *
    * This class lumps together the model, its constraints and the
    * matching condition to the next higher model.
    */
   struct TModel {
      Two_scale_model* model;                          ///< the model
      std::vector<Constraint<Two_scale>*> constraints; ///< model constraints
      const Matching<Two_scale>* matching_condition;   ///< matching condition

      TModel(Two_scale_model* m,
             const std::vector<Constraint<Two_scale>*>& c,
             Matching<Two_scale>* mc)
         : model(m)
         , constraints(c)
         , matching_condition(mc)
         {}
   };
   std::vector<TModel*> models;        ///< tower of models (from low to high scale)
   unsigned int max_iterations;        ///< maximum number of iterations
   Convergence_tester<Two_scale>* convergence_tester; ///< the convergence tester

   bool accuracy_goal_reached() const; ///< check if accuracy goal is reached
   void check_setup() const;           ///< check the setup
   void run_up();                      ///< run all models up
   void run_down();                    ///< run all models down
};

/**
 * Create empty two scale solver.  Sets maximum number of iterations
 * to 10 (default).
 */
inline RGFlow<Two_scale>::RGFlow()
   : models()
   , max_iterations(10)
   , convergence_tester(NULL)
{
}

inline RGFlow<Two_scale>::~RGFlow()
{
   for (size_t m = 0; m < models.size(); ++m)
      delete models[m];
}

inline void RGFlow<Two_scale>::solve()
{
   check_setup();

   // initial run
   run_up();
   run_down();

   unsigned int iter = max_iterations;
   while (iter && !accuracy_goal_reached()) {
      run_up();
      run_down();
      --iter;
   }

   if (iter == 0 && convergence_tester)
      throw NoConvergenceError(max_iterations);
}

inline void RGFlow<Two_scale>::check_setup() const
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
            WARNING("the matching condition of the last model [" << m
                    << "] is non-zero but will not be used");
      } else {
         if (model->matching_condition == NULL) {
            std::stringstream message;
            message << "RGFlow<Two_scale>::Error: matching condition "
                    << "of model " << m << " pointer is NULL";
            throw SetupError(message.str());
         }
      }
   }
}

inline void RGFlow<Two_scale>::run_up()
{
   for (size_t m = 0; m < models.size(); ++m) {
      TModel* model = models[m];
      // apply all constraints
      for (size_t c = 0; c < model->constraints.size(); ++c) {
         Constraint<Two_scale>* constraint = model->constraints[c];
         const double scale = constraint->get_scale();
         model->model->run_to(scale);
         constraint->update_scale();
         constraint->apply();
      }
      // apply matching condition if this is not the last model
      if (m != models.size() - 1) {
         const Matching<Two_scale>* mc = model->matching_condition;
         model->model->run_to(mc->get_scale());
         mc->match_low_to_high_scale_model();
      }
   }
}

inline void RGFlow<Two_scale>::run_down()
{
   for (long m = models.size() - 1; m >= 0; --m) {
      TModel* model = models[m];
      // apply all constraints:
      // If m is the last model, do not apply the highest mc, because
      // it was already appied when we ran up.
      // If m is the first model, do not apply the lowest mc, because
      // it will be appied when we run up next time.
      const long c_begin = (m == static_cast<long>(models.size() - 1) ?
                            model->constraints.size() - 2 :
                            model->constraints.size() - 1);
      const long c_end = (m == 0 ? 1 : 0);
      for (long c = c_begin; c >= c_end; --c) {
         Constraint<Two_scale>* constraint = model->constraints[c];
         const double scale = constraint->get_scale();
         model->model->run_to(scale);
         constraint->update_scale();
         constraint->apply();
      }
      // apply matching condition if this is not the first model
      if (m > 0) {
         const Matching<Two_scale>* mc = models[m - 1]->matching_condition;
         model->model->run_to(mc->get_scale());
         mc->match_high_to_low_scale_model();
      }
   }
}

/**
 * Add a model and the corresponding model constraints.  Note that the
 * order of the model registration is important: Models that are added
 * later are assumed to be valid at a higher scale.
 *
 * @param model model
 * @param constraints vector of model constraints
 */
inline void RGFlow<Two_scale>::add_model(Two_scale_model* model,
                                         const std::vector<Constraint<Two_scale>*>& constraints)
{
   models.push_back(new TModel(model, constraints, NULL));
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
inline void RGFlow<Two_scale>::add_model(Two_scale_model* model,
                                         Matching<Two_scale>* mc,
                                         const std::vector<Constraint<Two_scale>*>& constraints)
{
   models.push_back(new TModel(model, constraints, mc));
}

inline bool RGFlow<Two_scale>::accuracy_goal_reached() const
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
inline void RGFlow<Two_scale>::set_convergence_tester(Convergence_tester<Two_scale>* convergence_tester_)
{
   convergence_tester = convergence_tester_;
}

/**
 * Set the maximum number of iterations (excluding the initial run).
 * If 0 is given, there will be only an initial run.
 *
 * @param max_it maximum number of iterations
 */
inline void RGFlow<Two_scale>::set_max_iterations(unsigned int max_it)
{
   max_iterations = max_it;
}

#endif
