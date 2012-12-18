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

#include <vector>
#include <string>
#include <sstream>

template <class T> class Constraint;
template <class T> class Matching;
template <class T> class Convergence_tester;
template <class T> class Initial_guesser;
class Two_scale;
class Two_scale_model;

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
   /// add models and up- and downwards constraints
   void add_model(Two_scale_model*,
                  const std::vector<Constraint<Two_scale>*>&,
                  const std::vector<Constraint<Two_scale>*>&);
   /// add models, constraints and matching condition
   void add_model(Two_scale_model*,
                  Matching<Two_scale>* m = NULL,
                  const std::vector<Constraint<Two_scale>*>& constraints = std::vector<Constraint<Two_scale>*>());
   /// add models, up- and downward constraints and matching condition
   void add_model(Two_scale_model*,
                  Matching<Two_scale>* m,
                  const std::vector<Constraint<Two_scale>*>& upwards_constraints,
                  const std::vector<Constraint<Two_scale>*>& downwards_constraints);
   /// get number of used iterations
   unsigned int number_of_iterations_done() const;
   /// set convergence tester
   void set_convergence_tester(Convergence_tester<Two_scale>*);
   /// set increasing running precision
   void set_increasing_running_precision(bool);
   /// set initial guesser
   void set_initial_guesser(Initial_guesser<Two_scale>*);
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
      std::vector<Constraint<Two_scale>*> upwards_constraints; ///< upwards constraints
      std::vector<Constraint<Two_scale>*> downwards_constraints; ///< downwards constraints
      Matching<Two_scale>* matching_condition;         ///< matching condition

      TModel(Two_scale_model* m,
             const std::vector<Constraint<Two_scale>*>& uc,
             const std::vector<Constraint<Two_scale>*>& dc,
             Matching<Two_scale>* mc)
         : model(m)
         , upwards_constraints(uc)
         , downwards_constraints(dc)
         , matching_condition(mc)
         {}
   };
   std::vector<TModel*> models;        ///< tower of models (from low to high scale)
   unsigned int max_iterations;        ///< maximum number of iterations
   unsigned int iteration;             ///< iteration number (starting at 0)
   Convergence_tester<Two_scale>* convergence_tester; ///< the convergence tester
   Initial_guesser<Two_scale>* initial_guesser;       ///< does initial guess
   bool use_increasing_precision;      ///< use increasing RG running precision

   bool accuracy_goal_reached() const; ///< check if accuracy goal is reached
   void check_setup() const;           ///< check the setup
   void initial_guess();               ///< initial guess
   void run_up();                      ///< run all models up
   void run_down();                    ///< run all models down
   void apply_lowest_constaint();      ///< apply lowest constraint
   double get_precision();             ///< returns running precision
};

#endif
