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

#include <memory>
#include <vector>
#include <string>

/**
 * @file two_scale_solver.hpp
 * @brief contains the definition of the RGFlow<Two_scale> class
 */

namespace flexiblesusy {

template <class T> class Constraint;
template <class T> class Matching;
template <class T> class Convergence_tester;
template <class T> class Initial_guesser;
class Two_scale;
class Two_scale_model;
class Two_scale_running_precision;

/**
 * @class RGFlow<Two_scale>
 * @brief Boundary condition solver (two-scale algorithm)
 *
 * This boundary condition solver uses the two-scale algorithm to
 * solve the boundary value problem: It uses RG running to iteratively
 * run the models to the boundary condtion (constraint) scales and
 * imposes the constraints.
 *
 * To add constraints use the add() function.  Matching conditions are
 * added using the add_upwards() or add_downwards() functions,
 * depending on whether the low-scale model should be matched to the
 * high-scale one (add_upwards()) or vice versa (add_downwards()).
 * The added constraints and matching conditions are applied in their
 * given order.
 */

template<>
class RGFlow<Two_scale> {
public:
   RGFlow();
   ~RGFlow();

   /// add constraint
   void add(Constraint<Two_scale>*, Two_scale_model*);
   /// add upwards matching condition
   void add_upwards(Matching<Two_scale>*, Two_scale_model*, Two_scale_model*);
   /// add downwards matching condition
   void add_downwards(Matching<Two_scale>*, Two_scale_model*, Two_scale_model*);

   /// get number of used iterations
   unsigned int number_of_iterations_done() const;
   /// clear all internal data
   void reset();
   /// set convergence tester
   void set_convergence_tester(Convergence_tester<Two_scale>*);
   /// set running precision calculator
   void set_running_precision(Two_scale_running_precision*);
   /// set initial guesser
   void set_initial_guesser(Initial_guesser<Two_scale>*);
   /// solves the boundary value problem
   void solve();

private:
   struct Slider {
   public:
      virtual ~Slider() {}
      virtual void clear_problems() {}
      virtual void slide() {}
      virtual void set_precision(double) {}
   };

   struct Constraint_slider : public Slider {
   public:
      Constraint_slider(Two_scale_model* m, Constraint<Two_scale>* c)
         : model(m), constraint(c) {}
      virtual ~Constraint_slider() {}
      virtual void clear_problems();
      virtual void slide();
      virtual void set_precision(double);
   private:
      Two_scale_model* model;
      Constraint<Two_scale>* constraint;
   };

   struct Matching_up_slider : public Slider {
   public:
      Matching_up_slider(Two_scale_model* m1, Two_scale_model* m2, Matching<Two_scale>* mc)
         : model1(m1), model2(m2), matching(mc) {}
      virtual ~Matching_up_slider() {}
      virtual void clear_problems();
      virtual void slide();
      virtual void set_precision(double);
   private:
      Two_scale_model *model1, *model2;
      Matching<Two_scale>* matching;
   };

   struct Matching_down_slider : public Slider {
   public:
      Matching_down_slider(Two_scale_model* m1, Two_scale_model* m2, Matching<Two_scale>* mc)
         : model1(m1), model2(m2), matching(mc) {}
      virtual ~Matching_down_slider() {}
      virtual void clear_problems();
      virtual void slide();
      virtual void set_precision(double);
   private:
      Two_scale_model *model1, *model2;
      Matching<Two_scale>* matching;
   };

   std::vector<std::shared_ptr<Slider> > sliders; ///< sliders to be run up and down
   unsigned int iteration;             ///< iteration number (starting at 0)
   Convergence_tester<Two_scale>* convergence_tester; ///< the convergence tester
   Initial_guesser<Two_scale>* initial_guesser;       ///< does initial guess
   Two_scale_running_precision* running_precision_calculator; ///< RG running precision calculator
   double running_precision;           ///< RG running precision
   Two_scale_model* model_at_this_scale; ///< model at current scale

   RGFlow(const RGFlow&) {}
   bool accuracy_goal_reached() const; ///< check if accuracy goal is reached
   void check_setup() const;           ///< check the setup
   void clear_problems();              ///< clear model problems
   unsigned int get_max_iterations() const; ///< returns max. number of iterations
   void initial_guess();               ///< initial guess
   double get_precision();             ///< returns running precision
   void update_running_precision();    ///< update the RG running precision

   void run_sliders();                 ///< run all sliders
};

}

#endif
