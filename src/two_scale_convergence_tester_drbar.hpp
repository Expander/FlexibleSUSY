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

#ifndef TWO_SCALE_CONVERGENCE_TESTER_DRBAR_H
#define TWO_SCALE_CONVERGENCE_TESTER_DRBAR_H

#include "two_scale_convergence_tester.hpp"
#include "convergence_tester_drbar.hpp"
#include "error.hpp"
#include "logger.hpp"
#include "numerics2.hpp"

#include <cmath>
#include <limits>

namespace flexiblesusy {

template <template<class Method> class Model>
class Convergence_tester_DRbar<Model<Two_scale> > :
	public Convergence_tester<Two_scale> {
public:
   using Scale_getter = std::function<double()>;

   Convergence_tester_DRbar(const Model<Two_scale>*, double, const Scale_getter& sg = Scale_getter());
   virtual ~Convergence_tester_DRbar() {}

   virtual bool accuracy_goal_reached() override;
   virtual double get_accuracy_goal() const { return accuracy_goal; }
   virtual unsigned int max_iterations() const override { return max_it; }
   double get_current_accuracy() const { return current_accuracy; }
   /// set maximum number of iterations
   void set_max_iterations(unsigned it) { max_it = it; }

protected:
   /// get current iteration number
   unsigned get_iteration() const { return it_count; }
   /// get model at current iteration
   const Model<Two_scale>& get_current_iteration_model() const { return current_model; }
   /// get model state during last iteration
   const Model<Two_scale>& get_last_iteration_model() const { return last_iteration_model; }
   /// maximum relative difference to last iteration
   virtual double max_rel_diff() const = 0;

private:
   const Model<Two_scale>* model{nullptr};  ///< pointer to model
   Model<Two_scale> current_model{};        ///< model state at current iteration
   Model<Two_scale> last_iteration_model{}; ///< model state at last iteration
   Scale_getter scale_getter{};             ///< function to retrieve scale
   unsigned int it_count{0};                ///< iteration
   unsigned int max_it{40};                 ///< maximum number of iterations
   double accuracy_goal{1e-4};              ///< accuracy goal
   double current_accuracy{std::numeric_limits<double>::infinity()}; ///< current accuracy

   double scale_difference() const;         ///< absolute scale difference
   double rel_scale_difference() const;     ///< relative scale difference
   double run_to_scale();                   ///< runs models to comparison scale
};

template <template<class Method> class Model>
Convergence_tester_DRbar<Model<Two_scale> >::Convergence_tester_DRbar
(const Model<Two_scale>* model_, double accuracy_goal_, const Scale_getter& sg)
   : Convergence_tester<Two_scale>()
   , model(model_)
   , max_it(static_cast<int>(-log10(accuracy_goal_) * 10))
   , accuracy_goal(accuracy_goal_)
   , scale_getter(sg)
{
   if (!model)
      throw SetupError("Convergence_tester_DRbar<Model<Two_scale>>: "
                       "model pointer must not be zero!");
}

template <template<class Method> class Model>
bool Convergence_tester_DRbar<Model<Two_scale> >::accuracy_goal_reached()
{
   current_model = *model;
   bool precision_reached = false;

   if (it_count > 0) {
      const double scale_accuracy_goal = accuracy_goal * 16*M_PI*M_PI;
      if (rel_scale_difference() < scale_accuracy_goal) {
         run_to_scale();
	 current_accuracy = max_rel_diff();
	 precision_reached = current_accuracy < accuracy_goal;
	 VERBOSE_MSG("Convergence_tester_DRbar: current accuracy = "
		     << current_accuracy
		     << ", accuracy goal = " << accuracy_goal);
      } else {
         VERBOSE_MSG("scale has changed by " << scale_difference()
		     << " GeV (" << rel_scale_difference()
		     << "), skipping parameter comparison");
      }
   }

   // save old model parameters
   last_iteration_model = current_model;
   ++it_count;

   return precision_reached;
}

template <template<class Method> class Model>
double Convergence_tester_DRbar<Model<Two_scale> >::scale_difference() const
{
   return current_model.get_scale() - last_iteration_model.get_scale();
}

template <template<class Method> class Model>
double Convergence_tester_DRbar<Model<Two_scale> >::rel_scale_difference()
   const
{
   const double diff = scale_difference();
   const double last_scale = last_iteration_model.get_scale();
   if (!is_zero(last_scale))
      return diff / last_scale;
   return std::numeric_limits<double>::infinity();
}

template <template<class Method> class Model>
double Convergence_tester_DRbar<Model<Two_scale> >::run_to_scale()
{
   if (scale_getter) {
      const double scale = scale_getter();
      current_model.run_to(scale);
      current_model.calculate_DRbar_masses();
      last_iteration_model.run_to(scale);
      last_iteration_model.calculate_DRbar_masses();
   }
}

} // namespace flexiblesusy

#endif
