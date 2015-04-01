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

#ifndef MSSMD5O_MSSMRHN_SPECTRUM_GENERATOR_H
#define MSSMD5O_MSSMRHN_SPECTRUM_GENERATOR_H

#include "MSSMD5O_two_scale_model.hpp"
#include "MSSMD5O_two_scale_susy_scale_constraint.hpp"
#include "MSSMD5O_two_scale_low_scale_constraint.hpp"
#include "MSSMD5O_two_scale_convergence_tester.hpp"
#include "MSSMD5O_utilities.hpp"

#include "MSSMD5O_MSSMRHN_two_scale_matching.hpp"

#include "MSSMRHN_two_scale_model.hpp"
#include "MSSMRHN_two_scale_high_scale_constraint.hpp"
#include "MSSMRHN_two_scale_convergence_tester.hpp"
#include "MSSMRHN_utilities.hpp"

#include "MSSMD5O_MSSMRHN_two_scale_initial_guesser.hpp"

#include "coupling_monitor.hpp"
#include "error.hpp"
#include "numerics.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_composite_convergence_tester.hpp"

namespace flexiblesusy {

template <class T>
class MSSMD5O_MSSMRHN_spectrum_generator {
public:
   MSSMD5O_MSSMRHN_spectrum_generator()
      : solver()
      , model_1(), model_2()
      , high_scale_constraint_2()
      , susy_scale_constraint_1()
      , low_scale_constraint_1()
      , matching()
      , high_scale_2(0.), susy_scale_1(0.), low_scale_1(0.)
      , matching_scale(0)
      , parameter_output_scale_1(0.)
      , precision_goal(1.0e-4)
      , max_iterations(0)
      , calculate_sm_masses(false) {}
   ~MSSMD5O_MSSMRHN_spectrum_generator() {}

   double get_high_scale_2() const { return high_scale_2; }
   double get_susy_scale_1() const { return susy_scale_1; }
   double get_low_scale_1()  const { return low_scale_1;  }
   double get_matching_scale() const { return matching_scale; }
   const MSSMD5O<T>& get_model_1() const { return model_1; }
   const MSSMRHN<T>& get_model_2() const { return model_2; }
   const Problems<MSSMD5O_info::NUMBER_OF_PARTICLES>& get_problems() const {
      return model_1.get_problems();
   }
   int get_exit_code() const { return get_problems().have_problem(); };
   void set_parameter_output_scale_1(double s) { parameter_output_scale_1 = s; }
   void set_precision_goal(double precision_goal_) { precision_goal = precision_goal_; }
   void set_pole_mass_loop_order(unsigned l) { model_1.set_pole_mass_loop_order(l); }
   void set_ewsb_loop_order(unsigned l) { model_1.set_ewsb_loop_order(l); }
   void set_max_iterations(unsigned n) { max_iterations = n; }
   void set_calculate_sm_masses(bool flag) { calculate_sm_masses = flag; }

   void run(const QedQcd& oneset, const MSSMD5O_input_parameters& input_1, const MSSMRHN_input_parameters& input_2);
   void write_running_couplings_1(const std::string& filename = "MSSMD5O_rge_running.dat") const;
   void write_spectrum(const std::string& filename = "MSSMD5O_spectrum.dat") const;

private:
   RGFlow<T> solver;
   MSSMD5O<T> model_1;
   MSSMRHN<T> model_2;
   MSSMRHN_high_scale_constraint<T> high_scale_constraint_2;
   MSSMD5O_susy_scale_constraint<T> susy_scale_constraint_1;
   MSSMD5O_low_scale_constraint<T> low_scale_constraint_1;
   MSSMD5O_MSSMRHN_matching<T> matching;
   double high_scale_2, susy_scale_1, low_scale_1;
   double matching_scale;
   double parameter_output_scale_1; ///< output scale for running parameters
   double precision_goal; ///< precision goal
   unsigned max_iterations; ///< maximum number of iterations
   bool calculate_sm_masses; ///< calculate SM pole masses
};

/**
 * @brief Run's the RG solver with the given input parameters
 *
 * This function sets up the RG solver using a high-scale, susy-scale
 * and low-scale constraint.  Afterwards the solver is run until
 * convergence is reached or an error occours.  Finally the particle
 * spectrum (pole masses) is calculated.
 *
 * @param oneset Standard Model input parameters
 * @param input model input parameters
 */
template <class T>
void MSSMD5O_MSSMRHN_spectrum_generator<T>::run
(const QedQcd& oneset,
 const MSSMD5O_input_parameters& input_1, const MSSMRHN_input_parameters& input_2)
{
   high_scale_constraint_2.clear();
   susy_scale_constraint_1.clear();
   low_scale_constraint_1.clear();

   matching.reset();

   model_1.clear();
   model_1.set_input_parameters(input_1);
   model_1.do_calculate_sm_pole_masses(calculate_sm_masses);

   model_2.clear();
   model_2.set_input_parameters(input_2);

   // needed for constraint::initialize()
   high_scale_constraint_2.set_model(&model_2);
   susy_scale_constraint_1.set_model(&model_1);
   low_scale_constraint_1 .set_model(&model_1);

   low_scale_constraint_1.set_sm_parameters(oneset);
   matching.set_lower_input_parameters(input_1);
   high_scale_constraint_2.initialize();
   susy_scale_constraint_1.initialize();
   low_scale_constraint_1.initialize();

   std::vector<Constraint<T>*> upward_constraints_1;
   upward_constraints_1.push_back(&low_scale_constraint_1);

   std::vector<Constraint<T>*> downward_constraints_1;
   downward_constraints_1.push_back(&susy_scale_constraint_1);
   downward_constraints_1.push_back(&low_scale_constraint_1);

   std::vector<Constraint<T>*> upward_constraints_2;
   upward_constraints_2.push_back(&high_scale_constraint_2);

   std::vector<Constraint<T>*> downward_constraints_2;
   downward_constraints_2.push_back(&high_scale_constraint_2);

   MSSMD5O_convergence_tester<T> convergence_tester_1(&model_1, precision_goal);
   MSSMRHN_convergence_tester<T> convergence_tester_2(&model_2, precision_goal);
   if (max_iterations > 0) {
      convergence_tester_1.set_max_iterations(max_iterations);
      convergence_tester_2.set_max_iterations(max_iterations);
   }
   Composite_convergence_tester<T> convergence_tester;
   convergence_tester.add_convergence_tester(&convergence_tester_1);
   convergence_tester.add_convergence_tester(&convergence_tester_2);

   MSSMD5O_MSSMRHN_initial_guesser<T> initial_guesser
       (&model_1, &model_2, input_1, oneset,
	low_scale_constraint_1,
	susy_scale_constraint_1,
	high_scale_constraint_2,
	matching);
   Two_scale_increasing_precision precision(10.0, precision_goal);

   solver.reset();
   solver.set_convergence_tester(&convergence_tester);
   solver.set_running_precision(&precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_model(&model_1, &matching, upward_constraints_1, downward_constraints_1);
   solver.add_model(&model_2, upward_constraints_2, downward_constraints_2);

   high_scale_2 = susy_scale_1 = low_scale_1 = 0.;
   matching_scale = 0;

   try {
      solver.solve();
      high_scale_2 = high_scale_constraint_2.get_scale();
      susy_scale_1 = susy_scale_constraint_1.get_scale();
      low_scale_1 = low_scale_constraint_1.get_scale();
      matching_scale = matching.get_scale();

      model_1.run_to(susy_scale_1);
      model_1.calculate_spectrum();

      // run to output scale (if scale > 0)
      if (!is_zero(parameter_output_scale_1)) {
         model_1.run_to(parameter_output_scale_1);
      }
   } catch (const NoConvergenceError&) {
      model_1.get_problems().flag_no_convergence();
   } catch (const NonPerturbativeRunningError&) {
      model_1.get_problems().flag_no_perturbative();
   } catch (const Error& error) {
      model_1.get_problems().flag_thrown();
   } catch (const std::string& str) {
      model_1.get_problems().flag_thrown();
   } catch (const char* str) {
      model_1.get_problems().flag_thrown();
   } catch (const std::exception& error) {
      model_1.get_problems().flag_thrown();
   }
}

/**
 * Create a text file which contains the values of all model
 * parameters at all scales between the low-scale and the high-scale.
 *
 * @param filename name of output file
 */
template <class T>
void MSSMD5O_MSSMRHN_spectrum_generator<T>::write_running_couplings_1(const std::string& filename) const
{
   MSSMD5O<T> tmp_model(model_1);
   tmp_model.run_to(low_scale_1);

   MSSMD5O_parameter_getter parameter_getter;
   Coupling_monitor<MSSMD5O<T>, MSSMD5O_parameter_getter>
      coupling_monitor(tmp_model, parameter_getter);

   coupling_monitor.run(low_scale_1, matching_scale, 100, true);
   coupling_monitor.write_to_file(filename);
}

/**
 * Write spectrum (pole masses) to a text file
 *
 * @param filename output file name
 */
template <class T>
void MSSMD5O_MSSMRHN_spectrum_generator<T>::write_spectrum(const std::string& filename) const
{
   MSSMD5O_spectrum_plotter plotter;
   plotter.extract_spectrum<T>(model_1);
   plotter.write_to_file(filename);
}

} // namespace flexiblesusy

#endif
