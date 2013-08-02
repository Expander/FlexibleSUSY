#include "mssm_parameter_point.hpp"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale_sugra_constraint.hpp"
#include "mssm_two_scale_susy_scale_constraint.hpp"
#include "mssm_two_scale_low_scale_constraint.hpp"
#include "mssm_two_scale_convergence_tester.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"
#include "logger.hpp"
#include "coupling_monitor.hpp"

using namespace flexiblesusy;

int main()
{
   Mssm_parameter_point pp;

   Mssm<Two_scale> mssm;
   Mssm_sugra_constraint mssm_sugra_constraint(pp);
   Mssm_low_scale_constraint mssm_mz_constraint(pp);
   Mssm_susy_scale_constraint mssm_msusy_constraint(pp);
   Mssm_convergence_tester mssm_convergence_tester(&mssm, 1.0e-4);
   Mssm_initial_guesser initial_guesser(&mssm, pp, mssm_mz_constraint,
                                        mssm_msusy_constraint,
                                        mssm_sugra_constraint);
   Two_scale_increasing_precision two_scale_increasing_precision(10.0, 1.0e-5);

   std::vector<Constraint<Two_scale>*> mssm_upward_constraints;
   mssm_upward_constraints.push_back(&mssm_mz_constraint);
   mssm_upward_constraints.push_back(&mssm_sugra_constraint);

   std::vector<Constraint<Two_scale>*> mssm_downward_constraints;
   mssm_downward_constraints.push_back(&mssm_sugra_constraint);
   mssm_downward_constraints.push_back(&mssm_msusy_constraint);
   mssm_downward_constraints.push_back(&mssm_mz_constraint);

   RGFlow<Two_scale> solver;
   solver.set_convergence_tester(&mssm_convergence_tester);
   solver.set_running_precision(&two_scale_increasing_precision);
   solver.set_initial_guesser(&initial_guesser);
   solver.add_model(&mssm, mssm_upward_constraints, mssm_downward_constraints);

   INFO("Running: " << pp);
   try {
      solver.solve();
   } catch (Error& e) {
      ERROR("no solution found: " << e.what());
      exit(1);
   }

   mssm.calculate_spectrum();

   INFO("Solution found: ");
   mssm.print(std::cout);
}
