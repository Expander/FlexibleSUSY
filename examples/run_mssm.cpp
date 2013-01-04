#include "mssm_parameter_point.hpp"
#include "mssm_two_scale.hpp"
#include "mssm_two_scale_initial_guesser.hpp"
#include "mssm_two_scale_sugra_constraint.hpp"
#include "mssm_two_scale_msusy_constraint.hpp"
#include "mssm_two_scale_mz_constraint.hpp"
#include "mssm_two_scale_convergence_tester.hpp"
#include "two_scale_running_precision.hpp"
#include "two_scale_solver.hpp"
#include "logger.hpp"
#include "coupling_monitor.hpp"

int main()
{
   Mssm_parameter_point pp;

   Mssm<Two_scale> mssm;
   Mssm_sugra_constraint mssm_sugra_constraint(&mssm, pp.mxGuess, pp.m0, pp.m12, pp.a0, pp.signMu);
   Mssm_mz_constraint mssm_mz_constraint(&mssm, pp.tanBeta);
   Mssm_msusy_constraint mssm_msusy_constraint(&mssm, pp.get_soft_pars(), 1000.0, pp.signMu);
   Mssm_convergence_tester mssm_convergence_tester(&mssm, 1.0e-4);
   Mssm_initial_guesser initial_guesser(&mssm, pp.oneset, pp.mxGuess, pp.tanBeta, pp.signMu, pp.get_soft_pars(), false);
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

#if 0
   // create data: all gauge couplings at different scales
   const double gut_scale = mssm_sugra_constraint.get_scale();
   const double MZ = 91.1876;

   mssm.run_to(MZ);

   Coupling_monitor cm;
   Gauge_coupling_getter gcg;
   cm.run(mssm, gcg, MZ, gut_scale, 100, true);
   cm.write_to_file("mssm_running_coupling.dat");
#endif
}
