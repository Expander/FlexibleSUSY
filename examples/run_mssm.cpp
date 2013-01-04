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

class Sfermion_masses_getter {
public:
   template <class Rge>
   DoubleVector operator()(const Rge& rge) {
      SoftParsMssm softPars(static_cast<MssmSoftsusy>(rge).displaySoftPars());

      DoubleVector mQ      (softPars.displaySoftMassSquared(mQl).flatten());
      DoubleVector mU      (softPars.displaySoftMassSquared(mUr).flatten());
      DoubleVector mD      (softPars.displaySoftMassSquared(mUr).flatten());
      DoubleVector mL      (softPars.displaySoftMassSquared(mLl).flatten());
      DoubleVector mE      (softPars.displaySoftMassSquared(mEr).flatten());

      DoubleVector data(mQ.apply(sqrt));
      data.append(mU.apply(sqrt));
      data.append(mD.apply(sqrt));
      data.append(mL.apply(sqrt));
      data.append(mE.apply(sqrt));

      return data;
   }
};

class Gaugino_masses_getter {
public:
   template <class Rge>
   DoubleVector operator()(const Rge& rge) {
      SoftParsMssm softPars(static_cast<MssmSoftsusy>(rge).displaySoftPars());
      return softPars.displayGaugino();
   }
};

class Triliear_masses_getter {
public:
   template <class Rge>
   DoubleVector operator()(const Rge& rge) {
      SoftParsMssm softPars(static_cast<MssmSoftsusy>(rge).displaySoftPars());

      DoubleVector Au      (softPars.displayTrilinear(UA).flatten());
      DoubleVector Ad      (softPars.displayTrilinear(DA).flatten());
      DoubleVector Ae      (softPars.displayTrilinear(EA).flatten());

      DoubleVector data(Au);
      data.append(Ad);
      data.append(Ae);

      return data;
   }
};

class Soft_masses_getter {
public:
   template <class Rge>
   DoubleVector operator()(const Rge& rge) {
      Sfermion_masses_getter smg;
      Gaugino_masses_getter  gmg;
      Triliear_masses_getter tmg;

      SoftParsMssm softPars(static_cast<MssmSoftsusy>(rge).displaySoftPars());
      DoubleVector mH(2);
      mH(1) = sqrt(softPars.displayMh1Squared());
      mH(2) = sqrt(softPars.displayMh2Squared());

      DoubleVector data(smg(rge));
      data.append(gmg(rge));
      data.append(tmg(rge));
      data.append(mH);

      return data;
   }
};

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

   Coupling_monitor cm;
   Soft_masses_getter smg;

   mssm.run_to(MZ);
   cm.run(mssm, smg, MZ, gut_scale, 100, true);
   cm.write_to_file("mssm_soft_masses.dat");
#endif
}
