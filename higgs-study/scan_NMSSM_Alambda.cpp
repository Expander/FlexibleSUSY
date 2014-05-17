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

// File generated at Thu 15 May 2014 20:35:09

#include "NMSSM_input_parameters.hpp"
#include "NMSSM_spectrum_generator.hpp"

#include "NMSSMNoUni_input_parameters.hpp"
#include "NMSSMNoUni_spectrum_generator.hpp"

#include "error.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "options.hpp"

#include <iostream>

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   Options options(argc, argv);

   NMSSM_input_parameters input;
   NMSSMNoUni_input_parameters input_nouni;
   QedQcd oneset;
   oneset.toMz();

   NMSSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-3);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   NMSSMNoUni_spectrum_generator<algorithm_type> spectrum_generator_nouni;
   spectrum_generator_nouni.set_precision_goal(1.0e-3);
   spectrum_generator_nouni.set_max_iterations(0);         // 0 == automatic
   spectrum_generator_nouni.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator_nouni.set_parameter_output_scale(0); // 0 == susy scale

   const std::vector<double> range_TanBeta(
      float_range(options.tanb_start, options.tanb_stop, options.tanb_npoints));

   input.m0 = 5000.;
   input.m12 = 5000.;
   input.TanBeta = 10.;
   input.Azero = 5000.;
   input.LambdaInput = 0.2;

   input_nouni.m0 = input.m0;
   input_nouni.m12 = input.m12;
   input_nouni.TanBeta = input.TanBeta;
   input_nouni.Azero = input.Azero;
   input_nouni.LambdaInput = input.LambdaInput;
   input_nouni.ALambdaInput = 10000.;

   cout << "# NMSSM with "
        << "m0 = " << input.m0
        << ", m12 = " << input.m12
        << ", Azero = " << input.Azero
        << ", LambdaInput = " << input.LambdaInput
        << ", ALambdaInput = " << input_nouni.ALambdaInput
        << '\n';

   cout << "# "
        << std::setw(12) << std::left << "TanBeta" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "error" << ' '
        << std::setw(12) << std::left << "Mhh_nouni(1)/GeV" << ' '
        << std::setw(12) << std::left << "error_nouni"
        << '\n';

   for (auto tanBeta : range_TanBeta) {
      input.TanBeta       = tanBeta;
      input_nouni.TanBeta = tanBeta;

      spectrum_generator.run(oneset, input);
      spectrum_generator_nouni.run(oneset, input_nouni);

      const NMSSM<algorithm_type>& model = spectrum_generator.get_model();
      const NMSSM_physical& pole_masses = model.get_physical();
      const Problems<NMSSM_info::NUMBER_OF_PARTICLES>& problems
         = spectrum_generator.get_problems();
      const double higgs = pole_masses.Mhh(0);
      const bool error = problems.have_serious_problem();

      const NMSSMNoUni<algorithm_type>& model_nouni = spectrum_generator_nouni.get_model();
      const NMSSMNoUni_physical& pole_masses_nouni = model_nouni.get_physical();
      const Problems<NMSSMNoUni_info::NUMBER_OF_PARTICLES>& problems_nouni
         = spectrum_generator_nouni.get_problems();
      const double higgs_nouni = pole_masses_nouni.Mhh(0);
      const bool error_nouni = problems_nouni.have_serious_problem();

      cout << "  "
           << std::setw(12) << std::left << input.TanBeta << ' '
           << std::setw(12) << std::left << higgs << ' '
           << std::setw(12) << std::left << error << ' '
           << std::setw(12) << std::left << higgs_nouni << ' '
           << std::setw(12) << std::left << error_nouni;
      if (error) {
         cout << "\t# " << problems;
      }
      if (error_nouni) {
         cout << "\t#NoUni: " << problems_nouni;
      }
      cout << '\n';
   }

   return 0;
}
