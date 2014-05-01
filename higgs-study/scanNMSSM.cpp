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

// File generated at Thu 1 May 2014 12:16:53

#include "NMSSM_input_parameters.hpp"
#include "NMSSM_spectrum_generator.hpp"

#include "error.hpp"
#include "scan.hpp"
#include "lowe.h"
#include "options.hpp"

#include <iostream>

int main()
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   NMSSM_input_parameters input;
   QedQcd oneset;
   oneset.toMz();

   NMSSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-4);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   const std::vector<double> range_TanBeta(float_range(TANB_START, TANB_STOP, TANB_NPOINTS));
   const std::vector<double> range_m0(float_range(M0_START, M0_STOP, M0_NPOINTS));

   input.m0 = 5000.;
   input.m12 = 5000.;
   input.TanBeta = 10.;
   input.Azero = 5000.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;

   cout << "# @ModelName@ with "
        << "Azero = " << input.Azero
        << ", m12 = " << input.m12
        << ", LambdaInput = " << input.LambdaInput
        << ", SignvS = " << input.SignvS
        << '\n';

   cout << "# "
        << std::setw(12) << std::left << "TanBeta" << ' '
        << std::setw(12) << std::left << "m0" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "error"
        << '\n';

   for (auto tanBeta : range_TanBeta) {
      for (auto m0 : range_m0) {
         input.TanBeta = tanBeta;
         input.m0      = m0;

         spectrum_generator.run(oneset, input);

         const NMSSM<algorithm_type>& model = spectrum_generator.get_model();
         const NMSSM_physical& pole_masses = model.get_physical();
         const Problems<NMSSM_info::NUMBER_OF_PARTICLES>& problems
            = spectrum_generator.get_problems();
         const double higgs = pole_masses.Mhh(0);
         const bool error = problems.have_serious_problem();

         cout << "  "
              << std::setw(12) << std::left << input.TanBeta << ' '
              << std::setw(12) << std::left << input.m0 << ' '
              << std::setw(12) << std::left << higgs << ' '
              << std::setw(12) << std::left << error;
         if (error) {
            cout << "\t# " << problems;
         }
         cout << '\n';
      }
   }

   return 0;
}
