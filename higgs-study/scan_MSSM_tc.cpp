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

// File generated at Sun 27 Apr 2014 19:19:50

#include "MSSMNoFV_input_parameters.hpp"
#include "MSSMNoFV_spectrum_generator.hpp"

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

   MSSMNoFV_input_parameters input;
   QedQcd oneset;
   oneset.toMz();

   MSSMNoFV_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(1.0e-3);
   spectrum_generator.set_max_iterations(0);         // 0 == automatic
   spectrum_generator.set_calculate_sm_masses(0);    // 0 == no
   spectrum_generator.set_parameter_output_scale(0); // 0 == susy scale

   const unsigned loop_order = 2;
   spectrum_generator.set_pole_mass_loop_order(loop_order);
   spectrum_generator.set_ewsb_loop_order(loop_order);

   const std::vector<double> range_TanBeta(
      float_range(options.tanb_start, options.tanb_stop, options.tanb_npoints));
   const std::vector<double> range_tc(
      float_range(options.tc_start, options.tc_stop, options.tc_npoints));

   input.m0 = 5000.;
   input.m12 = 5000.;
   input.TanBeta = 10.;
   input.Azero = 5000.;
   input.SignMu = 1;

   cout << "# MSSMNoFV with "
        << "Azero = " << input.Azero
        << ", m0 = " << input.m0
        << ", m12 = " << input.m12
        << ", SignMu = " << input.SignMu
        << '\n';

   cout << "# scan range: "
        << "TanBeta = " << options.tanb_start
        << " ... " << options.tanb_stop
        << " with " << options.tanb_npoints << " points"
        << '\n';

   cout << "# "
        << std::setw(12) << std::left << "TanBeta" << ' '
        << std::setw(12) << std::left << "tc" << ' '
        << std::setw(12) << std::left << "Mhh(1)/GeV" << ' '
        << std::setw(12) << std::left << "error" << ' '
        << std::setw(12) << std::left << "MX/GeV" << ' '
        << '\n';

   for (auto tanBeta : range_TanBeta) {
      for (auto tc : range_tc) {
         input.TanBeta = tanBeta;
         input.tc      = tc;

         spectrum_generator.run(oneset, input);

         const MSSMNoFV<algorithm_type>& model = spectrum_generator.get_model();
         const MSSMNoFV_physical& pole_masses = model.get_physical();
         const Problems<MSSMNoFV_info::NUMBER_OF_PARTICLES>& problems
            = spectrum_generator.get_problems();
         const double higgs = pole_masses.Mhh(0);
         const bool error = problems.have_serious_problem();

         cout << "  "
              << std::setw(12) << std::left << input.TanBeta << ' '
              << std::setw(12) << std::left << input.tc << ' '
              << std::setw(12) << std::left << higgs << ' '
              << std::setw(12) << std::left << error << ' '
              << std::setw(12) << std::left << spectrum_generator.get_high_scale()
            ;
         if (error) {
            cout << "\t# " << problems;
         }
         cout << '\n';
      }
   }

   return 0;
}
