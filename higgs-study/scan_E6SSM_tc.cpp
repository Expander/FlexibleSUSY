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

#include "phdE6SSM_input_parameters.hpp"
#include "phdE6SSM_spectrum_generator.hpp"

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

   phdE6SSM_input_parameters input;
   QedQcd oneset;
   oneset.toMz();

   phdE6SSM_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(2.0e-3);
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
   input.LambdaInput = 0.1;
   input.KappaInput = 0.1;
   input.muPrimeInput = 10000.;
   input.BmuPrimeInput = 10000.;
   input.vSInput = 10000.;
   input.tc = 0.;

   cout << "# phdE6SSM with "
        << "m0 = " << input.m0
        << ", m12 = " << input.m12
        << ", TanBeta = " << input.TanBeta
        << ", Azero = " << input.Azero
        << ", LambdaInput = " << input.LambdaInput
        << ", KappaInput = " << input.KappaInput
        << ", muPrimeInput = " << input.muPrimeInput
        << ", BmuPrimeInput = " << input.BmuPrimeInput
        << ", vSInput = " << input.vSInput
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
        << std::setw(12) << std::left << "MSu(1)/GeV" << ' '
        << std::setw(12) << std::left << "MSu(6)/GeV" << ' '
        << '\n';

   for (auto tanBeta : range_TanBeta) {
      for (auto tc : range_tc) {
         input.TanBeta = tanBeta;
         input.tc      = tc;

         spectrum_generator.run(oneset, input);

         const phdE6SSM<algorithm_type>& model = spectrum_generator.get_model();
         const phdE6SSM_physical& pole_masses = model.get_physical();
         const Problems<phdE6SSM_info::NUMBER_OF_PARTICLES>& problems
            = spectrum_generator.get_problems();
         const double higgs = pole_masses.Mhh(0);
         const double sup1 = pole_masses.MSu(0);
         const double sup6 = pole_masses.MSu(5);
         const bool error = problems.have_serious_problem();

         cout << "  "
              << std::setw(12) << std::left << input.TanBeta << ' '
              << std::setw(12) << std::left << input.tc << ' '
              << std::setw(12) << std::left << higgs << ' '
              << std::setw(12) << std::left << error << ' '
              << std::setw(12) << std::left << spectrum_generator.get_high_scale() << ' '
              << std::setw(12) << std::left << sup1 << ' '
              << std::setw(12) << std::left << sup6 << ' '
            ;
         if (error) {
            cout << "\t# " << problems;
         }
         cout << '\n';
      }
   }

   return 0;
}
