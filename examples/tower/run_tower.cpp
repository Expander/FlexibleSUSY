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

#include "MSSMD5O_input_parameters.hpp"
#include "MSSMD5O_slha_io.hpp"

#include "MSSMRHN_input_parameters.hpp"
#include "MSSMRHN_slha_io.hpp"

#include "MSSMD5O_MSSMRHN_spectrum_generator.hpp"

#include "error.hpp"
#include "spectrum_generator_settings.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

int main(int argc, char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_exit())
      return options.status();

   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   MSSMD5O_slha_io slha_io_1;
   MSSMRHN_slha_io slha_io_2;
   Spectrum_generator_settings spectrum_generator_settings;
   QedQcd qedqcd;
   MSSMD5O_input_parameters input_1;
   MSSMRHN_input_parameters input_2;

   if (!slha_input_file.empty()) {
      try {
         slha_io_1.read_from_file(slha_input_file);
         slha_io_1.fill(qedqcd);
         slha_io_1.fill(input_1);
         slha_io_1.fill(spectrum_generator_settings);
         slha_io_2.read_from_file(slha_input_file);
         slha_io_2.fill(input_2);
      } catch (const Error& error) {
         ERROR(error.what());
         return EXIT_FAILURE;
      }
   }
   qedqcd.toMz(); // run SM fermion masses to MZ

   MSSMD5O_MSSMRHN_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      spectrum_generator_settings.get(Spectrum_generator_settings::precision));
   spectrum_generator.set_max_iterations(
      spectrum_generator_settings.get(Spectrum_generator_settings::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      spectrum_generator_settings.get(Spectrum_generator_settings::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_parameter_output_scale_1(
      slha_io_1.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      spectrum_generator_settings.get(Spectrum_generator_settings::ewsb_loop_order));

   spectrum_generator.run(qedqcd, input_1, input_2);

   const MSSMD5O<algorithm_type>& model_1
      = spectrum_generator.get_model_1();
   const auto& problems = spectrum_generator.get_problems();

   // output
   slha_io_1.set_spinfo(problems);
   slha_io_1.set_sminputs(qedqcd);
   slha_io_1.set_input(input_1);
   if (!problems.have_problem())
      slha_io_1.set_spectrum(model_1);

   if (slha_output_file.empty()) {
      slha_io_1.write_to_stream(std::cout);
   } else {
      slha_io_1.write_to_file(slha_output_file);
   }

   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;
}
