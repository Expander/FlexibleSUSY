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

#include "MSSM_input_parameters.hpp"
#include "MSSM_slha_io.hpp"

#include "MSSMRHN_input_parameters.hpp"
#include "MSSMRHN_slha_io.hpp"

#include "MSSM_MSSMRHN_spectrum_generator.hpp"

#include "error.hpp"
#include "program_options.hpp"
#include "lowe.h"
#include "command_line_options.hpp"

#include <iostream>
#include <cstdlib>

int main(int argc, const char* argv[])
{
   using namespace flexiblesusy;
   using namespace softsusy;
   typedef Two_scale algorithm_type;

   Command_line_options options(argc, argv);
   if (options.must_exit())
      return options.status();

   const std::string slha_input_file(options.get_slha_input_file());
   const std::string slha_output_file(options.get_slha_output_file());
   MSSM_slha_io    slha_io_1;
   MSSMRHN_slha_io slha_io_2;
   Program_options program_options;
   QedQcd oneset;
   MSSM_input_parameters    input_1;
   MSSMRHN_input_parameters input_2;

   if (!slha_input_file.empty()) {
      try {
         slha_io_1.read_from_file(slha_input_file);
         slha_io_1.fill(oneset);
         slha_io_1.fill(input_1);
         slha_io_1.fill(program_options);
         slha_io_2.read_from_file(slha_input_file);
         slha_io_2.fill(input_2);
      } catch (const Error& error) {
         ERROR(error.what());
         return EXIT_FAILURE;
      }
   }
   oneset.toMz(); // run SM fermion masses to MZ

   MSSM_MSSMRHN_spectrum_generator<algorithm_type> spectrum_generator;
   spectrum_generator.set_precision_goal(
      program_options.get(Program_options::precision));
   spectrum_generator.set_max_iterations(
      program_options.get(Program_options::max_iterations));
   spectrum_generator.set_calculate_sm_masses(
      program_options.get(Program_options::calculate_sm_masses) >= 1.0);
   spectrum_generator.set_input_scale_2(
      slha_io_2.get_input_scale());
   spectrum_generator.set_parameter_output_scale_1(
      slha_io_1.get_parameter_output_scale());
   spectrum_generator.set_pole_mass_loop_order(
      program_options.get(Program_options::pole_mass_loop_order));
   spectrum_generator.set_ewsb_loop_order(
      program_options.get(Program_options::ewsb_loop_order));

   spectrum_generator.run(oneset, input_1, input_2);

   const MSSM<algorithm_type>& model_1
      = spectrum_generator.get_model_1();
   const Problems<MSSM_info::NUMBER_OF_PARTICLES>& problems
      = spectrum_generator.get_problems();
   const MSSMRHN<algorithm_type>& model_2
      = spectrum_generator.get_model_2();

   // output
   slha_io_1.set_spinfo(problems);
   slha_io_1.set_sminputs(oneset);
   slha_io_1.set_minpar(input_1);
   slha_io_1.set_extpar(input_1);
   if (!problems.have_serious_problem())
      slha_io_1.set_spectrum(model_1);

   if (slha_output_file.empty()) {
      slha_io_1.write_to_stream(std::cout);
   } else {
      slha_io_1.write_to_file(slha_output_file);
   }

   const int exit_code = spectrum_generator.get_exit_code();

   return exit_code;
}
