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

#include "command_line_options.hpp"
#include "logger.hpp"
#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "error.hpp"
#include "wrappers.hpp"

#include "MSSMNoFVSLHA2_slha_io.hpp"
#include "MSSMNoFV_onshell.hpp"

using namespace flexiblesusy;

void print_usage(const char* program_name)
{
   std::cout <<
      "Usage: " << program_name << " [options]\n"
      "Options:\n"
      "  --slha-input-file=<source>\n      SLHA input source (file name or - for stdin)"
      "  --help,-h                         print this help message"
             << std::endl;
}

struct Gm2_cmd_line_options {
   std::string slha_input_source;
};

Gm2_cmd_line_options get_cmd_line_options(int argc, const char* argv[])
{
   Gm2_cmd_line_options options;

   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);

      if (Command_line_options::starts_with(option, "--slha-input-file=")) {
         options.slha_input_source = option.substr(18);
         continue;
      }

      if (option == "--help" || option == "-h") {
         print_usage(argv[0]);
         exit(EXIT_SUCCESS);
      }

      ERROR("Unrecognized command line option: " << option);
      exit(EXIT_FAILURE);
   }

   return options;
}

int main(int argc, const char* argv[])
{
   Gm2_cmd_line_options options(get_cmd_line_options(argc, argv));

   if (options.slha_input_source.empty()) {
      ERROR("No SLHA input source given!\n"
            "   Please provide one via the option --slha-input-file=");
      return EXIT_FAILURE;
   }

   MSSMNoFVSLHA2_mass_eigenstates model;
   MSSMNoFVSLHA2_slha_io slha_io;

   try {
      slha_io.read_from_source(options.slha_input_source);
      slha_io.fill(model);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   gm2os::MSSMNoFV_onshell osmodel(model);
   osmodel.calculate_DRbar_masses();

   try {
      osmodel.convert_to_onshell();
   } catch (const Error& e) {
      ERROR(e.what());
      return EXIT_FAILURE;
   }

   const double gm2_1l = gm2os::calculate_gm2_1loop(osmodel);
   const double gm2_2l = amu2LFSfapprox(osmodel)
                        + amuChipmPhotonic(osmodel) + amuChi0Photonic(osmodel);
   const double gm2_2l_tanb_approx =  + (tan_beta_cor(osmodel) - 1.) * gm2_1l;

   return 0;
}
