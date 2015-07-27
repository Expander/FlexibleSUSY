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

#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "error.hpp"
#include "config.h"

#include "gm2_slha_io.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <iostream>

#define ERROR(message) std::cerr << "Error: " << message << '\n';

using namespace flexiblesusy;

void print_usage(const char* program_name)
{
   std::cout <<
      "Usage: " << program_name << " [options]\n"
      "Options:\n"
      "  --slha-input-file=<source>      SLHA input source (file name or - for stdin)\n"
      "  --help,-h                       print this help message\n"
      "  --version,-v                    print version number"
             << std::endl;
}

struct Gm2_cmd_line_options {
   std::string slha_input_source;

   static bool starts_with(const std::string& str, const std::string& prefix) {
      return !str.compare(0, prefix.size(), prefix);
   }
};

Gm2_cmd_line_options get_cmd_line_options(int argc, const char* argv[])
{
   Gm2_cmd_line_options options;

   for (int i = 1; i < argc; ++i) {
      const std::string option(argv[i]);

      if (Gm2_cmd_line_options::starts_with(option, "--slha-input-file=")) {
         options.slha_input_source = option.substr(18);
         continue;
      }

      if (option == "--help" || option == "-h") {
         print_usage(argv[0]);
         exit(EXIT_SUCCESS);
      }

      if (option == "--version" || option == "-v") {
         std::cout << GM2CALC_VERSION << '\n';
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

   gm2calc::MSSMNoFV_onshell osmodel;
   gm2calc::GM2_slha_io slha_io;

   try {
      slha_io.read_from_source(options.slha_input_source);
      fill(slha_io, osmodel);
   } catch (const Error& error) {
      ERROR(error.what());
      return EXIT_FAILURE;
   }

   try {
      osmodel.convert_to_onshell();
   } catch (const Error& e) {
      ERROR(e.what());
      return EXIT_FAILURE;
   }

   const double gm2_1l = gm2calc::calculate_amu_1loop_non_tan_beta_resummed(osmodel);
   const double gm2_1l_TBresummed = gm2calc::calculate_amu_1loop(osmodel);
   const double gm2_2l_TBresummed = gm2calc::amu2LFSfapprox(osmodel)
      + gm2calc::amuChipmPhotonic(osmodel)
      + gm2calc::amuChi0Photonic(osmodel);
   const double gm2_2l_tanb_approx =  + (gm2calc::tan_beta_cor(osmodel) - 1.) * gm2_1l;

   const double gm2_best = gm2_1l_TBresummed + gm2_2l_TBresummed;

   std::cout << osmodel << '\n';

   std::cout <<
      "--------------------------------------\n"
      "g-2 (1-loop + 2-loop best) = " << gm2_best << '\n' <<
      "--------------------------------------\n"
      "g-2 (1-loop strict) = " << gm2_1l << '\n' <<
      "--------------------------------------\n"
      "g-2 (1-loop TB resummed) = " << gm2_1l_TBresummed << '\n' <<
      "--------------------------------------\n"
      "amuChi0 (TB resummed) = " << gm2calc::amuChi0(osmodel) << '\n' <<
      "amuChipm (TB resummed) = " << gm2calc::amuChipm(osmodel) << '\n' <<
      "--------------------------------------\n"
      "amu1Lapprox = " << amu1Lapprox(osmodel) << '\n' <<
      "--------------------------------------\n"
      "amuWHnu = " << gm2calc::amuWHnu(osmodel) << '\n' <<
      "amuBmuLmuR = " << gm2calc::amuBmuLmuR(osmodel) << '\n' <<
      "amuBHmuL = " << gm2calc::amuBHmuL(osmodel) << '\n' <<
      "amuWHmuL = " << gm2calc::amuWHmuL(osmodel) << '\n' <<
      "amuBHmuR = " << gm2calc::amuBHmuR(osmodel) << "\n\n" <<
      "--------------------------------------\n"
      "----- g-2 (2-loop) - corrections -----\n"
      "--------------------------------------\n"
      "g-2 (2-loop (TB resummed)) = " << gm2_2l_TBresummed << '\n' <<
      "2Loop / 1Loop = " << 100. * gm2_2l_TBresummed / gm2_1l_TBresummed << " %\n"
      "--------------------------------------\n"
      "amu2LSFsapprox = " << gm2calc::amu2LFSfapprox(osmodel) << '\n' <<
      "--------------------------------------\n"
      "amuWHnu2L = " << gm2calc::amuWHnu2L(osmodel) << '\n' <<
      "amuWHmuL2L = " << gm2calc::amuWHmuL2L(osmodel) << '\n' <<
      "amuBHmuL2L = " << gm2calc::amuBHmuL2L(osmodel) << '\n' <<
      "amuBHmuR2L = " << gm2calc::amuBHmuR2L(osmodel) << '\n' <<
      "amuBmuLmuR2L = " << gm2calc::amuBmuLmuR2L(osmodel) << '\n' <<
      "2L_FSfapprox / 1Loop = " <<
      100. * gm2calc::amu2LFSfapprox(osmodel) / gm2_1l << " %\n"
      "--------------------------------------\n"
      "TanBetaCorrection) = " << gm2_2l_tanb_approx << '\n' <<
      "2L_tanb / 1Loop = " << (100. * gm2_2l_tanb_approx / gm2_1l) << " %\n"
      "--------------------------------------\n"
      "amu2LPhotonic = " <<
      (gm2calc::amuChipmPhotonic(osmodel)
       + gm2calc::amuChi0Photonic(osmodel)) << '\n' <<
      "--------------------------------------\n"
      "amuChipmPhotonic = " << gm2calc::amuChipmPhotonic(osmodel) << '\n' <<
      "amuChi0Photonic = " << gm2calc::amuChi0Photonic(osmodel) << '\n' <<
      "2L_Photonic / 1Loop = " <<
      100. * (amuChipmPhotonic(osmodel)
              + gm2calc::amuChi0Photonic(osmodel)) / gm2_1l << " %\n"
      "--------------------------------------\n"
      "amu2LaSferm = " << gm2calc::amua2LSferm(osmodel) << '\n' <<
      "amua2LaCha = " << gm2calc::amua2LCha(osmodel) << '\n' <<
      "--------------------------------------\n"
      ;

   return 0;
}
