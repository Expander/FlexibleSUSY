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
   const double gm2_2l = gm2os::amu2LFSfapprox(osmodel)
      + gm2os::amuChipmPhotonic(osmodel)
      + gm2os::amuChi0Photonic(osmodel);
   const double gm2_2l_tanb_approx =  + (gm2os::tan_beta_cor(osmodel) - 1.) * gm2_1l;

   INFO(
      "--------------------------------------\n"
      " on-shell parameters \n"
      "--------------------------------------\n"
      "cos(TheatW) = " << (osmodel.get_MW() / osmodel.get_MZ()) << '\n' <<
      "BMu         = " << osmodel.get_BMu() << '\n' <<
      "SM vev      = " << Sqrt(Sqr(osmodel.get_vu()) + Sqr(osmodel.get_vd())) << '\n' <<
      "yu          = " << osmodel.get_Yu().diagonal().transpose() << '\n' <<
      "yd          = " << osmodel.get_Yd().diagonal().transpose() << '\n' <<
      "ye          = " << osmodel.get_Ye().diagonal().transpose() << '\n' <<
      "Mu          = " << osmodel.get_Mu() << '\n' <<
      "M1          = " << osmodel.get_MassB() << '\n' <<
      "M2          = " << osmodel.get_MassWB() << '\n' <<
      "ml2(2,2)    = " << osmodel.get_ml2(1,1) << '\n' <<
      "me2(2,2)    = " << osmodel.get_me2(1,1) << '\n'
      );

   INFO(
      "--------------------------------------\n"
      "g-2 (1-loop) = " << gm2_1l << '\n' <<
      "--------------------------------------\n"
      "amuChi0 = " << gm2os::amuChi0(osmodel) << '\n' <<
      "amuChipm = " << gm2os::amuChipm(osmodel) << '\n' <<
      "--------------------------------------\n"
      "amu1Lapprox = " << amu1Lapprox(osmodel) << '\n' <<
      "--------------------------------------\n"
      "amuWHnu = " << gm2os::amuWHnu(osmodel) << '\n' <<
      "amuBmuLmuR = " << gm2os::amuBmuLmuR(osmodel) << '\n' <<
      "amuBHmuL = " << gm2os::amuBHmuL(osmodel) << '\n' <<
      "amuWHmuL = " << gm2os::amuWHmuL(osmodel) << '\n' <<
      "amuBHmuR = " << gm2os::amuBHmuR(osmodel) << "\n\n" <<
      "--------------------------------------\n"
      "----- g-2 (2-loop) - corrections -----\n"
      "--------------------------------------\n"
      "g-2 (2-loop) = " << gm2_2l << '\n' <<
      "2Loop / 1Loop = " << 100. * gm2_2l / gm2_1l << " %\n"
      "--------------------------------------\n"
      "amu2LSFsapprox = " << gm2os::amu2LFSfapprox(osmodel) << '\n' <<
      "--------------------------------------\n"
      "amuWHnu2L = " << gm2os::amuWHnu2L(osmodel) << '\n' <<
      "amuWHmuL2L = " << gm2os::amuWHmuL2L(osmodel) << '\n' <<
      "amuBHmuL2L = " << gm2os::amuBHmuL2L(osmodel) << '\n' <<
      "amuBHmuR2L = " << gm2os::amuBHmuR2L(osmodel) << '\n' <<
      "amuBmuLmuR2L = " << gm2os::amuBmuLmuR2L(osmodel) << '\n' <<
      "2L_FSfapprox / 1Loop = " <<
      100. * gm2os::amu2LFSfapprox(osmodel) / gm2_1l << " %\n"
      "--------------------------------------\n"
      "TanBetaCorrection) = " << gm2_2l_tanb_approx << '\n' <<
      "2L_tanb / 1Loop = " << (100. * gm2_2l_tanb_approx / gm2_1l) << " %\n"
      "--------------------------------------\n"
      "amu2LPhotonic = " <<
      (gm2os::amuChipmPhotonic(osmodel)
       + gm2os::amuChi0Photonic(osmodel)) << '\n' <<
      "--------------------------------------\n"
      "amuChipmPhotonic = " << gm2os::amuChipmPhotonic(osmodel) << '\n' <<
      "amuChi0Photonic = " << gm2os::amuChi0Photonic(osmodel) << '\n' <<
      "2L_Photonic / 1Loop = " <<
      100. * (amuChipmPhotonic(osmodel)
              + gm2os::amuChi0Photonic(osmodel)) / gm2_1l << " %\n"
      "--------------------------------------\n"
      "amu2LaSferm = " << gm2os::amua2LSferm(osmodel) << '\n' <<
      "amua2LaCha = " << gm2os::amua2LCha(osmodel) << '\n' <<
      "--------------------------------------"
      );

   return 0;
}
