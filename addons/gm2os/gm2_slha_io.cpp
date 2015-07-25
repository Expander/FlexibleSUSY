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

#include "gm2_slha_io.hpp"
#include "slha_io.hpp"
#include "MSSMNoFV_onshell.hpp"

#include <cmath>
#include <limits>

#include <Eigen/Core>

namespace {
   int sign(double x) { return x < 0 ? -1 : 1; }
   double signedsqr(double x) { return sign(x) * x * x; }
}

namespace flexiblesusy {
namespace gm2os {

double read_scale(const SLHA_io& slha_io)
{
   char const * const drbar_blocks[] =
      { "Yu", "Yd", "Ye", "Ae", "Ad", "Au", "HMIX", "MSOFT" };

   double scale = 0.;

   for (unsigned i = 0; i < sizeof(drbar_blocks)/sizeof(*drbar_blocks); i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         scale = block_scale;
         break;
      }
   }

   if (is_zero(scale))
      ERROR("could not find renormalization scale in any SLHA block.");

   return scale;
}

void fill_drbar_parameters(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   {
      Eigen::Matrix<double,3,3> Ae(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AE", Ae);
      model.set_Ae(Ae);
   }
   {
      Eigen::Matrix<double,3,3> Au(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AU", Au);
      model.set_Au(Au);
   }
   {
      Eigen::Matrix<double,3,3> Ad(Eigen::Matrix<double,3,3>::Zero());
      slha_io.read_block("AD", Ad);
      model.set_Ad(Ad);
   }

   model.set_Mu(slha_io.read_entry("HMIX", 1));
   model.set_mHd2(slha_io.read_entry("MSOFT", 21));
   model.set_mHu2(slha_io.read_entry("MSOFT", 22));
   model.set_ml2(0, 0, signedsqr(slha_io.read_entry("MSOFT", 31)));
   model.set_ml2(1, 1, signedsqr(slha_io.read_entry("MSOFT", 32)));
   model.set_ml2(2, 2, signedsqr(slha_io.read_entry("MSOFT", 33)));
   model.set_me2(0, 0, signedsqr(slha_io.read_entry("MSOFT", 34)));
   model.set_me2(1, 1, signedsqr(slha_io.read_entry("MSOFT", 35)));
   model.set_me2(2, 2, signedsqr(slha_io.read_entry("MSOFT", 36)));
   model.set_mq2(0, 0, signedsqr(slha_io.read_entry("MSOFT", 41)));
   model.set_mq2(1, 1, signedsqr(slha_io.read_entry("MSOFT", 42)));
   model.set_mq2(2, 2, signedsqr(slha_io.read_entry("MSOFT", 43)));
   model.set_mu2(0, 0, signedsqr(slha_io.read_entry("MSOFT", 44)));
   model.set_mu2(1, 1, signedsqr(slha_io.read_entry("MSOFT", 45)));
   model.set_mu2(2, 2, signedsqr(slha_io.read_entry("MSOFT", 46)));
   model.set_md2(0, 0, signedsqr(slha_io.read_entry("MSOFT", 47)));
   model.set_md2(1, 1, signedsqr(slha_io.read_entry("MSOFT", 48)));
   model.set_md2(2, 2, signedsqr(slha_io.read_entry("MSOFT", 49)));
   model.set_MassB(slha_io.read_entry("MSOFT", 1));
   model.set_MassWB(slha_io.read_entry("MSOFT", 2));
   model.set_MassG(slha_io.read_entry("MSOFT", 3));

   const double tanb = slha_io.read_entry("HMIX", 2);
   const double MA2_drbar = slha_io.read_entry("HMIX", 4);
   const double MW = model.get_MW();
   const double MZ = model.get_MZ();
   const double cW = MW/MZ;
   const double sW = std::sqrt(1. - cW*cW);
   const double vev = 2. * model.get_MW() * sW / model.get_EL();
   const double sinb = tanb / std::sqrt(1 + tanb*tanb);
   const double cosb = 1.   / std::sqrt(1 + tanb*tanb);

   model.set_vd(vev * cosb);
   model.set_vu(vev * sinb);
   model.set_BMu(MA2_drbar * sinb * cosb);

   model.set_scale(read_scale(slha_io));
}

void fill_physical(const SLHA_io& slha_io, MSSMNoFVSLHA2_physical& physical)
{
   physical.MVWm = slha_io.read_entry("MASS", 24);
   physical.MVZ = slha_io.read_entry("SMINPUTS", 4);
   physical.MFd = slha_io.read_entry("SMINPUTS", 21);
   physical.MFs = slha_io.read_entry("SMINPUTS", 23);
   physical.MFb = slha_io.read_entry("SMINPUTS", 5);
   physical.MFu = slha_io.read_entry("SMINPUTS", 22);
   physical.MFc = slha_io.read_entry("SMINPUTS", 24);
   physical.MFt = slha_io.read_entry("SMINPUTS", 6);
   physical.MFe = slha_io.read_entry("SMINPUTS", 11);
   physical.MFm = slha_io.read_entry("SMINPUTS", 13);
   physical.MFtau = slha_io.read_entry("SMINPUTS", 7);
   physical.MSveL = slha_io.read_entry("MASS", 1000012);
   physical.MSvmL = slha_io.read_entry("MASS", 1000014);
   physical.MSvtL = slha_io.read_entry("MASS", 1000016);
   physical.MSd(0) = slha_io.read_entry("MASS", 1000001);
   physical.MSd(1) = slha_io.read_entry("MASS", 2000001);
   physical.MSu(0) = slha_io.read_entry("MASS", 1000002);
   physical.MSu(1) = slha_io.read_entry("MASS", 2000002);
   physical.MSe(0) = slha_io.read_entry("MASS", 1000011);
   physical.MSe(1) = slha_io.read_entry("MASS", 2000011);
   physical.MSm(0) = slha_io.read_entry("MASS", 1000013);
   physical.MSm(1) = slha_io.read_entry("MASS", 2000013);
   physical.MStau(0) = slha_io.read_entry("MASS", 1000015);
   physical.MStau(1) = slha_io.read_entry("MASS", 2000015);
   physical.MSs(0) = slha_io.read_entry("MASS", 1000003);
   physical.MSs(1) = slha_io.read_entry("MASS", 2000003);
   physical.MSc(0) = slha_io.read_entry("MASS", 1000004);
   physical.MSc(1) = slha_io.read_entry("MASS", 2000004);
   physical.MSb(0) = slha_io.read_entry("MASS", 1000005);
   physical.MSb(1) = slha_io.read_entry("MASS", 2000005);
   physical.MSt(0) = slha_io.read_entry("MASS", 1000006);
   physical.MSt(1) = slha_io.read_entry("MASS", 2000006);
   physical.Mhh(0) = slha_io.read_entry("MASS", 25);
   physical.Mhh(1) = slha_io.read_entry("MASS", 35);
   physical.MAh(1) = slha_io.read_entry("MASS", 36);
   physical.MHpm(1) = slha_io.read_entry("MASS", 37);
   physical.MChi(0) = slha_io.read_entry("MASS", 1000022);
   physical.MChi(1) = slha_io.read_entry("MASS", 1000023);
   physical.MChi(2) = slha_io.read_entry("MASS", 1000025);
   physical.MChi(3) = slha_io.read_entry("MASS", 1000035);
   physical.MCha(0) = slha_io.read_entry("MASS", 1000024);
   physical.MCha(1) = slha_io.read_entry("MASS", 1000037);
}

void fill_pole_masses(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   MSSMNoFVSLHA2_physical physical_hk;
   fill_physical(slha_io, physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

void fill_gm2_specific(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   const double alpha_MZ = std::abs(slha_io.read_entry("FlexibleSUSYGM2", 1));
   const double alpha_thompson = std::abs(slha_io.read_entry("FlexibleSUSYGM2", 2));

   if (alpha_MZ > std::numeric_limits<double>::epsilon())
      model.set_alpha_MZ(alpha_MZ);

   if (alpha_thompson > std::numeric_limits<double>::epsilon())
      model.set_alpha_thompson(alpha_thompson);
}

void fill(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_pole_masses(slha_io, model);
   fill_drbar_parameters(slha_io, model);
   fill_gm2_specific(slha_io, model);
}

} // namespace gm2os
} // namespace flexiblesusy
