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

#define MODELPARAMETER(p) model.get_##p()
#define LOCALPHYSICAL(p) physical.p
#define DEFINE_PARAMETER(p)                                            \
   typename std::remove_const<typename std::remove_reference<decltype(MODELPARAMETER(p))>::type>::type p;
#define DEFINE_PHYSICAL_PARAMETER(p) decltype(LOCALPHYSICAL(p)) p;

namespace {
   int sign(double x) { return x < 0 ? -1 : 1; }
   double signedsqr(double x) { return sign(x) * x * x; }
}

namespace flexiblesusy {
namespace gm2os {

double read_scale(const SLHA_io& slha_io)
{
   char const * const drbar_blocks[] =
      { "gauge", "Yu", "Yd", "Ye", "Te", "Td", "Tu", "HMIX", "MSQ2", "MSE2",
        "MSL2", "MSU2", "MSD2", "MSOFT" };

   double scale = 0.;

   for (unsigned i = 0; i < sizeof(drbar_blocks)/sizeof(*drbar_blocks); i++) {
      const double block_scale = slha_io.read_scale(drbar_blocks[i]);
      if (!is_zero(block_scale)) {
         if (!is_zero(scale) && !is_equal(scale, block_scale))
            WARNING("DR-bar parameters defined at different scales");
         scale = block_scale;
      }
   }

   return scale;
}

void fill_drbar_parameters(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   model.set_g1(slha_io.read_entry("gauge", 1) * Sqrt(5./3.));
   model.set_g2(slha_io.read_entry("gauge", 2));
   model.set_g3(slha_io.read_entry("gauge", 3));
   {
      DEFINE_PARAMETER(Yu);
      slha_io.read_block("Yu", Yu);
      model.set_Yu(Yu);
   }
   {
      DEFINE_PARAMETER(Yd);
      slha_io.read_block("Yd", Yd);
      model.set_Yd(Yd);
   }
   {
      DEFINE_PARAMETER(Ye);
      slha_io.read_block("Ye", Ye);
      model.set_Ye(Ye);
   }
   {
      DEFINE_PARAMETER(TYe);
      slha_io.read_block("Te", TYe);
      model.set_TYe(TYe);
   }
   {
      DEFINE_PARAMETER(TYd);
      slha_io.read_block("Td", TYd);
      model.set_TYd(TYd);
   }
   {
      DEFINE_PARAMETER(TYu);
      slha_io.read_block("Tu", TYu);
      model.set_TYu(TYu);
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
   {
      DEFINE_PHYSICAL_PARAMETER(ZH);
      slha_io.read_block("SCALARMIX", ZH);
      LOCALPHYSICAL(ZH) = ZH;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZA);
      slha_io.read_block("PSEUDOSCALARMIX", ZA);
      LOCALPHYSICAL(ZA) = ZA;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZP);
      slha_io.read_block("CHARGEMIX", ZP);
      LOCALPHYSICAL(ZP) = ZP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZN);
      slha_io.read_block("NMIX", ZN);
      LOCALPHYSICAL(ZN) = ZN;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UP);
      slha_io.read_block("VMIX", UP);
      LOCALPHYSICAL(UP) = UP;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(UM);
      slha_io.read_block("UMIX", UM);
      LOCALPHYSICAL(UM) = UM;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZD);
      slha_io.read_block("sdownmix", ZD);
      LOCALPHYSICAL(ZD) = ZD;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZS);
      slha_io.read_block("sstrmix", ZS);
      LOCALPHYSICAL(ZS) = ZS;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZB);
      slha_io.read_block("sbotmix", ZB);
      LOCALPHYSICAL(ZB) = ZB;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZU);
      slha_io.read_block("supmix", ZU);
      LOCALPHYSICAL(ZU) = ZU;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZC);
      slha_io.read_block("scharmmix", ZC);
      LOCALPHYSICAL(ZC) = ZC;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZT);
      slha_io.read_block("stopmix", ZT);
      LOCALPHYSICAL(ZT) = ZT;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZE);
      slha_io.read_block("selemix", ZE);
      LOCALPHYSICAL(ZE) = ZE;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZM);
      slha_io.read_block("smumix", ZM);
      LOCALPHYSICAL(ZM) = ZM;
   }
   {
      DEFINE_PHYSICAL_PARAMETER(ZTau);
      slha_io.read_block("staumix", ZTau);
      LOCALPHYSICAL(ZTau) = ZTau;
   }

   LOCALPHYSICAL(MVZ) = slha_io.read_entry("SMINPUTS", 4);
   LOCALPHYSICAL(MFd) = slha_io.read_entry("SMINPUTS", 21);
   LOCALPHYSICAL(MFs) = slha_io.read_entry("SMINPUTS", 23);
   LOCALPHYSICAL(MFb) = slha_io.read_entry("SMINPUTS", 5);
   LOCALPHYSICAL(MFu) = slha_io.read_entry("SMINPUTS", 22);
   LOCALPHYSICAL(MFc) = slha_io.read_entry("SMINPUTS", 24);
   LOCALPHYSICAL(MFt) = slha_io.read_entry("SMINPUTS", 6);
   LOCALPHYSICAL(MFe) = slha_io.read_entry("SMINPUTS", 11);
   LOCALPHYSICAL(MFm) = slha_io.read_entry("SMINPUTS", 13);
   LOCALPHYSICAL(MFtau) = slha_io.read_entry("SMINPUTS", 7);
   LOCALPHYSICAL(MSveL) = slha_io.read_entry("MASS", 1000012);
   LOCALPHYSICAL(MSvmL) = slha_io.read_entry("MASS", 1000014);
   LOCALPHYSICAL(MSvtL) = slha_io.read_entry("MASS", 1000016);
   LOCALPHYSICAL(MSd)(0) = slha_io.read_entry("MASS", 1000001);
   LOCALPHYSICAL(MSd)(1) = slha_io.read_entry("MASS", 2000001);
   LOCALPHYSICAL(MSu)(0) = slha_io.read_entry("MASS", 1000002);
   LOCALPHYSICAL(MSu)(1) = slha_io.read_entry("MASS", 2000002);
   LOCALPHYSICAL(MSe)(0) = slha_io.read_entry("MASS", 1000011);
   LOCALPHYSICAL(MSe)(1) = slha_io.read_entry("MASS", 2000011);
   LOCALPHYSICAL(MSm)(0) = slha_io.read_entry("MASS", 1000013);
   LOCALPHYSICAL(MSm)(1) = slha_io.read_entry("MASS", 2000013);
   LOCALPHYSICAL(MStau)(0) = slha_io.read_entry("MASS", 1000015);
   LOCALPHYSICAL(MStau)(1) = slha_io.read_entry("MASS", 2000015);
   LOCALPHYSICAL(MSs)(0) = slha_io.read_entry("MASS", 1000003);
   LOCALPHYSICAL(MSs)(1) = slha_io.read_entry("MASS", 2000003);
   LOCALPHYSICAL(MSc)(0) = slha_io.read_entry("MASS", 1000004);
   LOCALPHYSICAL(MSc)(1) = slha_io.read_entry("MASS", 2000004);
   LOCALPHYSICAL(MSb)(0) = slha_io.read_entry("MASS", 1000005);
   LOCALPHYSICAL(MSb)(1) = slha_io.read_entry("MASS", 2000005);
   LOCALPHYSICAL(MSt)(0) = slha_io.read_entry("MASS", 1000006);
   LOCALPHYSICAL(MSt)(1) = slha_io.read_entry("MASS", 2000006);
   LOCALPHYSICAL(Mhh)(0) = slha_io.read_entry("MASS", 25);
   LOCALPHYSICAL(Mhh)(1) = slha_io.read_entry("MASS", 35);
   LOCALPHYSICAL(MAh)(1) = slha_io.read_entry("MASS", 36);
   LOCALPHYSICAL(MHpm)(1) = slha_io.read_entry("MASS", 37);
   LOCALPHYSICAL(MChi)(0) = slha_io.read_entry("MASS", 1000022);
   LOCALPHYSICAL(MChi)(1) = slha_io.read_entry("MASS", 1000023);
   LOCALPHYSICAL(MChi)(2) = slha_io.read_entry("MASS", 1000025);
   LOCALPHYSICAL(MChi)(3) = slha_io.read_entry("MASS", 1000035);
   LOCALPHYSICAL(MCha)(0) = slha_io.read_entry("MASS", 1000024);
   LOCALPHYSICAL(MCha)(1) = slha_io.read_entry("MASS", 1000037);
   LOCALPHYSICAL(MVWm) = slha_io.read_entry("MASS", 24);
}

void fill_pole_masses(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   MSSMNoFVSLHA2_physical physical_hk;
   fill_physical(slha_io, physical_hk);
   physical_hk.convert_to_hk();
   model.get_physical() = physical_hk;
}

void fill(const SLHA_io& slha_io, MSSMNoFV_onshell& model)
{
   fill_pole_masses(slha_io, model);
   fill_drbar_parameters(slha_io, model);
}

} // namespace gm2os
} // namespace flexiblesusy
