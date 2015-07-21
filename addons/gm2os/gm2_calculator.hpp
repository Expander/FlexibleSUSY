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

#ifndef GM2_CALCULATOR_H
#define GM2_CALCULATOR_H

#include "MSSMNoFV_mass_eigenstates.hpp"
#include <Eigen/Core>

namespace flexiblesusy {
namespace gm2 {

class Gm2_calculator : public MSSMNoFV_mass_eigenstates {
public:
   Gm2_calculator(const MSSMNoFV_mass_eigenstates&);
   virtual ~Gm2_calculator() {}

   Eigen::Array<double,2,1> get_MSmu() const {return get_MSm();}
   Eigen::Matrix<double,2,2> get_USmu() const {return get_ZM();}
   Eigen::Array<double,2,1> get_MStau() const {return get_MStau();}
   Eigen::Matrix<double,2,2> get_UStau() const {return get_ZTau();}
   Eigen::Array<double,2,1> get_MSbot() const {return get_MSb();}
   Eigen::Matrix<double,2,2> get_USbot() const {return get_ZB();}
   Eigen::Array<double,2,1> get_MStop() const {return get_MSt();}
   Eigen::Matrix<double,2,2> get_UStop() const {return get_ZT();}
   double get_MA0() const {return get_MAh(1);}
   double get_MUDIM() const {return get_scale();}
   double get_EL0() const {return EL0;}
   double get_gY() const {return sqrt(0.6) * get_g1();}
   double get_EL() const;
   double get_TB() const {return get_vu() / get_vd();}
   double get_MW() const {return get_physical().MVWm;}
   double get_MZ() const {return get_physical().MVZ;}
   double get_ME() const {return get_physical().MFe;}
   double get_MM() const {return get_physical().MFm;}
   double get_ML() const {return get_physical().MFtau;}
   double get_MU() const {return get_physical().MFu;}
   double get_MC() const {return get_physical().MFc;}
   double get_MT() const {return get_physical().MFt;}
   double get_MD() const {return get_physical().MFd;}
   double get_MS() const {return get_physical().MFs;}
   double get_MB() const {return get_physical().MFb;}
   Eigen::Matrix<double,3,3> get_Ae() const;
   Eigen::Matrix<double,3,3> get_Au() const;
   Eigen::Matrix<double,3,3> get_Ad() const;

   void set_MSmu(const Eigen::Array<double,2,1>& MSmu_new)
   {get_physical().MSm = MSmu_new;}
   void set_MStau(const Eigen::Array<double,2,1>& MStau_new)
   {get_physical().MStau = MStau_new;}
   void set_MSbot(const Eigen::Array<double,2,1>& MSbot_new)
   {get_physical().MSb = MSbot_new;}
   void set_MStop(const Eigen::Array<double,2,1>& MStop_new)
   {get_physical().MSt = MStop_new;}
   void set_MW(double MW_new)
   {get_physical().MVWm = MW_new;}
   void set_MZ(double MZ_new)
   {get_physical().MVZ = MZ_new;}
   void set_ME(double ME_new)
   {get_physical().MFe = ME_new;}
   void set_MM(double MM_new)
   {get_physical().MFm = MM_new;}
   void set_ML(double ML_new)
   {get_physical().MFtau = ML_new;}
   void set_MU(double MU_new)
   {get_physical().MFu = MU_new;}
   void set_MC(double MC_new)
   {get_physical().MFc = MC_new;}
   void set_MT(double MT_new)
   {get_physical().MFt = MT_new;}
   void set_MD(double MD_new)
   {get_physical().MFu = MD_new;}
   void set_MS(double MS_new)
   {get_physical().MFs = MS_new;}
   void set_MB(double MB_new)
   {get_physical().MFb = MB_new;}
   void set_MUDIM(double MUDIM_new)
   {set_scale(MUDIM_new);}
   void set_MA0(double MA0_new)
   {get_physical().MAh(1) = MA0_new;}

   void convert_to_onshell();

private:
   double EL0; ///< electromagnetic gauge coupling in the Thompson limit
};

} // namespace gm2
} // namespace flexiblesusy

#endif
