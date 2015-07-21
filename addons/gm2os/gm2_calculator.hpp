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

   Eigen::Array<double,2,1> get_MSmu() const {return MSmu;}
   Eigen::Matrix<double,2,2> get_USmu() const {return USmu;}
   Eigen::Array<double,2,1> get_MStau() const {return MStau;}
   Eigen::Matrix<double,2,2> get_UStau() const {return UStau;}
   Eigen::Array<double,2,1> get_MSbot() const {return MSbot;}
   Eigen::Matrix<double,2,2> get_USbot() const {return USbot;}
   Eigen::Array<double,2,1> get_MStop() const {return MStop;}
   Eigen::Matrix<double,2,2> get_UStop() const {return UStop;}
   double get_MA0() const {return MA0;}
   double get_MUDIM() const {return MUDIM;}
   double get_EL0() const {return EL0;}
   double get_gY() const {return sqrt(3. / 5.) * get_g1();}
   double get_EL() const {return EL;}
   double get_TB() const {return TB;}
   double get_MW() const {return MW;}
   double get_MZ() const {return MZ;}
   double get_ME() const {return ME;}
   double get_MM() const {return MM;}
   double get_ML() const {return ML;}
   double get_MU() const {return MU;}
   double get_MC() const {return MC;}
   double get_MT() const {return MT;}
   double get_MD() const {return MD;}
   double get_MS() const {return MS;}
   double get_MB() const {return MB;}
   Eigen::Matrix<double,3,3> get_Ae() const {return Ae;}
   Eigen::Matrix<double,3,3> get_Au() const {return Au;}
   Eigen::Matrix<double,3,3> get_Ad() const {return Ad;}

   void set_MSmu(const Eigen::Array<double,2,1>& MSmu_new)
   {MSmu = MSmu_new;}
   void set_USmu(const Eigen::Matrix<double,2,2>& USmu_new)
   {USmu = USmu_new;}
   void set_MStau(const Eigen::Array<double,2,1>& MStau_new)
   {MStau = MStau_new;}
   void set_UStau(const Eigen::Matrix<double,2,2>& UStau_new)
   {UStau = UStau_new;}
   void set_MSbot(const Eigen::Array<double,2,1>& MSbot_new)
   {MSbot = MSbot_new;}
   void set_USbot(const Eigen::Matrix<double,2,2>& USbot_new)
   {USbot = USbot_new;}
   void set_MStop(const Eigen::Array<double,2,1>& MStop_new)
   {MStop = MStop_new;}
   void set_UStop(const Eigen::Matrix<double,2,2>& UStop_new)
   {UStop = UStop_new;}
   void set_TB(double TB_new)
   {TB = TB_new;}
   void set_MW(double MW_new)
   {MW = MW_new;}
   void set_MZ(double MZ_new)
   {MZ = MZ_new;}
   void set_EL(double EL_new)
   {EL = EL_new;}
   void set_ME(double ME_new)
   {ME = ME_new;}
   void set_MM(double MM_new)
   {MM = MM_new;}
   void set_ML(double ML_new)
   {ML = ML_new;}
   void set_MU(double MU_new)
   {MU = MU_new;}
   void set_MC(double MC_new)
   {MC = MC_new;}
   void set_MT(double MT_new)
   {MT = MT_new;}
   void set_MD(double MD_new)
   {MD = MD_new;}
   void set_MS(double MS_new)
   {MS = MS_new;}
   void set_MB(double MB_new)
   {MB = MB_new;}
   void set_MUDIM(double MUDIM_new)
   {MUDIM = MUDIM_new;}
   void set_MA0(double MA0_new)
   {MA0 = MA0_new;}
   void set_Ae(const Eigen::Matrix<double,3,3>& Ae_new)
   {Ae = Ae_new;}
   void set_Au(const Eigen::Matrix<double,3,3>& Au_new)
   {Au = Au_new;}
   void set_Ad(const Eigen::Matrix<double,3,3>& Ad_new)
   {Ad = Ad_new;}
   void convert_parameters();
   void calculate_DRbar_parameters();
   void convert_parameters_reverse();

private:
   Eigen::Matrix<double,3,3> Ae;
   Eigen::Matrix<double,3,3> Au;
   Eigen::Matrix<double,3,3> Ad;
   Eigen::Array<double,2,1> MSmu;
   Eigen::Matrix<double,2,2> USmu;
   Eigen::Array<double,2,1> MStau;
   Eigen::Matrix<double,2,2> UStau;
   Eigen::Array<double,2,1> MSbot;
   Eigen::Matrix<double,2,2> USbot;
   Eigen::Array<double,2,1> MStop;
   Eigen::Matrix<double,2,2> UStop;
   double MW;
   double MZ;
   double TB;
   double EL;
   double ME;
   double MM;
   double ML;
   double MU;
   double MC;
   double MT;
   double MD;
   double MS;
   double MB;
   double gY;
   double EL0;
   double MUDIM;
   double MA0;
};

} // namespace gm2
} // namespace flexiblesusy

#endif
