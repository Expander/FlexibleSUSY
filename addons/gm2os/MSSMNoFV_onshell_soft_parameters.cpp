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

// File generated at Wed 22 Jul 2015 18:12:20

#include "MSSMNoFV_onshell_soft_parameters.hpp"

#include <iostream>

namespace flexiblesusy {

MSSMNoFV_onshell_soft_parameters::MSSMNoFV_onshell_soft_parameters()
   : MSSMNoFV_onshell_susy_parameters()
   , TYd(Eigen::Matrix<double,3,3>::Zero()), TYe(Eigen::Matrix<double,3,3>
   ::Zero()), TYu(Eigen::Matrix<double,3,3>::Zero()), BMu(0), mq2(Eigen::Matrix
   <double,3,3>::Zero()), ml2(Eigen::Matrix<double,3,3>::Zero()), mHd2(0), mHu2
   (0), md2(Eigen::Matrix<double,3,3>::Zero()), mu2(Eigen::Matrix<double,3,3>
   ::Zero()), me2(Eigen::Matrix<double,3,3>::Zero()), MassB(0), MassWB(0),
   MassG(0)

{
   set_number_of_parameters(numberOfParameters);
}

MSSMNoFV_onshell_soft_parameters::MSSMNoFV_onshell_soft_parameters(
   const MSSMNoFV_onshell_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, const Eigen::Matrix<double,3,3>& TYu_, double BMu_, const
   Eigen::Matrix<double,3,3>& mq2_, const Eigen::Matrix<double,3,3>& ml2_,
   double mHd2_, double mHu2_, const Eigen::Matrix<double,3,3>& md2_, const
   Eigen::Matrix<double,3,3>& mu2_, const Eigen::Matrix<double,3,3>& me2_,
   double MassB_, double MassWB_, double MassG_

)
   : MSSMNoFV_onshell_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TYu(TYu_), BMu(BMu_), mq2(mq2_), ml2(ml2_), mHd2(
   mHd2_), mHu2(mHu2_), md2(md2_), mu2(mu2_), me2(me2_), MassB(MassB_), MassWB(
   MassWB_), MassG(MassG_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd MSSMNoFV_onshell_soft_parameters::beta() const
{
   Eigen::ArrayXd beta(numberOfParameters);
   beta.setZero();
   return beta;
}

MSSMNoFV_onshell_soft_parameters MSSMNoFV_onshell_soft_parameters::calc_beta() const
{
   return MSSMNoFV_onshell_soft_parameters();
}

void MSSMNoFV_onshell_soft_parameters::clear()
{
   MSSMNoFV_onshell_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TYu = Eigen::Matrix<double,3,3>::Zero();
   BMu = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;
}

Eigen::ArrayXd MSSMNoFV_onshell_soft_parameters::get() const
{
   Eigen::ArrayXd pars(MSSMNoFV_onshell_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(33) = TYd(0,0);
   pars(34) = TYd(0,1);
   pars(35) = TYd(0,2);
   pars(36) = TYd(1,0);
   pars(37) = TYd(1,1);
   pars(38) = TYd(1,2);
   pars(39) = TYd(2,0);
   pars(40) = TYd(2,1);
   pars(41) = TYd(2,2);
   pars(42) = TYe(0,0);
   pars(43) = TYe(0,1);
   pars(44) = TYe(0,2);
   pars(45) = TYe(1,0);
   pars(46) = TYe(1,1);
   pars(47) = TYe(1,2);
   pars(48) = TYe(2,0);
   pars(49) = TYe(2,1);
   pars(50) = TYe(2,2);
   pars(51) = TYu(0,0);
   pars(52) = TYu(0,1);
   pars(53) = TYu(0,2);
   pars(54) = TYu(1,0);
   pars(55) = TYu(1,1);
   pars(56) = TYu(1,2);
   pars(57) = TYu(2,0);
   pars(58) = TYu(2,1);
   pars(59) = TYu(2,2);
   pars(60) = BMu;
   pars(61) = mq2(0,0);
   pars(62) = mq2(0,1);
   pars(63) = mq2(0,2);
   pars(64) = mq2(1,0);
   pars(65) = mq2(1,1);
   pars(66) = mq2(1,2);
   pars(67) = mq2(2,0);
   pars(68) = mq2(2,1);
   pars(69) = mq2(2,2);
   pars(70) = ml2(0,0);
   pars(71) = ml2(0,1);
   pars(72) = ml2(0,2);
   pars(73) = ml2(1,0);
   pars(74) = ml2(1,1);
   pars(75) = ml2(1,2);
   pars(76) = ml2(2,0);
   pars(77) = ml2(2,1);
   pars(78) = ml2(2,2);
   pars(79) = mHd2;
   pars(80) = mHu2;
   pars(81) = md2(0,0);
   pars(82) = md2(0,1);
   pars(83) = md2(0,2);
   pars(84) = md2(1,0);
   pars(85) = md2(1,1);
   pars(86) = md2(1,2);
   pars(87) = md2(2,0);
   pars(88) = md2(2,1);
   pars(89) = md2(2,2);
   pars(90) = mu2(0,0);
   pars(91) = mu2(0,1);
   pars(92) = mu2(0,2);
   pars(93) = mu2(1,0);
   pars(94) = mu2(1,1);
   pars(95) = mu2(1,2);
   pars(96) = mu2(2,0);
   pars(97) = mu2(2,1);
   pars(98) = mu2(2,2);
   pars(99) = me2(0,0);
   pars(100) = me2(0,1);
   pars(101) = me2(0,2);
   pars(102) = me2(1,0);
   pars(103) = me2(1,1);
   pars(104) = me2(1,2);
   pars(105) = me2(2,0);
   pars(106) = me2(2,1);
   pars(107) = me2(2,2);
   pars(108) = MassB;
   pars(109) = MassWB;
   pars(110) = MassG;

   return pars;
}

void MSSMNoFV_onshell_soft_parameters::print(std::ostream& ostr) const
{
   MSSMNoFV_onshell_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "BMu = " << BMu << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';
}

void MSSMNoFV_onshell_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   MSSMNoFV_onshell_susy_parameters::set(pars);

   TYd(0,0) = pars(33);
   TYd(0,1) = pars(34);
   TYd(0,2) = pars(35);
   TYd(1,0) = pars(36);
   TYd(1,1) = pars(37);
   TYd(1,2) = pars(38);
   TYd(2,0) = pars(39);
   TYd(2,1) = pars(40);
   TYd(2,2) = pars(41);
   TYe(0,0) = pars(42);
   TYe(0,1) = pars(43);
   TYe(0,2) = pars(44);
   TYe(1,0) = pars(45);
   TYe(1,1) = pars(46);
   TYe(1,2) = pars(47);
   TYe(2,0) = pars(48);
   TYe(2,1) = pars(49);
   TYe(2,2) = pars(50);
   TYu(0,0) = pars(51);
   TYu(0,1) = pars(52);
   TYu(0,2) = pars(53);
   TYu(1,0) = pars(54);
   TYu(1,1) = pars(55);
   TYu(1,2) = pars(56);
   TYu(2,0) = pars(57);
   TYu(2,1) = pars(58);
   TYu(2,2) = pars(59);
   BMu = pars(60);
   mq2(0,0) = pars(61);
   mq2(0,1) = pars(62);
   mq2(0,2) = pars(63);
   mq2(1,0) = pars(64);
   mq2(1,1) = pars(65);
   mq2(1,2) = pars(66);
   mq2(2,0) = pars(67);
   mq2(2,1) = pars(68);
   mq2(2,2) = pars(69);
   ml2(0,0) = pars(70);
   ml2(0,1) = pars(71);
   ml2(0,2) = pars(72);
   ml2(1,0) = pars(73);
   ml2(1,1) = pars(74);
   ml2(1,2) = pars(75);
   ml2(2,0) = pars(76);
   ml2(2,1) = pars(77);
   ml2(2,2) = pars(78);
   mHd2 = pars(79);
   mHu2 = pars(80);
   md2(0,0) = pars(81);
   md2(0,1) = pars(82);
   md2(0,2) = pars(83);
   md2(1,0) = pars(84);
   md2(1,1) = pars(85);
   md2(1,2) = pars(86);
   md2(2,0) = pars(87);
   md2(2,1) = pars(88);
   md2(2,2) = pars(89);
   mu2(0,0) = pars(90);
   mu2(0,1) = pars(91);
   mu2(0,2) = pars(92);
   mu2(1,0) = pars(93);
   mu2(1,1) = pars(94);
   mu2(1,2) = pars(95);
   mu2(2,0) = pars(96);
   mu2(2,1) = pars(97);
   mu2(2,2) = pars(98);
   me2(0,0) = pars(99);
   me2(0,1) = pars(100);
   me2(0,2) = pars(101);
   me2(1,0) = pars(102);
   me2(1,1) = pars(103);
   me2(1,2) = pars(104);
   me2(2,0) = pars(105);
   me2(2,1) = pars(106);
   me2(2,2) = pars(107);
   MassB = pars(108);
   MassWB = pars(109);
   MassG = pars(110);
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
