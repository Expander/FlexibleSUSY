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

// File generated at Wed 22 Jul 2015 18:12:15

#include "MSSMNoFV_onshell_susy_parameters.hpp"

#include <iostream>

namespace flexiblesusy {

MSSMNoFV_onshell_susy_parameters::MSSMNoFV_onshell_susy_parameters()
   : Beta_function()
   , Yd(Eigen::Matrix<double,3,3>::Zero()), Ye(Eigen::Matrix<double,3,3>::Zero(
   )), Yu(Eigen::Matrix<double,3,3>::Zero()), Mu(0), g1(0), g2(0), g3(0), vd(0)
   , vu(0)
{
   set_number_of_parameters(numberOfParameters);
}

MSSMNoFV_onshell_susy_parameters::MSSMNoFV_onshell_susy_parameters(
   double scale_, double loops_, double thresholds_
   , const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<double,3,3>& Ye_
   , const Eigen::Matrix<double,3,3>& Yu_, double Mu_, double g1_, double g2_,
   double g3_, double vd_, double vu_

)
   : Beta_function()
   , Yd(Yd_), Ye(Ye_), Yu(Yu_), Mu(Mu_), g1(g1_), g2(g2_), g3(g3_), vd(vd_), vu
   (vu_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

Eigen::ArrayXd MSSMNoFV_onshell_susy_parameters::beta() const
{
   Eigen::ArrayXd beta(numberOfParameters);
   beta.setZero();
   return beta;
}

MSSMNoFV_onshell_susy_parameters MSSMNoFV_onshell_susy_parameters::calc_beta() const
{
   return MSSMNoFV_onshell_susy_parameters();
}

void MSSMNoFV_onshell_susy_parameters::clear()
{
   reset();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Mu = 0.;
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   vd = 0.;
   vu = 0.;

}


Eigen::ArrayXd MSSMNoFV_onshell_susy_parameters::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = Yd(0,0);
   pars(1) = Yd(0,1);
   pars(2) = Yd(0,2);
   pars(3) = Yd(1,0);
   pars(4) = Yd(1,1);
   pars(5) = Yd(1,2);
   pars(6) = Yd(2,0);
   pars(7) = Yd(2,1);
   pars(8) = Yd(2,2);
   pars(9) = Ye(0,0);
   pars(10) = Ye(0,1);
   pars(11) = Ye(0,2);
   pars(12) = Ye(1,0);
   pars(13) = Ye(1,1);
   pars(14) = Ye(1,2);
   pars(15) = Ye(2,0);
   pars(16) = Ye(2,1);
   pars(17) = Ye(2,2);
   pars(18) = Yu(0,0);
   pars(19) = Yu(0,1);
   pars(20) = Yu(0,2);
   pars(21) = Yu(1,0);
   pars(22) = Yu(1,1);
   pars(23) = Yu(1,2);
   pars(24) = Yu(2,0);
   pars(25) = Yu(2,1);
   pars(26) = Yu(2,2);
   pars(27) = Mu;
   pars(28) = g1;
   pars(29) = g2;
   pars(30) = g3;
   pars(31) = vd;
   pars(32) = vu;


   return pars;
}

void MSSMNoFV_onshell_susy_parameters::print(std::ostream& ostr) const
{
   ostr << "susy parameters:\n";
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Mu = " << Mu << '\n';
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "vd = " << vd << '\n';
   ostr << "vu = " << vu << '\n';

}

void MSSMNoFV_onshell_susy_parameters::set(const Eigen::ArrayXd& pars)
{
   Yd(0,0) = pars(0);
   Yd(0,1) = pars(1);
   Yd(0,2) = pars(2);
   Yd(1,0) = pars(3);
   Yd(1,1) = pars(4);
   Yd(1,2) = pars(5);
   Yd(2,0) = pars(6);
   Yd(2,1) = pars(7);
   Yd(2,2) = pars(8);
   Ye(0,0) = pars(9);
   Ye(0,1) = pars(10);
   Ye(0,2) = pars(11);
   Ye(1,0) = pars(12);
   Ye(1,1) = pars(13);
   Ye(1,2) = pars(14);
   Ye(2,0) = pars(15);
   Ye(2,1) = pars(16);
   Ye(2,2) = pars(17);
   Yu(0,0) = pars(18);
   Yu(0,1) = pars(19);
   Yu(0,2) = pars(20);
   Yu(1,0) = pars(21);
   Yu(1,1) = pars(22);
   Yu(1,2) = pars(23);
   Yu(2,0) = pars(24);
   Yu(2,1) = pars(25);
   Yu(2,2) = pars(26);
   Mu = pars(27);
   g1 = pars(28);
   g2 = pars(29);
   g3 = pars(30);
   vd = pars(31);
   vu = pars(32);

}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFV_onshell_susy_parameters& susy_pars)
{
   susy_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
