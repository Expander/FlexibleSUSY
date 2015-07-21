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

// File generated at Wed 5 Mar 2014 14:17:07

#include "logger.hpp"
#include "gm2_1loop.hpp"
#include "gm2_2loop.hpp"
#include "gm2_calculator.hpp"
#include "wrappers.hpp"

using namespace flexiblesusy;

/*
  * Setup model parameters.  The parameters are in agreement with gauge
  * coupling and soft parameter unification at the GUT scale.
  *
  * @param model model class to setup
*/

void setup_SPS1a(gm2::Gm2_calculator& model) {
   Eigen::Matrix<double,3,3> Yu;
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;
   double Mu;
   double TB;
   double EL;
   double g3;
   double MW;
   double MZ;
   double ME = 0.00051;
   double MM = 0.105658;
   double ML = 1.777;
   double MU = 0.04151;
   double MC = 1.5;
   double MT = 173.5;
   double MD = 0.0415;
   double MS = 0.15;
   double MB = 3.;
   Eigen::Matrix<double,3,3> Ae;
   Eigen::Matrix<double,3,3> Au;
   Eigen::Matrix<double,3,3> Ad;
   Eigen::Matrix<double,3,3> TYu;
   Eigen::Matrix<double,3,3> TYd;
   Eigen::Matrix<double,3,3> TYe;
   double BMu;
   Eigen::Matrix<double,3,3> mq2;
   Eigen::Matrix<double,3,3> ml2;
   double mHd2;
   double mHu2;
   Eigen::Matrix<double,3,3> md2;
   Eigen::Matrix<double,3,3> mu2;
   Eigen::Matrix<double,3,3> me2;
   double MassB;
   double MassWB;
   double MassG;
   double MUDIM;
   double MA0;

   // susy parameters
   Yu << 1.26136e-05,          0,          0,
                   0, 0.00667469,          0,
                   0,          0,   0.857849;

   Yd << 0.000242026,          0,          0,
                   0, 0.00529911,          0,
                   0,          0,   0.193602;

   Ye << 2.84161e-05,          0,          0,
                   0, 0.00587557,          0,
                   0,          0,    0.10199;

   Mu = 352.4;
   TB = 10.;
   MW = 80.385;
   MZ = 91.1876;
   EL = 0.308274;
   g3 = 1.06459;

   // soft parameters
   Ae << 1., 0., 0.,
         0., 1., 0.,
         0., 0., -254.2;

   Au << 1., 0., 0.,
         0., 1., 0.,
         0., 0., -510.;

   Ad << 1., 0., 0.,
         0., 1., 0.,
         0., 0., -772.7;

   TYu << -0.0144387,        0,        0, //f^uA^u
                   0, -7.64037,        0,
                   0,        0, -759.305;

   TYd << -0.336207,        0,        0, // analog
                  0, -7.36109,        0,
                  0,        0, -250.124;

   TYe << -0.00825134,        0,        0,
                    0, -1.70609,        0,
                    0,        0, -29.4466;

   BMu = 15338.71;

   mq2 << 291492.,           0,           0,//Massenterme squarks M_q^2
                    0, 291492.,           0,
                    0,           0,      245917.;

   ml2 << 38651.6,      0,      0,//analog
               0, 38651.6,      0,
               0,      0, 38337.6;

   mHd2 = 92436.9; // massenterme higgse m_1^2
   mHu2 = -380337; // m_2^2

   md2 << 269880.,      0,      0,
               0, 269880.,      0,
               0,      0, 267186.;

   mu2 << 272171.,      0,      0,
               0, 272171.,      0,
               0,      0, 180455.;

   me2 << 18550.4,       0,       0,
                0, 18550.4,       0,
                0,       0, 17849.;

   MassB = 99.1; //massenterme superpartner eichbosonen, M_1
   MassWB = 192.7;// M_2
   MassG = 1114.45;// M_3

   MUDIM = 454.7;// scale-parameter
   MA0 = 393.6;// one Higgsmass

   // set parameters
   model.set_TB(TB);
   model.set_MW(MW);
   model.set_MZ(MZ);
   model.set_EL(EL);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_Ae(Ae);
   model.set_Ad(Ad);
   model.set_Au(Au);
   model.set_ME(ME);
   model.set_MM(MM);
   model.set_ML(ML);
   model.set_MU(MU);
   model.set_MC(MC);
   model.set_MT(MT);
   model.set_MD(MD);
   model.set_MS(MS);
   model.set_MB(MB);
   model.set_Mu(Mu);
   model.set_g3(g3);
   model.set_TYu(TYu);
   model.set_TYd(TYd);
   model.set_TYe(TYe);
   model.set_BMu(BMu);
   model.set_mq2(mq2);
   model.set_ml2(ml2);
   model.set_mHd2(mHd2);
   model.set_mHu2(mHu2);
   model.set_md2(md2);
   model.set_mu2(mu2);
   model.set_me2(me2);
   model.set_MassB(MassB);
   model.set_MassWB(MassWB);
   model.set_MassG(MassG);
   model.set_MUDIM(MUDIM);
   model.set_MA0(MA0);

   // calculate tree-level masses
   model.calculate_DRbar_parameters();
}

int main()
{
   MSSMNoFV_mass_eigenstates model;
   gm2::Gm2_calculator calculator(model);
   setup_SPS1a(calculator);

   const double gm2_1l = gm2::calculate_gm2_1loop(calculator);
   const double gm2_2l = amu2LFSfapprox(calculator)
                        + amuChipmPhotonic(calculator) + amuChi0Photonic(calculator);
   const double gm2_2l_tanb_approx =  + (tan_beta_cor(calculator) - 1.) * gm2_1l;

   INFO("--------------------------------------");
   INFO("g-2 (1-loop) = " << gm2_1l);
   INFO("--------------------------------------"); 
   INFO("amuChi0 = " << gm2::amuChi0(calculator));
   INFO("amuChipm = " << gm2::amuChipm(calculator));
   INFO("--------------------------------------");
   INFO("amu1Lapprox = " << amu1Lapprox(calculator));
   INFO("--------------------------------------"); 
   INFO("amuWHnu = " << gm2::amuWHnu(calculator));
   INFO("amuBmuLmuR = " << gm2::amuBmuLmuR(calculator));
   INFO("amuBHmuL = " << gm2::amuBHmuL(calculator));
   INFO("amuWHmuL = " << gm2::amuWHmuL(calculator));
   INFO("amuBHmuR = " << gm2::amuBHmuR(calculator));
   INFO("--------------------------------------");
   INFO("----- g-2 (2-loop) - corrections -----");
   INFO("--------------------------------------");
   INFO("g-2 (2-loop) = " << gm2_2l);
   INFO("2Loop / 1Loop = " << 100. * gm2_2l / gm2_1l << " %");
   INFO("--------------------------------------");
   INFO("amu2LSFsapprox = " << amu2LFSfapprox(calculator));
   INFO("--------------------------------------"); 
   INFO("amuWHnu2L = " << amuWHnu2L(calculator));
   INFO("amuWHmuL2L = " << amuWHmuL2L(calculator));
   INFO("amuBHmuL2L = " << amuBHmuL2L(calculator));
   INFO("amuBHmuR2L = " << amuBHmuR2L(calculator));
   INFO("amuBmuLmuR2L = " << amuBmuLmuR2L(calculator));
   INFO("2L_FSfapprox / 1Loop = " << 100. * amu2LFSfapprox(calculator) / gm2_1l << " %");
   INFO("--------------------------------------"); 
   INFO("TanBetaCorrection) = " << gm2_2l_tanb_approx);
   INFO("2L_tanb / 1Loop = " << 100. * gm2_2l_tanb_approx / gm2_1l << " %");
   INFO("--------------------------------------");
   INFO("amu2LPhotonic = " << amuChipmPhotonic(calculator) + amuChi0Photonic(calculator));
   INFO("--------------------------------------");
   INFO("amuChipmPhotonic = " << amuChipmPhotonic(calculator));
   INFO("amuChi0Photonic = " << amuChi0Photonic(calculator));
   INFO("2L_Photonic / 1Loop = " << 100. * (amuChipmPhotonic(calculator)
                                    + amuChi0Photonic(calculator)) / gm2_1l << " %");
   INFO("--------------------------------------");
   INFO("amu2LaSferm = " << amua2LSferm(calculator));
   INFO("amua2LaCha = " << amua2LCha(calculator));
   INFO("--------------------------------------");

   return 0;
}
