#ifndef TEST_VCMSSM_H
#define TEST_VCMSSM_H

#include "CMSSM_mass_eigenstates.hpp"
#include "VCMSSM_mass_eigenstates.hpp"

#include "ew_input.hpp"
#include "wrappers.hpp"

void setup_VCMSSM_const(flexiblesusy::VCMSSM_mass_eigenstates& m,
                         const flexiblesusy::VCMSSM_input_parameters& input)
{
   using namespace flexiblesusy;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5.0 * ALPHAMZ / (3.0 * (1.0 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = Sqrt(4.0 * Pi * alpha1);
   const double g2 = Sqrt(4.0 * Pi * alpha2);
   const double g3 = Sqrt(4.0 * Pi * ALPHASMZ);
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero + 100.;
   const double root2 = Sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * Sin(ArcTan(input.TBGuess));
   const double vd = vev * Cos(ArcTan(input.TBGuess));
   const double scale = Electroweak_constants::MZ;
   const double mu_guess = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> mm0(Eigen::Matrix<double,3,3>::Zero()),
      mm0b(Eigen::Matrix<double,3,3>::Zero()),
      Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero());

   mm0(0,0) = Sqr(m0);
   mm0(1,1) = Sqr(m0 + 1);
   mm0(2,2) = Sqr(m0 + 2);
   mm0b(0,0) = Sqr(m0 + 3);
   mm0b(1,1) = Sqr(m0 + 4);
   mm0b(2,2) = Sqr(m0 + 5);

   Yu(2,2) = 165.0   * root2 / vu;
   Yd(2,2) = 2.9     * root2 / vd;
   Ye(2,2) = 1.77699 * root2 / vd;

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_Mu(mu_guess);
   m.set_MassB(M12);
   m.set_MassWB(M12);
   m.set_MassG(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0b);
   m.set_mu2(mm0);
   m.set_me2(mm0b);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(-0.5 * Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_BMu((a0 + m0) * m.get_Mu());
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vMSSM(vev);
}

void setup_VCMSSM(flexiblesusy::VCMSSM_mass_eigenstates& m,
                  flexiblesusy::VCMSSM_input_parameters& input)
{
   input.m0 = 125.;
   input.m12 = 500.;
   input.Azero = 0.;
   input.SignMu = 1.;
   input.TBGuess = 5.0;

   setup_VCMSSM_const(m, input);
}

void match_CMSSM_to_VCMSSM(flexiblesusy::CMSSM_mass_eigenstates& cmssm,
                           const flexiblesusy::VCMSSM_mass_eigenstates& vcmssm)
{
   using namespace flexiblesusy;

   const VCMSSM_input_parameters vcmssm_input(vcmssm.get_input());
   CMSSM_input_parameters cmssm_input;
   cmssm_input.m0 = vcmssm_input.m0;
   cmssm_input.m12 = vcmssm_input.m12;
   cmssm_input.Azero = vcmssm_input.Azero;
   cmssm_input.SignMu = vcmssm_input.SignMu;
   cmssm.set_input_parameters(cmssm_input);

   cmssm.set_scale(vcmssm.get_scale());
   cmssm.set_loops(vcmssm.get_loops());
   cmssm.set_g1(vcmssm.get_g1());
   cmssm.set_g2(vcmssm.get_g2());
   cmssm.set_g3(vcmssm.get_g3());
   cmssm.set_Yu(vcmssm.get_Yu());
   cmssm.set_Yd(vcmssm.get_Yd());
   cmssm.set_Ye(vcmssm.get_Ye());
   cmssm.set_Mu(vcmssm.get_Mu());
   cmssm.set_MassB(vcmssm.get_MassB());
   cmssm.set_MassWB(vcmssm.get_MassWB());
   cmssm.set_MassG(vcmssm.get_MassG());
   cmssm.set_mq2(vcmssm.get_mq2());
   cmssm.set_ml2(vcmssm.get_ml2());
   cmssm.set_mu2(vcmssm.get_mu2());
   cmssm.set_md2(vcmssm.get_md2());
   cmssm.set_me2(vcmssm.get_me2());
   cmssm.set_mHd2(vcmssm.get_mHd2());
   cmssm.set_mHu2(vcmssm.get_mHu2());
   cmssm.set_TYu(vcmssm.get_TYu());
   cmssm.set_TYd(vcmssm.get_TYd());
   cmssm.set_TYe(vcmssm.get_TYe());
   cmssm.set_BMu(vcmssm.get_BMu());
   cmssm.set_vu(vcmssm.get_vu());
   cmssm.set_vd(vcmssm.get_vd());
}

#endif
