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

#include "gm2_2loop.hpp"
#include "gm2_1loop.hpp"
#include "MSSMNoFV_onshell.hpp"
#include "ffunctions.hpp"

#include <complex>
#include <cmath>

namespace flexiblesusy {
namespace gm2os {

//tanbeta2 corrections

double tan_beta_cor(const MSSMNoFV_onshell& model)
{
   const double delta_mu = delta_mu_correction(model);

   return 1. / (1. + delta_mu);
}

double delta_down_lepton_correction(const MSSMNoFV_onshell& model, int gen) {
   double delta_lep;
   const double mu = model.get_Mu();
   const double TB = model.get_TB();
   const double g2 = model.get_g2();
   const double gY = model.get_gY();
   const double M1 = model.get_MassB();
   const double M2 = model.get_MassWB();
   const double MW = model.get_MW();
   const double MZ = model.get_MZ();
   const double SW = sqrt(1. - sqr(MW / MZ));

   double m1 = ( sqrt(0.5 * (sqr(M2) + sqr(mu) + 2. * sqr(MW)
               - sqrt(sqr(sqr(M2) + sqr(mu) + 2. * sqr(MW)) - sqr(2. * M2 * mu)))) );
   double m2 = ( sqrt(0.5 * (sqr(M2) + sqr(mu) + 2. * sqr(MW)
               + sqrt(sqr(sqr(M2) + sqr(mu) + 2. * sqr(MW)) - sqr(2. * M2 * mu)))) );
   double m_sneu_lep = sqrt(model.get_ml2(gen, gen) - 0.5 * sqr(MZ));
   double m_slep_L = sqrt(model.get_ml2(gen, gen) - sqr(MZ) * (sqr(SW) - 0.5));
   double m_slep_R = sqrt(model.get_me2(gen, gen) + sqr(MZ * SW));

   delta_lep = ( - mu * TB * oneOver16PiSqr
            * (sqr(g2) * M2 * (Iabc(m1, m2, m_sneu_lep) + 0.5 * Iabc(m1, m2, m_slep_L))
            + sqr(gY) * M1 * (Iabc(mu, M1, m_slep_R) - 0.5 * Iabc(mu, M1, m_slep_L)
                              - Iabc(M1, m_slep_L, m_slep_R))) );

   return delta_lep;
}

double delta_mu_correction(const MSSMNoFV_onshell& model)
{
   return delta_down_lepton_correction(model, 1);
}

double delta_tau_correction(const MSSMNoFV_onshell& model)
{
   return delta_down_lepton_correction(model, 2);
}

/**
 * Returns the \f$\Delta_b\f$ corrections from arxiv:0901.2065,
 * Eq. (103) and Eqs. (31)-(35)
 */
double delta_bottom_correction(const MSSMNoFV_onshell& model)
{
   const double TB = model.get_TB();
   const double At = model.get_Au(2,2);
   const double yt = model.get_Yu(2,2);
   const double gY = model.get_gY();
   const double g2 = model.get_g2();
   const double alpha_S = sqr(model.get_g3()) / (4*M_PI);
   const double mu = model.get_Mu();
   const double M1 = model.get_MassB();
   const double M2 = model.get_MassWB();
   const double M3 = model.get_MassG();
   const double mstL2 = std::abs(model.get_mq2(2,2));
   const double mstR2 = std::abs(model.get_mu2(2,2));
   const double msbL2 = std::abs(model.get_mq2(1,1));
   const double msbR2 = std::abs(model.get_md2(1,1));

   const double eps_0 =
      - 2. * alpha_S / (3*M_PI) * std::conj(mu) / M3
        * H2(msbL2 / std::norm(M3), msbR2 / std::norm(M3))

      + sqr(gY) / (96.*sqr(M_PI)) * std::conj(mu) / M1
        * (H2(msbL2 / std::norm(M1), std::norm(mu) / std::norm(M1))
           + 2. * H2(msbR2 / std::norm(M1), std::norm(mu) / std::norm(M1)))

      + sqr(gY) / (144.*sqr(M_PI)) * std::conj(mu) / M1
        * H2(msbL2 / std::norm(M1), msbR2 / std::norm(M1))

      + 3. * sqr(g2) / (32.*sqr(M_PI)) * std::conj(mu) / M2
        * H2(msbL2 / std::norm(M2), std::norm(mu) / std::norm(M2))
      ;

   const double eps_Y_vM =
      oneOver16PiSqr * std::conj(At) / mu
      * H2(mstL2 / std::norm(mu), mstR2 / std::norm(mu))
      // term ~ self-energy neglected
      ;

   const double eps_Y = - oneOver16PiSqr * std::conj(At) / mu
      * H2(mstL2 / std::norm(mu), mstR2 / std::norm(mu))
      + eps_Y_vM;

   const double eps_3_tilde = eps_0 + sqr(yt) * eps_Y;

   const double delta_b = TB * eps_3_tilde;

   return delta_b;
}

// fermion/sfermion corrections, log-approximations

double LogNorm(const MSSMNoFV_onshell& model) {
   // function to find minimum special masses to normalize logarithms

   return fmin(std::abs(model.get_MassB()),
           fmin(std::abs(model.get_MassWB()),
            fmin(std::abs(model.get_Mu()),
             fmin(sqrt(model.get_me2(1, 1)), sqrt(model.get_ml2(1, 1))))));
}

double Deltag1(const MSSMNoFV_onshell& model) {
   const double gY = model.get_gY();
   const Eigen::Matrix<double,3,3> mu2(model.get_mu2());
   const Eigen::Matrix<double,3,3> md2(model.get_md2());
   const Eigen::Matrix<double,3,3> mq2(model.get_mq2());
   const Eigen::Matrix<double,3,3> me2(model.get_me2());
   const Eigen::Matrix<double,3,3> ml2(model.get_ml2());
   const double LogScale = LogNorm(model);

   return ( sqr(gY) * oneOver16PiSqr * 4. / 3.
            * (8. / 3. * log(sqrt(mu2(0, 0)) / LogScale)
             + 4. / 3. * log(sqrt(mu2(2, 2)) / LogScale)
             + 2. / 3. * log(sqrt(md2(0, 0)) / LogScale)
             + 1. / 3. * log(sqrt(md2(2, 2)) / LogScale)
             + 1. / 3. * log(sqrt(mq2(0, 0)) / LogScale)
             + 1. / 6. * log(sqrt(mq2(2, 2)) / LogScale)
             + log(sqrt(me2(2, 2)) / LogScale)
             + 0.5 * log(sqrt(ml2(2, 2)) / LogScale)) );
}

double DeltaYukHiggsino(const MSSMNoFV_onshell& model) {
   const double ytau = model.get_Ye(2, 2);
   const double ytop = model.get_Yu(2, 2);
   const double ybot = model.get_Yd(2, 2);
   const Eigen::Matrix<double,3,3> mu2(model.get_mu2());
   const Eigen::Matrix<double,3,3> md2(model.get_md2());
   const Eigen::Matrix<double,3,3> mq2(model.get_mq2());
   const Eigen::Matrix<double,3,3> me2(model.get_me2());
   const Eigen::Matrix<double,3,3> ml2(model.get_ml2());
   const double LogScale = LogNorm(model);

   return ( oneOver16PiSqr * 0.5
            * (3. * sqr(ytop) * log(sqrt(mu2(2, 2)) / LogScale)
             + 3. * sqr(ybot) * log(sqrt(md2(2, 2)) / LogScale)
             + 3. * (sqr(ytop) + sqr(ybot)) * log(sqrt(mq2(2, 2)) / LogScale)
             + sqr(ytau) * (log(sqrt(me2(2, 2)) / LogScale)
                           + log(sqrt(ml2(2, 2)) / LogScale))) );
}

double DeltaYukBinoHiggsino(const MSSMNoFV_onshell& model) {
   const double ytop = model.get_Yu(2, 2);
   const Eigen::Matrix<double,3,3> mu2(model.get_mu2());
   const Eigen::Matrix<double,3,3> mq2(model.get_mq2());
   const double LogScale = LogNorm(model);

   return ( oneOver16PiSqr * sqr(ytop) * ( - 8. * log(sqrt(mu2(2, 2)) / LogScale)
                            + 2. * log(sqrt(mq2(2, 2)) / LogScale)) );
}

double Deltag2(const MSSMNoFV_onshell& model) {
   const double g2 = model.get_g2();
   const Eigen::Matrix<double,3,3> mq2(model.get_mq2());
   const Eigen::Matrix<double,3,3> ml2(model.get_ml2());
   const double LogScale = LogNorm(model);

   return ( sqr(g2) * oneOver16PiSqr * 4. / 3.
            * (3. * log(sqrt(mq2(0, 0)) / LogScale)
             + 1.5 * log(sqrt(mq2(2, 2)) / LogScale)
             + 0.5 * log(sqrt(ml2(2, 2)) / LogScale)) );
}

double DeltaYukWinoHiggsino(const MSSMNoFV_onshell& model) {
   const double ytop = model.get_Yu(2, 2);
   const Eigen::Matrix<double,3,3> mq2(model.get_mq2());
   const double LogScale = LogNorm(model);

   return oneOver16PiSqr * - 6. * sqr(ytop) * log(sqrt(mq2(2, 2)) / LogScale);
}

double DeltaTanBeta(const MSSMNoFV_onshell& model) {;
   const double ytau = model.get_Ye(2, 2);
   const double ytop = model.get_Yu(2, 2);
   const double ybot = model.get_Yd(2, 2);
   const double LogScale = LogNorm(model);
   const double MUDIM = model.get_MUDIM();

   return ( oneOver16PiSqr * (sqr(ytau) - 3. * sqr(ytop) + 3. * sqr(ybot))
             * log(MUDIM / LogScale) );
}

double amuWHnu2L(const MSSMNoFV_onshell& model) {
   const double test1 = .75;

   return ( amuWHnu(model)
            * (.02 * test1 + Deltag2(model) + DeltaYukHiggsino(model)
              + DeltaYukWinoHiggsino(model) + DeltaTanBeta(model)) );
}

double amuWHmuL2L(const MSSMNoFV_onshell& model) {
   const double test2 = .75;

   return ( amuWHmuL(model)
            * (.02 * test2 + Deltag2(model) + DeltaYukHiggsino(model)
             + DeltaYukWinoHiggsino(model) + DeltaTanBeta(model))  );
}

double amuBHmuL2L(const MSSMNoFV_onshell& model) {
   const double test3 = .75;

   return ( amuBHmuL(model)
            * (.02 * test3 + Deltag1(model) + DeltaYukHiggsino(model)
              + DeltaYukBinoHiggsino(model) + DeltaTanBeta(model))  );
}

double amuBHmuR2L(const MSSMNoFV_onshell& model) {
   const double test4 = 2.;

   return ( amuBHmuR(model)
            * (.02 * test4 + Deltag1(model) + DeltaYukHiggsino(model)
              + DeltaYukBinoHiggsino(model) + DeltaTanBeta(model))  );
}

double amuBmuLmuR2L(const MSSMNoFV_onshell& model) {
   const double test5 = 1.5;

   return ( amuBmuLmuR(model)
            * (.02 * test5 + Deltag1(model) + DeltaTanBeta(model)) );
}

/**
 * No tan(beta) resummation
 */
double amu2LFSfapprox_non_tan_beta_resummed(const MSSMNoFV_onshell& model) {

   return ( amuWHnu2L(model) + amuWHmuL2L(model) + amuBHmuL2L(model)
           + amuBHmuR2L(model) + amuBmuLmuR2L(model) );
}

/**
 * Includes tan(beta) resummation
 */
double amu2LFSfapprox(const MSSMNoFV_onshell& model) {

   return ( amuWHnu2L(model) + amuWHmuL2L(model) + amuBHmuL2L(model)
	    + amuBHmuR2L(model) + amuBmuLmuR2L(model) ) * tan_beta_cor(model);
}

// photonic corrections, all

double amuChipmPhotonic(const MSSMNoFV_onshell& model) {
   double result = 0.;
   const double MM = model.get_MM();
   const Eigen::Array<double,2,1> AAC_(AAC(model));
   const Eigen::Array<double,2,1> BBC_(BBC(model));
   const Eigen::Array<double,2,1> MCha(model.get_MCha());
   const double MSvmL = model.get_MSvmL();
   const Eigen::Array<double,2,1> x__k(x_k(model));
   const double mu_DREG = model.get_MUDIM();

   for(int k=0; k<2; k++) {
      result += ( (AAC_(k) * F1C(x__k(k)) / 12.
                   + MCha(k) / 3. * BBC_(k) / model.get_MM() * F2C(x__k(k)))
                    * 16. * log(MM / MSvmL)
                  - 47. * AAC_(k) * F3C(x__k(k)) / 72.
                  - 61. * MCha(k) / 9. * BBC_(k) / model.get_MM() * F4C(x__k(k))
                  - (0.5 * AAC_(k) * F1C(x__k(k))
                   + MCha(k) * BBC_(k) / model.get_MM() * F2C(x__k(k)))
                    * log(sqr(MSvmL / mu_DREG)) );
   }

   return  sqr(model.get_EL0()) * sqr(oneOver16PiSqr) * sqr(MM / MSvmL) * result;
}

double amuChi0Photonic(const MSSMNoFV_onshell& model) {
   double result = 0.;
   const double MM = model.get_MM();
   const Eigen::Matrix<double,4,2> AAN_(AAN(model));
   const Eigen::Matrix<double,4,2> BBN_(BBN(model));
   const Eigen::Array<double,4,1> MNeu(model.get_MChi());
   const Eigen::Array<double,2,1> MSmu(model.get_MSmu());
   const Eigen::Matrix<double,4,2> x__im(x_im(model));
   const double mu_DREG = model.get_MUDIM();

   for(int i=0; i<4; ++i) {
      for(int m=0; m<2; ++m) {
         result +=  1. / sqr(MSmu(m))
                    * ((- 1. / 12. * AAN_(i, m) * F1N(x__im(i, m))
                       - MNeu(i) / 6. * BBN_(i, m) / model.get_MM() * F2N(x__im(i, m)))
                        * 16. * log(MM / MSmu(m))
                      + 35. / 72. * AAN_(i, m) * F3N(x__im(i, m))
                      + 8. * MNeu(i) / 9. * BBN_(i, m) / model.get_MM() * F4N(x__im(i, m))
                      + (0.25 * AAN_(i, m) * F1N(x__im(i, m)))
                       * log(sqr(MSmu(m) / mu_DREG)) );
      }
   }

   return sqr(model.get_EL0()) * sqr(oneOver16PiSqr) * sqr(MM) * result;
}

// amu2Loop_a corrections

/**
 * include resummation of 1/(1+Deltamu) within muon Yukawa coupling
 * but not in other Yukawa couplings (in particular not in tau,
 * bottom)
 *
 * @note the result is < 0
 */
double tan_alpha(const MSSMNoFV_onshell& model) {
   const double TB = model.get_TB();
   const double MZ = model.get_MZ();
   const double MA0 = model.get_MA0();
   const double tan2beta = 2. * TB / (1. - sqr(TB));
   const double tan2alpha = tan2beta * (sqr(MA0) + sqr(MZ)) / (sqr(MA0) - sqr(MZ));

   return - 1. / tan2alpha - sqrt(1. / sqr(tan2alpha) + 1.); // alpha < 0 !
}

/**
 * includes tan(beta) resummation
 */
Eigen::Matrix<std::complex<double>,3,3> lambda_mu_cha(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<std::complex<double>,3,3> result;
   const Eigen::Array<double,2,1> MCha(model.get_MCha());
   const double MW(model.get_MW());
   const double TB(model.get_TB());
   const double CB = 1. / sqrt(1. + sqr(TB));
   const double SB = sqrt(1. - sqr(CB));
   const double TA(tan_alpha(model));
   const double CA = 1. / sqrt(1. + sqr(TA));
   const double SA = - sqrt(1. - sqr(CA));
   const Eigen::Matrix<std::complex<double>,2,2> U(model.get_UM());
   const Eigen::Matrix<std::complex<double>,2,2> V(model.get_UP());
   const double one_over_cb_eff = sqrt(2.) * model.get_Ye(1,1)
      * model.get_MW() / model.get_MM() / model.get_g2();

   for(int k=0; k<2; ++k) {
      result(k, 0) = ( std::sqrt(2.) * MW / MCha(k)
                      * (U(k, 0) * V(k, 1) * CA + U(k, 1) * V(k, 0) * (-SA)) );
      result(k, 1) = ( std::sqrt(2.) * MW / MCha(k)
                      * (U(k, 0) * V(k, 1) * SA + U(k, 1) * V(k, 0) * CA) );
      result(k, 2) = ( std::sqrt(2.) * MW / MCha(k)
                      * (U(k, 0) * V(k, 1) * (-CB) + U(k, 1) * V(k, 0) * (-SB)) );
   }
   result(2, 0) = -SA * one_over_cb_eff;
   result(2, 1) = CA * one_over_cb_eff;
   result(2, 2) = SB * one_over_cb_eff;

   return result;
}

Eigen::Matrix<std::complex<double>,2,2> lambda_stop(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<std::complex<double>,2,2> result;
   const double TB(model.get_TB());
   const double CB = 1. / sqrt(1. + sqr(TB));
   const double SB = sqrt(1. - sqr(CB));
   const double TA(tan_alpha(model));
   const double CA = 1. / sqrt(1. + sqr(TA));
   const double SA = - sqrt(1. - sqr(CA));
   const double MT(model.get_MT());
   const Eigen::Array<double,2,1> MStop(model.get_MStop());
   const Eigen::Matrix<double,2,2> UStop(model.get_UStop());
   const double At(model.get_Au(2, 2));
   const double Mu(model.get_Mu());

   for(int i=0; i<2; ++i) {
      result(i, 0) = 2. * MT / (sqr(MStop(i)) * SB) * (Mu * SA + At * CA)
                      * std::conj(UStop(i, 0)) * UStop(i, 1);
      result(i, 1) = 2. * MT / (sqr(MStop(i)) * SB) * (Mu * (-CA) + At * SA)
                      * std::conj(UStop(i, 0)) * UStop(i, 1);
   }

   return result;
}

/**
 * includes tan(beta) resummation
 */
Eigen::Matrix<std::complex<double>,2,2> lambda_sbot(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<std::complex<double>,2,2> result;
   const double TA(tan_alpha(model));
   const double CA = 1. / sqrt(1. + sqr(TA));
   const double SA = - sqrt(1. - sqr(CA));
   const Eigen::Array<double,2,1> MSbot(model.get_MSbot());
   const Eigen::Matrix<double,2,2> USbot(model.get_USbot());
   const double Ab(model.get_Ad(2, 2));
   const double Mu(model.get_Mu());
   const double mb_over_cb_eff = sqrt(2.) * model.get_Yd(2,2)
      * model.get_MW() / model.get_g2();

   for(int i=0; i<2; ++i) {
      result(i, 0) = 2. * mb_over_cb_eff / (sqr(MSbot(i))) * (- Mu * CA + Ab * (-SA))
                      * std::conj(USbot(i, 0)) * USbot(i, 1);
      result(i, 1) = 2. * mb_over_cb_eff / (sqr(MSbot(i))) * (- Mu * SA + Ab * CA)
                      * std::conj(USbot(i, 0)) * USbot(i, 1);
   }

   return result;
}

/**
 * includes tan(beta) resummation
 */
Eigen::Matrix<std::complex<double>,2,2> lambda_stau(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<std::complex<double>,2,2> result;
   const double TA(tan_alpha(model));
   const double CA = 1. / sqrt(1. + sqr(TA));
   const double SA = - sqrt(1. - sqr(CA));
   const Eigen::Array<double,2,1> MStau(model.get_MStau());
   const Eigen::Matrix<double,2,2> UStau(model.get_UStau());
   const double Al(model.get_Ae(2, 2));
   const double Mu(model.get_Mu());
   const double mtau_over_cb_eff = sqrt(2.) * model.get_Ye()(2,2)
      * model.get_MW() / model.get_g2();

   for(int i=0; i<2; ++i) {
      result(i, 0) = 2. * mtau_over_cb_eff / (sqr(MStau(i))) * (- Mu * CA + Al * (-SA))
                      * std::conj(UStau(i, 0)) * UStau(i, 1);
      result(i, 1) = 2. * mtau_over_cb_eff / (sqr(MStau(i))) * (- Mu * SA + Al * CA)
                      * std::conj(UStau(i, 0)) * UStau(i, 1);
   }

   return result;
}

double amua2LSferm(const MSSMNoFV_onshell& model) {
   double result = 0.;
   const double MM(model.get_MM());
   const double MW(model.get_MW());
   const double SW = sqrt(1. - sqr(MW / model.get_MZ()));
   const double EL(model.get_EL());
   const Eigen::Array<double,2,1> m_stop(model.get_MStop());
   const Eigen::Array<double,2,1> m_sbot(model.get_MSbot());
   const Eigen::Array<double,2,1> m_stau(model.get_MStau());
   const Eigen::Array<double,2,1> M_higgs(model.get_Mhh());
   Eigen::Array<std::complex<double>,2,1> lambda_mu;
   lambda_mu(0) = lambda_mu_cha(model)(2, 0);
   lambda_mu(1) = lambda_mu_cha(model)(2, 1);
   const Eigen::Matrix<std::complex<double>,2,2> lambdastop(lambda_stop(model));
   const Eigen::Matrix<std::complex<double>,2,2> lambdasbot(lambda_sbot(model));
   const Eigen::Matrix<std::complex<double>,2,2> lambdastau(lambda_stau(model));

   double N_c = 3.;
   double Q = 2. / 3.;

   for(int i=0; i<2; ++i) {
      for(int s=0; s<2; ++s) {
         result += ( N_c * sqr(Q) * real(lambda_mu(s)
                   * lambdastop(i, s))
                   * f_sferm(sqr(m_stop(i) / M_higgs(s))) );\
      }
   }

   Q = - 1. / 3.;
   for(int i=0; i<2; ++i) {
      for(int s=0; s<2; ++s) {
         result += ( N_c * sqr(Q) * real(lambda_mu(s)
                   * lambdasbot(i, s))
                   * f_sferm(sqr(m_sbot(i) / M_higgs(s))) );
      }
   }

   N_c = 1.;
   Q = - 1.;
   for(int i=0; i<2; ++i) {
      for(int s=0; s<2; ++s) {
         result += ( N_c * sqr(Q) * real(lambda_mu(s)
                   * lambdastau(i, s))
                   * f_sferm(sqr(m_stau(i) / M_higgs(s))) );
      }
   }

   return result * 2. * sqr(oneOver16PiSqr * sqr(EL) * MM / (MW * SW));
}

double amua2LCha(const MSSMNoFV_onshell& model) {
   double result = 0.;
   const double MM(model.get_MM());
   const double MW(model.get_MW());
   const double MA0(model.get_MA0());
   const Eigen::Array<double,2,1> M_higgs(model.get_Mhh());
   const double SW = sqrt(1. - sqr(MW / model.get_MZ()));
   const double EL(model.get_EL());
   const Eigen::Array<double,2,1> m_cha(model.get_MCha());
   const Eigen::Matrix<std::complex<double>,3,3> lambda(lambda_mu_cha(model));

   for(int k=0; k<2; ++k) {
      result += real(lambda(2, 2) * lambda(k, 2)) * f_PS(sqr(m_cha(k) / MA0));
      for(int s=0; s<2; ++s) {
         result += ( real(lambda(2, s) * lambda(k, s))
                   * f_S(sqr(m_cha(k) / M_higgs(s))) );
      }
   }

   return result * 2. * sqr(oneOver16PiSqr * sqr(EL) * MM / (MW * SW));
}

} // namespace gm2os
} // namespace flexiblesusy
