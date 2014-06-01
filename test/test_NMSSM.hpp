
#ifndef TEST_NMSSM_H
#define TEST_NMSSM_H

#include "wrappers.hpp"
#include "ew_input.hpp"
#include "nmssmsoftsusy.h"
#include "NMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void setup_NMSSM_const(NMSSM<Two_scale>& m, NmssmSoftsusy& s, const NMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double lambda = input.LambdaInput;
   const double kappa = 0.01;
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double vS = 1000.;
   const double susyMu = 0;
   const double BMu = 0;
   const double scale = Electroweak_constants::MZ;

   DoubleMatrix Yu_SS(3,3), Yd_SS(3,3), Ye_SS(3,3);
   Yu_SS(3,3) = 165.0   * root2 / (vev * sinBeta);
   Yd_SS(3,3) = 2.9     * root2 / (vev * cosBeta);
   Ye_SS(3,3) = 1.77699 * root2 / (vev * cosBeta);
   DoubleMatrix ID(3, 3), mm0_SS(3, 3);
   for (int i=1; i<=3; i++) ID(i, i) = 1.0;
   mm0_SS = ID * sqr(m0);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_Lambdax(lambda);
   m.set_Kappa(kappa);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_ms2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_TLambdax(a0 * lambda);
   m.set_TKappa(a0 * kappa);
   // m.set_Mu(susyMu);
   // m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);
   m.set_vS(vS);

   s.setMu(scale);
   s.setLoops(1);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setYukawaMatrix(YU, Yu_SS);
   s.setYukawaMatrix(YD, Yd_SS);
   s.setYukawaMatrix(YE, Ye_SS);
   s.setLambda(lambda);
   s.setKappa(kappa);
   s.setGauginoMass(1, M12);
   s.setGauginoMass(2, M12);
   s.setGauginoMass(3, M12);
   s.setSoftMassMatrix(mQl, mm0_SS);
   s.setSoftMassMatrix(mUr, mm0_SS);
   s.setSoftMassMatrix(mDr, mm0_SS);
   s.setSoftMassMatrix(mLl, mm0_SS);
   s.setSoftMassMatrix(mEr, mm0_SS);
   s.setMh1Squared(sqr(m0));
   s.setMh2Squared(sqr(m0));
   s.setMsSquared(sqr(m0));
   s.setTrilinearMatrix(UA, a0 * Yu_SS);
   s.setTrilinearMatrix(DA, a0 * Yd_SS);
   s.setTrilinearMatrix(EA, a0 * Ye_SS);
   s.setTrialambda(a0 * lambda);
   s.setTriakappa(a0 * kappa);
   s.setSusyMu(susyMu);
   s.setM3Squared(BMu);
   s.setHvev(vev);
   s.setTanb(tanBeta);
   s.setSvev(vS);
   s.setMw(s.displayMwRun());
}

void setup_NMSSM(NMSSM<Two_scale>& m, NmssmSoftsusy& s, NMSSM_input_parameters& input)
{
   input.m0 = 200.;
   input.m12 = 200.;
   input.TanBeta = 10.;
   input.Azero = -500.;
   input.LambdaInput = 0.1;
   input.SignvS = 1;

   setup_NMSSM_const(m, s, input);
}

void ensure_n_loop_ewsb(NMSSM<Two_scale>& m, int loop_level)
{
   const double precision = m.get_ewsb_iteration_precision();
   m.set_ewsb_loop_order(loop_level);
   m.solve_ewsb();

   if (loop_level == 1) {
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1() - m.tadpole_hh(0).real(), precision);
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2() - m.tadpole_hh(1).real(), precision);
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3() - m.tadpole_hh(2).real(), precision);
   } else if (loop_level == 2) {
      double two_loop_tadpole[3];
      m.tadpole_hh_2loop(two_loop_tadpole);
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1() - Re(m.tadpole_hh(0)) - two_loop_tadpole[0], precision);
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_2() - Re(m.tadpole_hh(1)) - two_loop_tadpole[1], precision);
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_3() - Re(m.tadpole_hh(2)) - two_loop_tadpole[2], precision);
   }
}

void ensure_n_loop_ewsb(NmssmSoftsusy& s, int loop_level)
{
   const int signMu = s.displaySusyMu() >= 0.0 ? 1 : -1;
   const double mtrun = s.displayDrBarPars().mt;
   softsusy::numRewsbLoops = loop_level;
   s.rewsb(signMu, mtrun);
}

#endif
