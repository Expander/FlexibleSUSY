
#ifndef TEST_NUTNMSSM_H
#define TEST_NUTNMSSM_H

#include "wrappers.hpp"
#include "nmssmsoftsusy.h"
#include "NUTNMSSM_two_scale_model.hpp"

using namespace flexiblesusy;
using namespace softsusy;

void copy_parameters(const NUTNMSSM<Two_scale>& nmssm, NmssmSoftsusy& softsusy)
{
   // copy base class parameters
   softsusy.setLoops(nmssm.get_loops());
   softsusy.setMu(nmssm.get_scale());
   softsusy.setThresholds(nmssm.get_thresholds());

   // copy susy parameters
   softsusy.setGaugeCoupling(1, nmssm.get_g1());
   softsusy.setGaugeCoupling(2, nmssm.get_g2());
   softsusy.setGaugeCoupling(3, nmssm.get_g3());

   const double vu = nmssm.get_vu();
   const double vd = nmssm.get_vd();
   const double tanBeta = vu / vd;
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));
   const double vs = nmssm.get_vS();

   softsusy.setTanb(tanBeta);
   softsusy.setHvev(vev);
   softsusy.setSvev(vs);

   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         softsusy.setYukawaElement(YU, i, k, nmssm.get_Yu()(i-1,k-1));
         softsusy.setYukawaElement(YD, i, k, nmssm.get_Yd()(i-1,k-1));
         softsusy.setYukawaElement(YE, i, k, nmssm.get_Ye()(i-1,k-1));
      }
   }

   softsusy.setLambda(nmssm.get_Lambdax());
   softsusy.setKappa(nmssm.get_Kappa());

   // copy soft parameters
   softsusy.setGauginoMass(1, nmssm.get_MassB());
   softsusy.setGauginoMass(2, nmssm.get_MassWB());
   softsusy.setGauginoMass(3, nmssm.get_MassG());

   softsusy.setMh1Squared(nmssm.get_mHd2());
   softsusy.setMh2Squared(nmssm.get_mHu2());
   softsusy.setMsSquared(nmssm.get_ms2());

   for (int i = 1; i <= 3; i++) {
      for (int k = 1; k <= 3; k++) {
         softsusy.setSoftMassElement(mQl, i, k,  nmssm.get_mq2()(i-1,k-1));
         softsusy.setSoftMassElement(mUr, i, k,  nmssm.get_mu2()(i-1,k-1));
         softsusy.setSoftMassElement(mDr, i, k,  nmssm.get_md2()(i-1,k-1));
         softsusy.setSoftMassElement(mLl, i, k,  nmssm.get_ml2()(i-1,k-1));
         softsusy.setSoftMassElement(mEr, i, k,  nmssm.get_me2()(i-1,k-1));

         softsusy.setTrilinearElement(UA, i, k, nmssm.get_TYu()(i-1,k-1));
         softsusy.setTrilinearElement(DA, i, k, nmssm.get_TYd()(i-1,k-1));
         softsusy.setTrilinearElement(EA, i, k, nmssm.get_TYe()(i-1,k-1));
      }
   }

   softsusy.setTrialambda(nmssm.get_TLambdax());
   softsusy.setTriakappa(nmssm.get_TKappa());

   softsusy.setSusyMu(0.);
   softsusy.setM3Squared(0.);
   softsusy.setMspSquared(0.);
   softsusy.setXiS(0.);
   softsusy.setXiF(0.);
   softsusy.setMupr(0.);

   softsusy.setMw(softsusy.displayMwRun());

   softsusy::Z3 = true;
}

void setup_NUTNMSSM_const(NUTNMSSM<Two_scale>& m, NmssmSoftsusy& s, const NUTNMSSM_input_parameters& input)
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

#endif
