
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

#endif
