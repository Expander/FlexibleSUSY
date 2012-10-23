
/** \file SMsusy.cpp
   - Project:     cE6SSM specgen
   - Author:      Peter Athron
*/

#include "sm.hpp"
#include "lowe.h"

//Peter:: New ESSM versions of constructors
StandardModel::StandardModel()
   : SMu(3, 3), SMd(3, 3), SMe(3, 3), g(3), smu(0.0), tanb(0.0), hVev(0.0), lambda(3), kappa(3), MN_SUSY(3), h_N(0.0), mu_0(0.0), gdash_1(0.0), g_11(0.0)
{
   setPars(numSusySMPars);
   setMu(0.0);
   setLoops(2);
   setThresholds(0);
}

StandardModel::StandardModel(const StandardModel &s)
   : SMu(s.SMu), SMd(s.SMd), SMe(s.SMe), g(s.g), smu(s.smu), tanb(s.tanb), hVev(s.hVev),  lambda(s.lambda), kappa(s.kappa), MN_SUSY(s.MN_SUSY), h_N(s.h_N), mu_0(s.mu_0), gdash_1(s.gdash_1), g_11(s.g_11)
{
   setPars(numSusySMPars);
   setMu(s.displayMu());
   setLoops(s.displayLoops());
   setThresholds(s.displayThresholds());
}

StandardModel::StandardModel(const DoubleMatrix & SMu, const DoubleMatrix & SMd, const
               DoubleMatrix & SMe, const DoubleVector & v, double m,
               double tb, double MU, int l, int t, double hv, DoubleVector lam, DoubleVector kap, DoubleVector MN_SUSY, double h_N, double mu_0, double gd, double g_11)
   : SMu(SMu), SMd(SMd), SMe(SMe), g(v), smu(m), tanb(tb), hVev(hv), lambda(lam), kappa(kap), MN_SUSY(MN_SUSY), h_N(h_N), mu_0(mu_0), gdash_1(gd), g_11(g_11)
{
   setPars(numSusySMPars);
   setMu(MU);
   setLoops(l);
   setThresholds(t);
}

const StandardModel & StandardModel::operator=(const StandardModel & s)
{
   if (this == &s) return *this;
   SMu = s.SMu;
   SMd = s.SMd;
   SMe = s.SMe;
   smu = s.smu;
   tanb = s.tanb;
   g = s.g;
   setMu(s.displayMu());
   setLoops(s.displayLoops());
   setThresholds(s.displayThresholds());
   hVev = s.hVev;
   kappa = s.kappa;
   lambda = s.lambda;
   h_N = s.h_N;
   mu_0 = s.mu_0;
   gdash_1 = s.gdash_1;
   g_11 = s.g_11;
   MN_SUSY = s.MN_SUSY;
   return *this;
}

void StandardModel::setSomePars(const StandardModel & s)
{
   SMu = s.SMu;
   SMd = s.SMd;
   SMe = s.SMe;
   g = s.g;
}

void StandardModel::setSMYukawaElement(yukawa k, int i, int j, double f)
{
   switch (k) {
   case YU:
      SMu(i, j) = f;
      break;
   case YD:
      SMd(i, j) = f;
      break;
   case YE:
      SMe(i, j) = f;
      break;
   default:
      ostringstream ii;
      ii << "StandardModel::set called with illegal " << int(k) << "\n";
      throw ii.str();
      break;
   }
}

void StandardModel::setSMYukawaMatrix(yukawa k, const DoubleMatrix & m)
{
   switch (k) {
   case YU:
      SMu = m;
      break;
   case YD:
      SMd = m;
      break;
   case YE:
      SMe = m;
      break;
   default:
      ostringstream ii;
      ii << "StandardModel::set called with illegal " << int(k) << "\n";
      throw ii.str();
      break;
   }
}

double StandardModel::displayYukawaElement(yukawa k, int i, int j) const
{
   switch (k) {
   case YU:
      return SMu.display(i, j);
      break;
   case YD:
      return SMd.display(i, j);
      break;
   case YE:
      return SMe.display(i, j);
      break;
   default:
      ostringstream ii;
      ii << "StandardModel::display called with illegal " << int(k) << "\n";
      throw ii.str();
      break;
   }
   return 0.0;
}

DoubleMatrix StandardModel::displayYukawaMatrix(yukawa k) const
{
   switch (k) {
   case YU:
      return SMu;
      break;
   case YD:
      return SMd;
      break;
   case YE:
      return SMe;
      break;
   default:
      ostringstream ii;
      ii << "StandardModel::display called with illegal " << int(k) << "\n";
      throw ii.str();
      break;
   }
}

//Peter:: edited to include new Essm parameters

const DoubleVector StandardModel::display() const
{
   DoubleVector y(numSusySMPars);
   int i, j, k = 0;
   for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++) {
         k++;
         y(k) = SMu.display(i, j);
         y(k + 9) = SMd.display(i, j);
         y(k + 18) = SMe.display(i, j);
      }
   k = 27;
   for (i = 1; i <= 3; i++) {
      k++;
      y(k) = g.display(i);
   }
   y(31) = smu;
   y(32) = tanb;
   y(33) = hVev;

   //Peter:: New Essm parameters
   k = 33;
   for (i = 1; i <= 3; i++) {
      k++;
      y(k) = lambda.display(i);
   }

   for (i = 1; i <= 3; i++) {
      k++;
      y(k) = kappa.display(i);
   }

   y(40) = h_N;
   y(41) = mu_0;
   y(42) = gdash_1;
   y(43) = g_11;
   y(44) = MN_SUSY.display(1);
   y(45) = MN_SUSY.display(2);
   y(46) = MN_SUSY.display(3);
   return y;
}
//Peter:: edited to include new Essm parameters

void StandardModel::set(const DoubleVector & y)
{
   int i, j, k = 0;
   for (i = 1; i <= 3; i++)
      for (j = 1; j <= 3; j++) {
         k++;
         SMu(i, j) = y.display(k);
         SMd(i, j) = y.display(k + 9);
         SMe(i, j) = y.display(k + 18);
      }
   k = 27;
   for (i = 1; i <= 3; i++) {
      k++;
      g(i) = y.display(k);
   }
   smu = y.display(31);
   tanb = y.display(32);
   hVev = y.display(33);
   k = 33;
   for (i = 1; i <= 3; i++) {
      k++;
      lambda(i) = y.display(k);
   }
   for (i = 1; i <= 3; i++) {
      k++;
      kappa(i) = y.display(k);
   }

   h_N = y.display(40);
   mu_0 = y.display(41);
   gdash_1 = y.display(42);
   g_11 = y.display(43);
   k = 43;
   for (i = 1; i <= 3; i++) {
      k++;
      MN_SUSY(i) = y.display(k);
   }

}






double StandardModel::displayTanb() const
{
   return tanb;
}

//Peter:: new Sm display functions

double StandardModel::displayh_N() const
{
   return h_N;
}
double StandardModel::displaymu_0() const
{
   return mu_0;
}
double StandardModel::displaygdash_1() const
{
   return gdash_1;
}
double StandardModel::displayg_11() const
{
   return g_11;
}
double StandardModel::displayMN_SUSY(int i) const
{
   return MN_SUSY.display(i);
}
DoubleVector StandardModel::displaykappa() const
{
   return kappa;
}
DoubleVector StandardModel::displaylambda() const
{
   return lambda;
}

ostream & operator <<(ostream &left, const StandardModel &s)
{
   left << "Supersymmetric parameters at Q: " << s.displayMu() << endl;
   left << " Y^U" << s.displayYukawaMatrix(YU) << " Y^D" <<
        s.displayYukawaMatrix(YD) << " Y^E" << s.displayYukawaMatrix(YE);
   left << "higgs VEV: " << s.displayHvev()
        << " tan beta: " << s.displayTanb() << " smu: " << s.displaySusyMu() <<
        "\n";
   left << "g1: " << s.displayGaugeCoupling(1) << " g2: " <<
        s.displayGaugeCoupling(2) << " g3: " <<
        s.displayGaugeCoupling(3) << endl;
   left << "thresholds: " << s.displayThresholds()
        << " #loops: " << s.displayLoops() << '\n';
   return left;
}

void StandardModel::setSusy(const StandardModel & s)
{
   setLoops(s.displayLoops());
   setThresholds(s.displayThresholds());
   setMu(s.displayMu());
   setSMYukawaMatrix(YU, s.displayYukawaMatrix(YU));
   setSMYukawaMatrix(YD, s.displayYukawaMatrix(YD));
   setSMYukawaMatrix(YE, s.displayYukawaMatrix(YE));
   setHvev(s.displayHvev());
   setTanb(s.displayTanb());
   setSusyMu(s.displaySusyMu());
   setAllGauge(s.displayGauge());
   seth_N(s.displayh_N());
   setgash_1(s.displaygdash_1());
   setg_11(s.displayg_11());
   setmu_0(s.displaymu_0());
   setkappa(s.displaykappa());
   setlambda(s.displaylambda());
   DoubleVector MN(3);

   int gen;
   for (gen = 1; gen < 4; gen++) MN(gen) = s.displayMN_SUSY(gen);
   setMN_SUSY(MN);

}

istream & operator >>(istream &left, StandardModel &s)
{
   char c[70];
   DoubleMatrix u(3, 3), d(3, 3), e(3, 3);
   double g1, g2, g3, smu, mu, tanb, hv;
   int loops, thresh;
   left >> c >> c >> c >> c >> mu;
   left >> c >> u >> c >> d >> c >> e >> c >> c >> hv;
   left >> c >> c >> tanb >> c >> smu;
   left >> c >> g1 >> c >> g2 >> c >> g3;
   left >> c >> thresh >> c >> loops;
   s.setSMYukawaMatrix(YU, u);
   s.setSMYukawaMatrix(YD, d);
   s.setSMYukawaMatrix(YE, e);
   s.setHvev(hv);
   s.setTanb(tanb);
   s.setGaugeCoupling(1, g1);
   s.setGaugeCoupling(2, g2);
   s.setGaugeCoupling(3, g3);
   s.setThresholds(thresh);
   s.setSusyMu(smu);
   s.setMu(mu);
   s.setLoops(loops);
   return left;
}



// Outputs derivatives (DRbar scheme) in the form of ds. a contains the
// matrices calculated that are handy for computation.
// W=  LL Y^E H1 ER + QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1
// is the superpotential. Consistent with Allanach, Dedes, Dreiner
// hep-ph/9902251 and Barger, Berger and Ohmann hep-ph/9209232, 9311269
// EXCEPT for the sign of smu, which is opposite. These equations are also
// valid for W=  - LL Y^E H1 ER - QL Y^D H1 DR + QL Y^U H2 UR + smu H2 H1, the
// new SOFTSUSY convention
StandardModel StandardModel::beta(sBrevity & a) const
{
   // Wave function renormalisations: convention for g**(i, j) is that i is the
   // LOWER index and j the upper in our paper hep-ph/9902251
   // Peter: In third family approximation we are only using these as
   // numbers not matrices
   static DoubleMatrix gEE(3, 3), gLL(3, 3), gQQ(3, 3), gDD(3, 3), gUU(3, 3);
   double gH1H1, gH2H2;
   static DoubleVector dg(1, 3);
   static double dgdash_1, dg_11;
   // keep this option in order to interface with RPVSUSY
   // Peter:: keep this to do the gauge couplings but then ignore the
   // matrices gEE etc and replace with 3rd gen numbers
   anomalousDimension(gEE, gLL, gQQ, gUU, gDD, dg, dgdash_1, dg_11, gH1H1, gH2H2, a);

   // To keep this a const function
   const DoubleMatrix &SMu1 = SMu.display(), &SMd1 = SMd.display(), &SMe1 = SMe.display();

   // contain derivatives of up, down quarks and leptons
   static DoubleMatrix dSMu(3, 3), dSMd(3, 3), dSMe(3, 3);
   // Peter:: initialise all entries to zero since we are going to only
   // use the third generation i don't want any funny preset values for
   // the other betas screwing up the evolution
   int just_mak, ing_sure;


   for (just_mak = 1; just_mak < 4; just_mak++) {
      for (ing_sure = 1; ing_sure < 4; ing_sure++) {
         dSMu(just_mak, ing_sure) = 0;
         dSMd(just_mak, ing_sure) = 0;
         dSMe(just_mak, ing_sure) = 0;
      }
   }

// mu parameter derivatives
   double dmu;
   // DoubleMatrix mu_tilde(3,3), dmu_tilde(3,3);

   // RGEs of SUSY parameters//Peter:: for MSSM
   /* du = u1 * (gUU + gH2H2) + gQQ * u1;
      dd = d1 * (gDD + gH1H1) + gQQ * d1;
      de = e1 * (gEE + gH1H1) + gLL * e1;*/
   //Peter:: Yukawas for ESSM
   //Peter:: note a 1/16pi^2 factor is contained in gUU etc
   DoubleVector lambda = displaylambda();
   DoubleVector kappa = displaykappa();
   //double lambda = lambdan(3) ; double kappa = kappan(3);
   double h_N = displayh_N();
   DoubleVector MN_SUSY(3) ;
   MN_SUSY(1)  = displayMN_SUSY(1);
   MN_SUSY(2)  = displayMN_SUSY(2);
   MN_SUSY(3)  = displayMN_SUSY(3);
   bool eta_N = 0;
   // Peter:: to implement this threshold codition for the nuetrino
   // coupling in the rges we need to be able to access the the current
   // scale.  To a firts approximation we may simlt ignore these terms.
   if (displayMu() > displayMN_SUSY(3)) eta_N = 1;
   else eta_N = 0;
   double dmu_0, dh_N;
   DoubleVector dlambda(3), dkappa(3), dMN_SUSY(3);

   // Peter:: since we are using the third family approximation the
   // matrices gEE ect need to be changed appropriately
   double  oneO16Pisq = 1.0 / (16.0 * PI * PI);

   dSMu(3, 3) = oneO16Pisq * SMu1.display(3, 3) * (-17.0 / (20.0) * sqr(displayGaugeCoupling(1)) - 9.0 / (4.0) * sqr(displayGaugeCoupling(2)) - 8.0 * sqr(displayGaugeCoupling(3)) + 4.5 * sqr(SMu1.display(3, 3)) + 1.5 * sqr(SMd1.display(3, 3)) + sqr(SMe1.display(3, 3))) ;


   dSMd(3, 3) = oneO16Pisq * SMd1.display(3, 3) * (-0.25 * sqr(displayGaugeCoupling(1)) - 9.0 / (4.0) * sqr(displayGaugeCoupling(2)) - 8.0 * sqr(displayGaugeCoupling(3)) + 1.5 * sqr(SMu1.display(3, 3)) + 4.5 * sqr(SMd1.display(3, 3)) + sqr(SMe1.display(3, 3))) ;

   dSMe(3, 3) = oneO16Pisq * SMe1.display(3, 3) * (-9.0 / (4.0) * sqr(displayGaugeCoupling(1)) - 9.0 / (4.0) * sqr(displayGaugeCoupling(2)) + 3.0 * sqr(SMu1.display(3, 3)) + 3.0 * sqr(SMd1.display(3, 3)) + 2.5 * sqr(SMe1.display(3, 3))) ;



   dmu = smu * (gH1H1 + gH2H2);

   //cout << "de = " << de << endl;

   dh_N =   1.0 / (16.0 * sqr(PI)) * h_N * (lambda(3) * lambda(3) + 3 * SMu1.display(3, 3) * SMu1.display(3, 3) + sqr(SMe1.display(3, 3)) + 2 * h_N * h_N * (1 + eta_N) - 3 * displayGaugeCoupling(2) * displayGaugeCoupling(2) - 0.6 * displayGaugeCoupling(1) * displayGaugeCoupling(1) - 0.4 * displaygdash_1() * displaygdash_1());

   dmu_0 = 1.0 / (16.0 * sqr(PI)) * mu_0 * (-3 * displayGaugeCoupling(2) * displayGaugeCoupling(2) - 0.6 * displayGaugeCoupling(1) * displayGaugeCoupling(1) - 0.4 * displaygdash_1() * displaygdash_1());
   int gen;

// Peter:: origionally wrote this out from rges in ESSM theory and
// Phenomonology paper. There the charges, Q_tilde etc apear
// explicitly.  In the rge copy Roman has put in numerical value for
// them so i have updated the beta correspondingly

   for (gen = 1; gen < 3; gen++) {
      dlambda(gen) = 1.0 / (16.0 * sqr(PI)) * lambda(gen) * (2 * lambda(gen) * lambda(gen) + 2 * (lambda(1) * lambda(1) + lambda(2) * lambda(2) + lambda(3) * lambda(3)) + 3 * (kappa(1) * kappa(1) + kappa(2) * kappa(2) + kappa(3) * kappa(3))  - 3 * displayGaugeCoupling(2) * displayGaugeCoupling(2) - 0.6 * displayGaugeCoupling(1) * displayGaugeCoupling(1) - 1.9 * displaygdash_1() * displaygdash_1()) ;
   }

   dlambda(3) = 1.0 / (16.0 * sqr(PI)) * lambda(3) * (2 * lambda(3) * lambda(3) + 2 * (lambda(1) * lambda(1) + lambda(2) * lambda(2) + lambda(3) * lambda(3)) + 3 * (kappa(1) * kappa(1) + kappa(2) * kappa(2) + kappa(3) * kappa(3))  - 3 * displayGaugeCoupling(2) * displayGaugeCoupling(2) - 0.6 * displayGaugeCoupling(1) * displayGaugeCoupling(1) - 1.9 * displaygdash_1() * displaygdash_1() + 3 * SMu1.display(3, 3) * SMu1.display(3, 3) + 3 * SMd1.display(3, 3) * SMd1.display(3, 3) + SMe1.display(3, 3) * SMe1.display(3, 3) + h_N * h_N * eta_N);


   for (gen = 1; gen < 4; gen++) {
      dkappa(gen) = 1.0 / (16.0 * sqr(PI)) * kappa(gen) * (2 * kappa(gen) * kappa(gen) + 2 * (lambda(1) * lambda(1) + lambda(2) * lambda(2) + lambda(3) * lambda(3)) + 3 * (kappa(1) * kappa(1) + kappa(2) * kappa(2) + kappa(3) * kappa(3)) - (float)16 / (float)3 * displayGaugeCoupling(3) * displayGaugeCoupling(3) - (float)4 / (float)15 * displayGaugeCoupling(1) * displayGaugeCoupling(1) - 1.9 * displaygdash_1() * displaygdash_1());
   }


   //Peter:: below is unedited thus far.

   dMN_SUSY(1) = 0;
   dMN_SUSY(2) = 0;
   dMN_SUSY(3) =  oneO16Pisq * 4 * h_N * h_N * MN_SUSY(3);

   //cout << "dMN_SUSY(3) = " <<  dMN_SUSY(3) << endl;

   double cosb2 = sqr(cos(atan(tanb))), sinb2 = 1.0 - cosb2;

   // Following is from hep-ph/9308335
   double dt = displayTanb() * (gH1H1 - gH2H2);

   double dHvev = -hVev * (cosb2 * gH1H1 + sinb2 * gH2H2);

   StandardModel ds(dSMu, dSMd, dSMe, dg, dmu, dt, displayMu(), displayLoops(),
             displayThresholds(), dHvev, dlambda, dkappa, dMN_SUSY, dh_N, dmu_0, dgdash_1, dg_11);

   // cout << "ds = " << ds.display() << endl;
   return ds;
}

void setSMBetas(DoubleMatrix & babBeta, DoubleVector &cuBeta, DoubleVector
                & cdBeta, DoubleVector & ceBeta, DoubleVector & bBeta)
{
   // 1 loop gauge beta fns //Peter:: for MSSM
   //  bBeta(1) = 33.0 / 5.0; bBeta(2) = 1.0; bBeta(3) = -3.0;

   //Peter:: 1 loop gauge beta fns for SM
   bBeta(1) = 41.0 / 10.0;
   bBeta(2) = -19.0 / 6.0;
   bBeta(3) = -7.0;

   //ESSM ones
   babBeta(1, 1) = 199.0 / 25.0;
   babBeta(1, 2) = 27.0 / 5.0;
   babBeta(1, 3) = 88.0 / 5.0;
   babBeta(2, 1) = 9.0 / 5.0;
   babBeta(2, 2) = 25.0;
   babBeta(2, 3) = 24.0;
   babBeta(3, 1) = 11.0 / 5.0;
   babBeta(3, 2) = 9.0;
   babBeta(3, 3) = 14.0;
   cuBeta(1) = 26.0 / 5.0;
   cuBeta(2) = 6.0;
   cuBeta(3) = 4.0;
   cdBeta(1) = 14.0 / 5.0;
   cdBeta(2) = 6.0;
   cdBeta(3) = 4.0;
   ceBeta(1) = 18.0 / 5.0;
   ceBeta(2) = 2.0;
   ceBeta(3) = 0.0;

}

// outputs one-loop anomlous dimensions gii given matrix inputs
// Note that we use the convention (for matrices in terms of gamma's)
// gamma^Li_Lj = M_ij for LH fields and
// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
// conjugates of the RH fields): CHECKED 23/5/02

void StandardModel::getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
                          DoubleMatrix & gQQ, DoubleMatrix & gDD,
                          DoubleMatrix & gUU, double & gH1H1, double &
                          gH2H2, sBrevity & a) const
{
   // For calculational brevity
   DoubleMatrix &u2 = a.u2, &d2 = a.d2, &e2 = a.e2, &u2t = a.u2t, &d2t = a.d2t,
                 &e2t = a.e2t;
   double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
   DoubleVector &gsq = a.gsq;

   static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));



   gEE = oneO16Pisq * (2.0 * e2t - 1.2 * gsq(1));
   gLL = oneO16Pisq * (e2 - (0.3 * gsq(1) + 1.5 * gsq(2)));
   gQQ = oneO16Pisq * (d2 + u2 - (gsq(1) / 30.0 + 1.5 * gsq(2) + 8 *
                                  gsq(3) / 3.0));
   gUU = oneO16Pisq * (2.0 * u2t - (8 * gsq(1) / 15.0 + 8 * gsq(3) /
                                    3.0));
   gDD = oneO16Pisq * (2.0 * d2t -
                       (2 * gsq(1) / 15.0 + 8 * gsq(3) / 3.0));
   gH1H1 = oneO16Pisq * (3.0 * ddT + eeT - (0.3 * gsq(1) + 1.5 *
                                            gsq(2)));
   gH2H2 = oneO16Pisq * (3.0 * uuT - (0.3 * gsq(1) + 1.5 * gsq(2)));
}

// adds two-loop anomalous dimension contribution to gii given matrix inputs
// g^Li_Lj = m_{ij} for LH fields
// g^Ei_Ej = m_{ji} for RH fields CHECKED: 23/5/02
void StandardModel::getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
                          DoubleMatrix & gQQ, DoubleMatrix & gDD,
                          DoubleMatrix & gUU, double & gH1H1, double &
                          gH2H2, sBrevity & a) const
{
   // For calculational brevity
   DoubleMatrix &dt = a.dt, &ut = a.ut, &u2 = a.u2, &d2 = a.d2, &e2 = a.e2,
                 &u2t = a.u2t, &d2t = a.d2t, &e2t = a.e2t, &u1 = a.u1, &d1 = a.d1;
   double &uuT = a.uuT, &ddT = a.ddT, &eeT = a.eeT;
   DoubleVector &gsq = a.gsq, &g4 = a.g4;

   // Everything gets the (1/16pi^2)^2 factor at the bottom
   DoubleMatrix ee(3, 3), ll(3, 3), qq(3, 3), dd(3, 3), uu(3, 3);

   // Two-loop pure gauge anom dimensions
   double h1h1 = (3.75 * g4(2) + 2.07 * g4(1) + 0.9 * gsq(2) * gsq(1));
   double h2h2 = h1h1;
   ll = h1h1;
   ee = (234. * g4(1) / 25.0);
   qq = (-8.0 * g4(3) / 9.0 + 3.75 * g4(2) + 199.0 * g4(1) / 900.0 + 8.0 *
         gsq(3) * gsq(2) + 8 * gsq(3) * gsq(1) / 45.0 + 0.1 * gsq(1) *
         gsq(2));
   dd = (-8.0 * g4(3) / 9.0 + 202.0 / 225.0 * g4(1) + 32.0 / 45.0 *
         gsq(3) * gsq(1));
   uu = (-8.0 * g4(3) / 9.0 + 856.0 / 225.0 * g4(1) + 128.0 / 45.0 *
         gsq(3) * gsq(1));

   ll = ll + 1.2 * gsq(1) * e2;
   ee = ee + (6.0 * gsq(2) - 1.2 * gsq(1)) * e2t;
   qq = qq + 0.4 * gsq(1) * (d2 + 2.0 * u2);
   dd = dd + (6.0 * gsq(2) + 0.4 * gsq(1)) * d2t;
   uu = uu + (6.0 * gsq(2) - 0.4 * gsq(1)) * u2t;

   h1h1 = h1h1 + (16 * gsq(3) - 0.4 * gsq(1)) * ddT + 1.2 * gsq(1) *
          eeT;
   h2h2 = h2h2 + (16 * gsq(3) + 0.8 * gsq(1)) * uuT;

   // Two-loop pure Yukawa contributions
   double s = (eeT + 3.0 * ddT), t = (d2 * u2).trace();

   ll = ll - (2.0 * e2 * e2 + s * e2);
   ee = ee - (2.0 * e2t * e2t + 2.0 * s * e2t);
   qq = qq - (2.0 * d2 * d2 + d2 * s + 2.0 * u2 * u2 + 3.0 * uuT * u2);
   dd = dd - (2.0 * d2t * d2t + 2.0 * (dt * u2 * d1) + 2 * s * d2t);
   uu = uu - (2.0 * u2t * u2t + 2.0 * (ut * d2 * u1) + 6.0 * uuT * u2t);
   h1h1 = h1h1 - (3.0 * (e2 * e2).trace() + 9.0 * (d2t * d2t).trace() +
                  3.0 * t);
   h2h2 = h2h2 - (9.0 * (u2 * u2).trace() + 3.0 * t);

   const static double twolp = 4.010149318236068e-5; // 1/(16 pi^2)^2

   gLL = gLL + twolp * ll;
   gEE = gEE + twolp * ee;
   gQQ = gQQ + twolp * qq;
   gDD = gDD + twolp * dd;
   gUU = gUU + twolp * uu;
   gH1H1 = gH1H1 + twolp * h1h1;
   gH2H2 = gH2H2 + twolp * h2h2;
}

// Outputs wave function renormalisation for SUSY parameters and gauge beta
// functions up to 2 loops.
void StandardModel::anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
                                DoubleMatrix & gQQ, DoubleMatrix & gUU,
                                DoubleMatrix & gDD, DoubleVector & dg,
                                double & dgdash_1, double & dg_11, double & gH1H1, double & gH2H2,
                                sBrevity & a)  const
{
   // Constants for gauge running
   static DoubleVector bBeta(3), cuBeta(3), cdBeta(3), ceBeta(3);
   static DoubleMatrix babBeta(3, 3);
   if (bBeta(1) < 1.0e-5) // Constants not set yet
      setSMBetas(babBeta, cuBeta, cdBeta, ceBeta, bBeta);

   //  sBrevity a contains all of the shortcutted matrices etc;
   a.calculate(SMu.display(), SMd.display(), SMe.display(), g.display());

   // For calculational brevity
   DoubleVector &g3 = a.g3;

// 1 loop contributions:
   if (displayLoops() > 0) {
      static const double oneO16Pisq = 1.0 / (16.0 * sqr(PI));

      getOneLpAnom(gEE, gLL, gQQ, gDD, gUU, gH1H1, gH2H2, a);
      dg = oneO16Pisq * g3 * bBeta;

   }
}

// Outputs derivatives vector y[n] for SUSY parameters: interfaces to
// integration routines
DoubleVector StandardModel::beta() const
{
   static sBrevity a;
   // cout <<"In susy beta(() " << endl;
   // calculate the derivatives
   static StandardModel ds;

   ds = beta(a);

   return ds.display(); // convert to a long vector
}



// r should be valid AT mt
void StandardModel::setDiagYukawas(const QedQcd & r, double vev)
{

   double v1, v2; // Higgs VEVs

   v1 = vev * cos(atan(displayTanb()));
   v2 = vev * sin(atan(displayTanb()));

   DoubleMatrix u1(3, 3), d1(3, 3), e1(3, 3);

   double invv2, invv1;
   invv2 = 1.0 / v2;
   invv1 = 1.0 / v1;
   u1(1, 1) = r.displayMass(mUp) * invv2;
   u1(2, 2) = r.displayMass(mCharm) * invv2;
   u1(3, 3) = r.displayMass(mTop) * invv2;

   d1(1, 1) = r.displayMass(mDown) * invv1;
   d1(2, 2) = r.displayMass(mStrange) * invv1;
   d1(3, 3) = r.displayMass(mBottom) * invv1;
   e1(1, 1) = r.displayMass(mElectron) * invv1;
   e1(2, 2) = r.displayMass(mMuon) * invv1;
   e1(3, 3) = r.displayMass(mTau) * invv1;

   setSMYukawaMatrix(YU, u1);
   setSMYukawaMatrix(YD, d1);
   setSMYukawaMatrix(YE, e1);
}

// mix = 0 for all mixing in downs...at present this is the only possibility.
// Takes diagonal quark Yukawa matrices and mixes them up according to the CKM
// matrix assuming:
// mix=2, all mixing is in down sector
// mix=1, all mixing is in up sector
void StandardModel::quarkMixing(const DoubleMatrix & CKM, int mix)
{
   switch (mix) {
   case 1:
      setSMYukawaMatrix(YU, CKM.transpose() * displayYukawaMatrix(YU) * CKM);
      break;
   case 2:
      setSMYukawaMatrix(YD, CKM * displayYukawaMatrix(YD) * CKM.transpose());
      break;

   default:
      ostringstream ii;
      ii << "Error. StandardModel::quarkMixing called with mix=" << mix;
      throw ii.str();
   }
}

void StandardModel::getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
                                  CKM, int mix, double vev)
{
   setDiagYukawas(r, vev);
   quarkMixing(CKM, mix);
}

// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, at
// present, only diagonal masses are handled.
void StandardModel::getMasses(QedQcd & r, double vev) const
{
   double v1, v2;
   v1 = vev * cos(atan(displayTanb()));
   v2 = vev * sin(atan(displayTanb()));

   DoubleMatrix u1(displayYukawaMatrix(YU)), d1(displayYukawaMatrix(YD)),
                e1(displayYukawaMatrix(YE));
   r.setMass(mUp, v2 * u1(1, 1));
   r.setMass(mCharm, v2 * u1(2, 2));
   r.setMass(mTop, v2 * u1(3, 3));
   r.setMass(mDown, v1 * d1(1, 1));
   r.setMass(mStrange, v1 * d1(2, 2));
   r.setMass(mBottom, v1 * d1(3, 3));
   r.setMass(mElectron, v1 * e1(1, 1));
   r.setMass(mMuon, v1 * e1(2, 2));
   r.setMass(mTau, v1 * e1(3, 3));
}

// Rotates to quark mass basis, returning the mixing matrices defined as
// yu_diag = vul yu vur^+
// yd_diag = vdl yd vdr^+
// All matrices should be 3 by 3
void StandardModel::diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr,
			    DoubleMatrix & vul, DoubleMatrix & vur) const {
  DoubleMatrix u(3, 3), v(3, 3);
  DoubleVector ydDiag(3), yuDiag(3);
  displayYukawaMatrix(YU).diagonalise(u, v, yuDiag);
  vul = u.transpose(); vur = v.transpose();

  displayYukawaMatrix(YD).diagonalise(u, v, ydDiag);
  vdl = u.transpose(); vdr = v.transpose();
}
