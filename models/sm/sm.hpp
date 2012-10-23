
/** \file SMsusy.h
    - Project:     cE6SSM specgen
    - Author:      Peter Athron
*/

// This is a SM object created to do SM evolution of gauge and yukawas
// at one loop to allow me to do thresholds as Roman does.  I have
// ripped off the Essm object and I will not remove or change most of
// it as I am only interested in the beta functions. This menas object
// has lots of redundent stuff and functions which should not be
// called.


#ifndef SM_H
#define SM_H

#include "rge.h"
#include "susy.h"

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class QedQcd;
class sBrevity;

const static int numSusySMPars = 46;

class StandardModel: public RGE {
private:
   // Renamed these because the Standard Model Yukawas are different
   // from the susy ones (by factors of sin or cos beta) and I want to
   // avoid confusion
   DoubleMatrix SMu, SMd, SMe; ///< Yukawa matrices for ups, downs and leptons
   DoubleVector g; ///< Gauge couplings
   /// Bilinear Higgs superpotential parameter and ratio of Higgs VEVs,
   /// \f$ v_1/v_2 \f$
   double smu, tanb, hVev;
   //Peter:: New ESSM susy parameters
   DoubleVector lambda, kappa, MN_SUSY;
   double h_N, mu_0, gdash_1, g_11;

public:
   StandardModel(); ///< Constructor fills object with zeroes by default
   /// Constructor sets object to be equal to another
   StandardModel(const StandardModel &);
   /// Constructor given Yukawa matrices u,d,e, gauge couplings v, mu
   /// parameter=m, tan beta=tb, renormalisation scale MU, number of loops in
   /// RG evolution l and thresholds parameter t
   //  MssmSusy(const DoubleMatrix & u, const DoubleMatrix & d, const
   //DoubleMatrix & e, const DoubleVector & v, double m,
   //   double tb, double MU, int l, int t, double h);

   //Peter:: New constructor which takes new ESSM parameters

   StandardModel(const DoubleMatrix & SMu, const DoubleMatrix & SMd, const
          DoubleMatrix & SMe, const DoubleVector & v, double m,
          double tb, double MU, int l, int t, double h, DoubleVector lambda, DoubleVector kappa, DoubleVector MN_SUSY, double h_N, double mu_0, double gdash_1, double g_11);

   virtual ~StandardModel() {}; ///< Default destructor

   /// sets object to be equal to another
   const StandardModel & operator=(const StandardModel & s);
   /// sets object to be equal to another
   void setSusy(const StandardModel &s);

   /// Sets DRbar running Higgs vev
   void setHvev(double h);
   /// Copies Yukawa matrices and gauge couplings from s only
   void setSomePars(const StandardModel & s);
   /// Sets one element of a Yukawa matrix
   void setSMYukawaElement(yukawa, int, int, double);
   /// Sets whole Yukawa matrix
   void setSMYukawaMatrix(yukawa, const DoubleMatrix &);
   /// Set a single gauge coupling
   void setGaugeCoupling(int, double);
   /// Set all gauge couplings
   void setAllGauge(const DoubleVector  &);
   /// Sets superpotential mu parameter
   void setSusyMu(double);
   /// Sets all RGE parameters to elements of vector
   void set(const DoubleVector &);
   /// Sets tan beta
   void setTanb(double);

   //Peter:: setting new ESSM parameters.
   void seth_N(double);
   void setgash_1(double);
   void setg_11(double);
   void setmu_0(double);
   void setkappa(DoubleVector);
   void setlambda(DoubleVector);
   void setMN_SUSY(DoubleVector);
   /// Returns DRbar running Higgs vev
   double displayHvev() const;
   /// Returns whole object as a const
   inline StandardModel displaySusy() const;
   /// Returns a single Yukawa matrix element
   double displayYukawaElement(yukawa, int, int) const;
   /// Returns a whole Yukawa matrix
   DoubleMatrix displayYukawaMatrix(yukawa) const;
   /// Returns a single gauge coupling
   double displayGaugeCoupling(int) const;
   /// Returns all gauge couplings
   DoubleVector  displayGauge() const;
   /// Returns superpotential mu parameter
   double displaySusyMu() const;
   /// Returns all parameters as elements of a vector
   const DoubleVector display() const;
   /// Returns tan beta
   double displayTanb() const;
   double displayh_N() const;
   double displaymu_0() const;
   double displaygdash_1() const;
   double displayg_11() const;
   double displayMN_SUSY(int i) const;
   DoubleVector displaykappa() const;
   DoubleVector displaylambda() const;
   /// outputs object QedQcd & r valid at 1 GeV from SUSY data at mt, from
   /// diagonal elements of Yukawa couplings and Higgs VEV vev.
   void getMasses(QedQcd & r, double vev) const;
   /// This turns diagonal Yukawa couplings at MZ into CKM mixed ones
   /// Takes diagonal quark Yukawa matrices and mixes them up
   /// according to the CKM matrix assuming:
   /// mix=2, all mixing is in down sector
   /// mix=1, all mixing is in up sector
   void quarkMixing(const DoubleMatrix & CKM, int mix);
   /// Sets diagonal Yukawa couplings according to data in QedQcd input and
   /// Higgs VEV parameter vev=\f$v_1^2+v_2^2\f$
   void setDiagYukawas(const QedQcd &, double vev);
   /// Defines mixed Yukawa matrices from data input in form of CKM matrix and
   /// r, vev. If mix=2, all mixing is in down sector
   /// mix=1, all mixing is in up sector
   void getQuarkMixedYukawas(const QedQcd & r, const DoubleMatrix &
                             CKM, int mix, double vev);
   /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
   DoubleVector beta() const;
   /// Calculate beta functions of SUSY preserving parameters of RPC MSSM
   StandardModel beta(sBrevity &) const;
   /// Outputs one-loop anomlous dimensions gii given matrix inputs.
   /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
   /// respectively. Note that we use the convention (for matrices in terms of
   /// gamma's): gamma^Li_Lj = M_ij for LH fields and
   /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
   /// conjugates of the RH fields). a should already be defined.
   void getOneLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
                     DoubleMatrix & gQQ, DoubleMatrix & gDD,
                     DoubleMatrix & gUU, double & gH1H1, double &
                     gH2H2, sBrevity & a) const;
   /// Outputs two-loop anomlous dimensions gii given matrix inputs.
   /// for RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1 and H2
   /// respectively. Note that we use the convention (for matrices in terms of
   /// gamma's): gamma^Li_Lj = M_ij for LH fields and
   /// gamma^Rj_Ri = M_ij for RH fields (since they are really the complex
   /// conjugates of the RH fields). a should already be defined.
   void getTwoLpAnom(DoubleMatrix & gEE, DoubleMatrix & gLL,
                     DoubleMatrix & gQQ, DoubleMatrix & gDD,
                     DoubleMatrix & gUU, double & gH1H1, double &
                     gH2H2, sBrevity & a) const;
   /// Outputs wave function renormalisation for SUSY parameters and gauge beta
   /// functions up to 2 loops. Also calculates and outputs a.
   /// IO parameters: RH leptons, LH leptons, LH quarks, RH downs, RH ups, H1
   /// and H2 respectively.
   /// g^Li_Lj = m_{ij} for LH fields
   /// g^Ei_Ej = m_{ji} for RH fields
   //Peter:: edited to include new ESSM gage coupling
   void anomalousDimension(DoubleMatrix & gEE, DoubleMatrix & gLL,
                           DoubleMatrix & gQQ, DoubleMatrix & gUU,
                           DoubleMatrix & gDD, DoubleVector & dg,
                           double & dgdash_1, double & dg_11, double & gH1H1,
                           double & gH2H2, sBrevity & a) const;
   /// Rotates to quark mass basis, returning the mixing matrices defined as
   /// \f$(Y_U)_{diag}\f$ = vul \f$Y_U\f$ vur^+,
   /// \f$(Y_D)_{diag}\f$ = vdl \f$Y_D\f$ vdr^+
   /// All matrices should be 3 by 3
   void diagQuarkBasis(DoubleMatrix & vdl, DoubleMatrix & vdr,
                       DoubleMatrix & vul, DoubleMatrix & vur) const;
};
/// Formatted output
ostream & operator <<(ostream &, const StandardModel &);
/// Formatted input
istream & operator >>(istream &left, StandardModel &s);
/// Outputs beta function coefficients for MSSM gauge coupling evolution in
/// arguments.
void setSMBetas(DoubleMatrix &, DoubleVector  &, DoubleVector  &, DoubleVector
                &, DoubleVector  &);

inline StandardModel StandardModel::displaySusy() const
{
   return *this;
}

inline void StandardModel::setGaugeCoupling(int i, double f)
{
   g(i) = f;
}

inline void StandardModel::setAllGauge(const DoubleVector & v)
{
   if (v.displayStart() != 1 || v.displayEnd() != 3) {
      ostringstream ii;
      ii <<
         "Initialising SUSY params gauge function with vector NOT 1..3\n" <<
         v;
      throw ii.str();
   }
   g = v;
}

inline double StandardModel::displayHvev() const
{
   return hVev;
}

inline void StandardModel::setHvev(double h)
{
   hVev = h;
}
inline void StandardModel::setSusyMu(double f)
{
   smu = f;
}
inline void StandardModel::setTanb(double f)
{
   tanb = f;
}
inline void StandardModel::seth_N(double f)
{
   h_N = f;
}
inline  void  StandardModel::setgash_1(double f)
{
   gdash_1 = f;
}
inline  void  StandardModel::setg_11(double f)
{
   g_11 = f;
}
inline void  StandardModel::setmu_0(double f)
{
   mu_0 = f;
}
inline void  StandardModel::setkappa(DoubleVector v)
{
   kappa = v;
}
inline void  StandardModel::setlambda(DoubleVector v)
{
   lambda = v;
}
inline void  StandardModel::setMN_SUSY(DoubleVector v)
{
   MN_SUSY = v;
}

inline DoubleVector StandardModel::displayGauge() const
{
   return g;
}
inline double StandardModel::displayGaugeCoupling(int i) const
{
   return g.display(i);
}
inline double StandardModel::displaySusyMu() const
{
   return smu;
}

#endif
