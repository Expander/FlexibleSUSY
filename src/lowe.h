
/** \file lowe.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief QedQcd object contains Standard Model quark and lepton 
   masses. It integrates them using 3 loop qcd x 1 loop qed effective theory.

*/

#ifndef LOWE_H
#define LOWE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "def.h"
#include "utils.h"
#include "linalg.h"
#include "rge.h"
#include "ckm.hpp"
#include "pmns.hpp"

namespace softsusy {
const double MUP = 2.4e-3; ///< default running quark mass from PDG
const double MDOWN = 4.75e-3; ///< default running quark mass from PDG
const double MSTRANGE = 0.104; ///< default running quark mass from PDG
const double MCHARM = 1.27; ///< default running quark mass from PDG
const double MBOTTOM = 4.20; ///< default running quark mass from PDG
const double MTOP = 165.0; ///< default running quark mass from PDG
/// default pole lepton mass from PDG
const double MELECTRON = 5.10998902e-4; 
const double MMUON = 1.05658357e-1; ///< default pole lepton mass from PDG
const double MTAU = 1.77699; ///< default pole lepton mass from PDG
const double ALPHASMZ = 0.1184; ///< default running mass from PDG
const double ALPHAMZ = 1.0 / 127.916; ///< default running alpha(MZ) from PDG

const double PMTOP = 173.18; ///< default pole mass from CDF/D0 Run II 1207.1069
const double PMBOTTOM = 4.9; ///< default pole mass from PDG
/// default central values of CKM matrix elements from PDG 2006 in radians
const double THETA12CKM = 0.229206; ///< From Vus/Vud in global CKM fit, PDG
const double THETA13CKM = 0.003960; ///< From Vub in global CKM fit, PDG
const double THETA23CKM = 0.042223; ///< From Vcb/Vtb in global CKM fit, PDG

/// used to give order of quark masses stored
typedef enum {mUp=1, mCharm, mTop, mDown, mStrange, mBottom, mElectron,
	      mMuon, mTau} mass;
/// order of gauge couplings stored in QedQcd
typedef enum {ALPHA=1, ALPHAS} leGauge;

enum QedQcd_input_parmeters : unsigned { alpha_em_MSbar_at_MZ, GFermi,
      alpha_s_MSbar_at_MZ, MZ_pole, mb_mb, MT_pole, MTau_pole, MMuon_pole, Mv3_pole,
      MW_pole, ME_pole, Mv1_pole, MM_pole, Mv2_pole, MD, MU, MS, MC,
      NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS};

extern const char* QedQcd_input_parmeter_names[NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS];

/// Returns beta functions of alpha, alpha_s only
DoubleVector gaugeDerivs(double, const DoubleVector &);

/// Quark and lepton masses and gauge couplings in QEDxQCD effective theory
class QedQcd: public RGE 
{
private:
  DoubleVector a;   ///< gauge couplings
  DoubleVector mf;  ///< fermion running masses
  DoubleVector mnu; ///< neutrino pole masses
  double mtPole, mbPole; ///< pole masses of third family quarks
  double mbMb; ///< mb(mb) in the MSbar scheme with only QCD corrections
  double mtauPole; ///< tau pole mass
  double mmuonPole; ///< muon pole mass
  double mwPole; ///< W boson pole mass
  double mzPole; ///< Z boson pole mass
  double gfermi; ///< Fermi constant
  flexiblesusy::CKM_parameters ckm; ///< CKM parameters (in the MS-bar scheme at MZ)
  flexiblesusy::PMNS_parameters pmns; ///< PMNS parameters (in the MS-bar scheme at MZ)
  double qhiggs; ///< scale at which Mh is calculated

public:
  QedQcd(); ///< Initialises with default values defined in lowe.h
  QedQcd(const QedQcd &); ///< Initialises object with another
  const QedQcd& operator=(const QedQcd & m); ///< Sets two objects equal
  virtual ~QedQcd() {};
  
  void setPoleMt(double mt) { mtPole = mt; }; ///< set pole top mass
  void setPoleMb(double mb) { mbPole = mb; }; ///< set pole bottom mass
  void setPoleMtau(double mtau) { mtauPole = mtau; }; ///< set pole tau mass
  void setPoleMmuon(double m) { mmuonPole = m; }; ///< set pole tau mass
  void setMbMb(double mb)   { mbMb = mb;   }; ///< set mb(mb)
  void setPoleMW(double mw) { mwPole = mw; } ///< set W boson pole mass
  void setPoleMZ(double mz) { mzPole = mz; } ///< set Z boson pole mass
  /// sets a running quark mass
  void setMass(mass mno, double m) { mf(mno) = m; }; 
  /// sets a neutrino pole mass
  void setNeutrinoPoleMass(int i, double m) { mnu(i) = m; }
  /// sets QED or QCD structure constant
  void setAlpha(leGauge ai, double ap) { a(ai) = ap; }; 
  /// sets CKM parameters (in the MS-bar scheme at MZ)
  void setCKM(const flexiblesusy::CKM_parameters& ckm_) { ckm = ckm_; }
  /// sets PMNS parameters (in the MS-bar scheme at MZ)
  void setPMNS(const flexiblesusy::PMNS_parameters& pmns_) { pmns = pmns_; }
  /// sets Fermi constant
  void setFermiConstant(double gf) { gfermi = gf; }
  /// For exporting beta functions to Runge-Kutta
  void set(const DoubleVector &); 
  /// sets all input parameters
  void set_input(const Eigen::ArrayXd&);
  void setQHiggs(double Q) { qhiggs = Q; }
  
  /// Display pole top mass
  double displayPoleMt() const { return mtPole; };
  /// Display pole tau mass
  double displayPoleMtau() const { return mtauPole; };
  /// Display pole muon mass
  double displayPoleMmuon() const { return mmuonPole; };
  /// Returns bottom "pole" mass
  double displayPoleMb() const { return mbPole; };
  /// Returns W boson pole mass
  double displayPoleMW() const { return mwPole; }
  /// Returns Z boson pole mass
  double displayPoleMZ() const { return mzPole; }
  /// Returns a vector of running fermion masses
  const DoubleVector & displayMass() const { return mf; };
  /// Returns a single running mass
  double displayMass(mass mno) const { return mf.display(mno); };
  /// Returns a single neutrino pole mass
  double displayNeutrinoPoleMass(int i) const { return mnu.display(i); }
  /// Returns a single gauge structure constant
  double displayAlpha(leGauge ai) const { return a.display(ai); };
  /// Returns Fermi constant
  double displayFermiConstant() const { return gfermi; }
  double displayQHiggs() const { return qhiggs; }
  /// Obgligatory: returns vector of all running parameters
  const DoubleVector display() const;
  /// returns vector of all input parameters
  Eigen::ArrayXd display_input() const;
  /// returns vector of all parameter names
  static std::vector<std::string> display_input_parameter_names();
  /// Returns mb(mb) MSbar
  double displayMbMb() const { return mbMb; }
  /// returns CKM parameters
  flexiblesusy::CKM_parameters displayCKM() const { return ckm; }
  /// Returns real CKM matrix
  Eigen::Matrix<double,3,3> get_real_ckm() const { return ckm.get_real_ckm(); }
  /// Returns complex CKM matrix
  Eigen::Matrix<std::complex<double>,3,3> get_complex_ckm() const { return ckm.get_complex_ckm(); }
  /// returns PMNS parameters
  flexiblesusy::PMNS_parameters displayPMNS() const { return pmns; }
  /// Returns real PMNS matrix
  Eigen::Matrix<double,3,3> get_real_pmns() const { return pmns.get_real_pmns(); }
  /// Returns complex PMNS matrix
  Eigen::Matrix<std::complex<double>,3,3> get_complex_pmns() const { return pmns.get_complex_pmns(); }
  
  int flavours(double) const;  /// returns number of active flavours
  
  double qedBeta() const;   ///< QED beta function
  double qcdBeta() const;   ///< QCD beta function
  void massBeta(DoubleVector &) const; ///< beta functions of masses
  /// Beta functions of both beta-functions and all MSbar masses
  DoubleVector beta() const; 
  
  /// Does not run the masses, just gauge couplings from start to end
  void runGauge(double start, double end);
  /// calculates pole bottom mass given alpha_s(Mb)^{MSbar} from running b mass
  double extractPoleMb(double asMb);
  /// Done at pole mb: extracts running mb(polemb)
  double extractRunningMb(double asMb);
  /// calculates running bottom mass given alpha_s(Mb)^{MSbar} from pole m_b
  void calcRunningMb();
  /// Calculates the pole mass from the running mass, which should be defined
  /// at mb
  void calcPoleMb();

  /// Evolves object to running top mass
  void toMt();
  /// Evolves object to MZ
  void toMz();
  /// Evolves object to given scale
  void to(double);
  /// This will calculate the three gauge couplings of the Standard Model at
  /// the scale m2.
  /// It's a simple one-loop calculation only and no
  /// thresholds are assumed. Range of validity is electroweak to top scale.
  // alpha1 is in the GUT normalisation. sinth = sin^2 thetaW(Q) in MSbar
  // scheme
  DoubleVector  getGaugeMu(const double m2, const
		     double sinth) const;
};

/// Input numbers into the object: by file stream
ostream & operator <<(ostream &, const QedQcd &);
/// Formatted output from QedQcd object
istream & operator >>(istream &left, QedQcd &m);

/// Reads in a QedQed-type object and returns it in oneset.
/// Call with fname "" if you want it to come from standard input
/// "massIn" is an example of a data initialisation file: 
void readIn(QedQcd & oneset, const char fname[80]); 
/// Input pole mass of top and alphaS(mt), outputs running mass mt(mt)
/// including one-loop standard model correction only
double getRunMt(double poleMt, double asmt);
/// Given a value of mt, and alphas(MZ), find alphas(mt) to 1 loops in qcd:
/// it's a very good approximation at these scales, better than 10^-3 accuracy
double getAsmt(double mtop, double alphasMz);
/// Given pole mass and alphaS(MZ), returns running top mass -- one loop qcd
double getRunMtFromMz(double poleMt, double asMZ);

inline QedQcd::QedQcd(const QedQcd &m)
  : RGE(), a(m.a), mf(m.mf), mnu(m.mnu), mtPole(m.mtPole), mbPole(m.mbPole), mbMb(m.mbMb), 
   mtauPole(m.mtauPole), mmuonPole(m.mmuonPole), mwPole(m.mwPole), mzPole(m.mzPole), gfermi(m.gfermi),
   ckm(m.ckm), pmns(m.pmns), qhiggs(0.)
{
  setPars(11); 
  setMu(m.displayMu());
  setLoops(m.displayLoops());
  setThresholds(m.displayThresholds());
}

/// Returns diagonal fermion mass matrices given input object r
void massFermions(const QedQcd & r, DoubleMatrix & mDon, 
		  DoubleMatrix & mUpq, DoubleMatrix & mEle);

bool operator ==(const QedQcd&, const QedQcd&);

} // namespace softsusy

#endif
