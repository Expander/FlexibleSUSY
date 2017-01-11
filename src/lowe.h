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

/** \file lowe.h
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
   \brief QedQcd object contains Standard Model quark and lepton
   masses. It integrates them using 3 loop qcd x 1 loop qed effective theory.
*/

#ifndef LOWE_H
#define LOWE_H

#include "betafunction.hpp"
#include "ckm.hpp"
#include "pmns.hpp"
#include <array>
#include <iosfwd>
#include <Eigen/Core>

namespace softsusy {

/// used to give order of quark masses stored
enum mass {mUp=1, mCharm, mTop, mDown, mStrange, mBottom, mElectron,
           mMuon, mTau};
/// order of gauge couplings stored in QedQcd
enum leGauge {ALPHA=1, ALPHAS};

enum QedQcd_input_parmeters : unsigned {
   alpha_em_MSbar_at_MZ,
   alpha_s_MSbar_at_MZ,
   GFermi,
   MZ_pole, MW_pole,
   Mv1_pole, Mv2_pole, Mv3_pole,
   MElectron_pole, MMuon_pole, MTau_pole,
   MU_2GeV, MS_2GeV, MT_pole,
   MD_2GeV, mc_mc, mb_mb,
   NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS
};

extern const std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> QedQcd_input_parmeter_names;

/// Quark and lepton masses and gauge couplings in QEDxQCD effective theory
class QedQcd: public flexiblesusy::Beta_function
{
private:
  Eigen::ArrayXd a; ///< gauge couplings
  Eigen::ArrayXd mf; ///< fermion running masses
  Eigen::ArrayXd input; ///< SLHA input parmeters
  double mbPole;    ///< pole masses of third family quarks
  flexiblesusy::CKM_parameters ckm; ///< CKM parameters (in the MS-bar scheme at MZ)
  flexiblesusy::PMNS_parameters pmns; ///< PMNS parameters (in the MS-bar scheme at MZ)

  double qedBeta() const;   ///< QED beta function
  double qcdBeta() const;   ///< QCD beta function
  Eigen::ArrayXd massBeta() const; ///< beta functions of masses
  Eigen::ArrayXd gaugeDerivs(double, const Eigen::ArrayXd&);
  Eigen::ArrayXd smGaugeDerivs(double, const Eigen::ArrayXd&);
  /// Does not run the masses, just gauge couplings from start to end
  void runGauge(double start, double end);
  Eigen::ArrayXd runSMGauge(double, const Eigen::ArrayXd&);
  void runto_safe(double, double); ///< throws if non-perturbative error occurs

  int flavours(double) const;  /// returns number of active flavours

  /// calculates pole bottom mass given alpha_s(Mb)^{MSbar} from running b mass
  double extractPoleMb(double asMb);

public:
  QedQcd();
  QedQcd(const QedQcd&) = default;
  QedQcd(QedQcd&&) = default;
  QedQcd& operator=(const QedQcd&) = default;
  QedQcd& operator=(QedQcd&&) = default;
  virtual ~QedQcd() {}

  // Beta_function interface
  virtual Eigen::ArrayXd get() const override;
  virtual void set(const Eigen::ArrayXd&) override;
  virtual Eigen::ArrayXd beta() const override;

  void setPoleMt(double mt) { input(MT_pole) = mt; }; ///< set pole top mass
  void setPoleMb(double mb) { mbPole = mb; }; ///< set pole bottom mass
  void setPoleMtau(double mtau) { input(MTau_pole) = mtau; }; ///< set pole tau mass
  void setPoleMmuon(double m) { input(MMuon_pole) = m; } ///< set pole muon mass
  void setPoleMel(double m) { input(MElectron_pole) = m; } ///< set pole electron mass
  void setMbMb(double mb)   { input(mb_mb) = mb;   }; ///< set mb(mb)
  void setMcMc(double mc)   { input(mc_mc) = mc;   }  ///< set mc(mc)
  void setMu2GeV(double mu) { input(MU_2GeV) = mu; } ///< set mu(2 GeV)
  void setMd2GeV(double md) { input(MD_2GeV) = md; } ///< set md(2 GeV)
  void setMs2GeV(double ms) { input(MS_2GeV) = ms; } ///< set ms(2 GeV)
  void setPoleMW(double mw) { input(MW_pole) = mw; } ///< set W boson pole mass
  void setPoleMZ(double mz) { input(MZ_pole) = mz; } ///< set Z boson pole mass
  /// sets a running quark mass
  void setMass(mass mno, double m) { mf(mno - 1) = m; }
  /// sets a neutrino pole mass
  void setNeutrinoPoleMass(int i, double m) { input(Mv1_pole + i - 1) = m; }
  /// sets QED or QCD structure constant
  void setAlpha(leGauge ai, double ap) { a(ai - 1) = ap; }
  /// set input value of alpha_em(MZ)
  void setAlphaEmInput(double a) { input(alpha_em_MSbar_at_MZ) = a; }
  /// set input value of alpha_s(MZ)
  void setAlphaSInput(double a) { input(alpha_s_MSbar_at_MZ) = a; }
  /// sets CKM parameters (in the MS-bar scheme at MZ)
  void setCKM(const flexiblesusy::CKM_parameters& ckm_) { ckm = ckm_; }
  /// sets PMNS parameters (in the MS-bar scheme at MZ)
  void setPMNS(const flexiblesusy::PMNS_parameters& pmns_) { pmns = pmns_; }
  /// sets Fermi constant
  void setFermiConstant(double gf) { input(GFermi) = gf; }
  /// sets all input parameters
  void set_input(const Eigen::ArrayXd&);

  /// Displays input parameters
  Eigen::ArrayXd displayInput() const { return input; }
  /// Display pole top mass
  double displayPoleMt() const { return input(MT_pole); };
  /// Display pole tau mass
  double displayPoleMtau() const { return input(MTau_pole); };
  /// Display pole muon mass
  double displayPoleMmuon() const { return input(MMuon_pole); };
  /// Display pole electron mass
  double displayPoleMel() const { return input(MElectron_pole); };
  /// Returns bottom "pole" mass
  double displayPoleMb() const { return mbPole; };
  /// Returns W boson pole mass
  double displayPoleMW() const { return input(MW_pole); }
  /// Returns Z boson pole mass
  double displayPoleMZ() const { return input(MZ_pole); }
  /// Returns a vector of running fermion masses
  const Eigen::ArrayXd& displayMass() const { return mf; }
  /// Returns a single running mass
  double displayMass(mass mno) const { return mf(mno - 1); }
  /// Returns a single neutrino pole mass
  double displayNeutrinoPoleMass(int i) const { return input(Mv1_pole + i - 1); }
  /// Returns a single gauge structure constant
  double displayAlpha(leGauge ai) const { return a(ai - 1); };
  /// Returns gauge structure constants
  Eigen::ArrayXd displayAlphas() const { return a; }
  /// Returns input value alpha_em(MZ)
  double displayAlphaEmInput() const { return input(alpha_em_MSbar_at_MZ); }
  /// Returns input value alpha_s(MZ)
  double displayAlphaSInput() const { return input(alpha_s_MSbar_at_MZ); }
  /// Returns Fermi constant
  double displayFermiConstant() const { return input(GFermi); }
  /// Obgligatory: returns vector of all running parameters
  /// returns vector of all input parameters
  Eigen::ArrayXd display_input() const;
  /// returns vector of all parameter names
  static std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> display_input_parameter_names();
  /// Returns mb(mb) MSbar
  double displayMbMb() const { return input(mb_mb); }
  /// Returns mc(mc) MSbar
  double displayMcMc() const { return input(mc_mc); }
  /// Returns mu(2 GeV)
  double displayMu2GeV() const { return input(MU_2GeV); }
  /// Returns md(2 GeV)
  double displayMd2GeV() const { return input(MD_2GeV); }
  /// Returns ms(2 GeV)
  double displayMs2GeV() const { return input(MS_2GeV); }
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

  /// Evolves object to MZ
  void toMz();
  /// Evolves object to given scale.  This implementation can be called multiple times
  void to(double, double tol = 1e-5, unsigned max_iterations = 20);
  /// This will calculate the three gauge couplings of the Standard Model at
  /// the scale m2.
  /// It's a simple one-loop calculation only and no
  /// thresholds are assumed. Range of validity is electroweak to top scale.
  // alpha1 is in the GUT normalisation. sinth = sin^2 thetaW(Q) in MSbar
  // scheme
  Eigen::ArrayXd getGaugeMu(double m2, double sinth) const;
};

/// Formatted output from QedQcd object
std::ostream & operator<<(std::ostream &, const QedQcd &);

bool operator ==(const QedQcd&, const QedQcd&);

} // namespace softsusy

#endif
