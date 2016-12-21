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

/** \file lowe.cpp
   - Project:     SOFTSUSY
   - Author:      Ben Allanach, Alexander Voigt
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305
   - Webpage:     http://hepforge.cedar.ac.uk/softsusy/
*/

#include "lowe.h"
#include "ew_input.hpp"
#include "error.hpp"
#include "utils.h"
#include "wrappers.hpp"

#include <iostream>
#include <cmath>

namespace softsusy {

namespace {

// Given a value of mt, and alphas(MZ), find alphas(mt) to 1 loops in qcd:
// it's a very good approximation at these scales, better than 10^-3 accuracy
double getAsmt(double mtop, double alphasMz, double mz) {
  using std::log;
  return alphasMz /
      (1.0 - 23.0 * alphasMz / (6.0 * M_PI) * log(mz / mtop));
}

// Input pole mass of top and alphaS(mt), outputs running mass mt(mt)
// including one-loop standard model correction only
double getRunMt(double poleMt, double asmt) {
  return poleMt / (1.0 + (4.0 / (3.0 * M_PI)) * asmt);
}

// Given pole mass and alphaS(MZ), returns running top mass -- one loop qcd
double getRunMtFromMz(double poleMt, double asMZ, double mz) {
  return getRunMt(poleMt, getAsmt(poleMt, asMZ, mz));
}

} // anonymous namespace

const std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> QedQcd_input_parmeter_names = {
   "alpha_em_MSbar_at_MZ",
   "alpha_s_MSbar_at_MZ",
   "GFermi",
   "MZ_pole", "MW_pole",
   "Mv1_pole", "Mv2_pole", "Mv3_pole",
   "MElectron_pole", "MMuon_pole", "MTau_pole",
   "MU_2GeV", "MS_2GeV", "MT_pole",
   "MD_2GeV", "mc_mc", "mb_mb"
};

Eigen::ArrayXd QedQcd::gaugeDerivs(double x, const Eigen::ArrayXd& y)
{
  set_scale(std::exp(x));
  setAlpha(ALPHA, y(0));
  setAlpha(ALPHAS, y(1));

  Eigen::ArrayXd dydx(2);
  dydx(0) = qedBeta();
  dydx(1) = qcdBeta();

  return dydx;
}

// SM beta functions for the gauge couplings, neglecting Yukawa
// contributions, from arXiv:1208.3357 [hep-ph].
Eigen::ArrayXd QedQcd::smGaugeDerivs(double x, const Eigen::ArrayXd& y)
{
  const double oneO4Pi = 1.0 / (4.0 * PI);
  const double a1 = y(0);
  const double a2 = y(1);
  const double a3 = y(2);
  const int nG = 3;

  set_scale(std::exp(x));

  Eigen::ArrayXd dydx(3);

  dydx(0) = oneO4Pi * a1 * a1 * (0.2 + 8.0 * nG / 3.0 + oneO4Pi * (0.36 * a1
    + 1.8 * a2 + nG * (38.0 * a1 / 15.0 + 1.2 * a2 + 88.0 * a3 / 15.0)));
  dydx(1) = oneO4Pi * a2 * a2 * (-43.0 / 3.0 + 8.0 * nG / 3.0 + oneO4Pi *
    (0.6 * a1 - 259.0 * a2 / 3.0 + nG * (0.4 * a1 + 98.0 * a2 / 3.0 + 8.0
    * a3)));
  dydx(2) = oneO4Pi * a3 * a3 * (-22.0 + 8.0 * nG / 3.0 + oneO4Pi * (-204.0
    * a3 + nG * (11.0 * a1 / 15.0 + 3.0 * a2 + 152.0 * a3 / 3.0)));

  return dydx;
}

QedQcd::QedQcd()
  : a(2)
  , mf(9)
  , input(static_cast<unsigned>(NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS))
  , mbPole(PMBOTTOM)
  , ckm()
  , pmns()
{
  set_number_of_parameters(11);
  mf(0) = MUP; mf(1) = MCHARM;
  mf(3) = MDOWN; mf(4) = MSTRANGE; mf(5) = MBOTTOM;
  mf(6) = MELECTRON; mf(7) = MMUON; mf(8) = MTAU;
  a(0) = ALPHAMZ;  a(1) = ALPHASMZ;
  mf(2) = getRunMtFromMz(PMTOP, ALPHASMZ, flexiblesusy::Electroweak_constants::MZ);
  input(alpha_em_MSbar_at_MZ) = ALPHAMZ;
  input(alpha_s_MSbar_at_MZ) = ALPHASMZ;
  input(MT_pole) = PMTOP;
  input(mb_mb) = MBOTTOM;
  input(MTau_pole) = MTAU;
  input(MMuon_pole) = MMUON;
  input(MElectron_pole) = MELECTRON;
  input(MW_pole) = flexiblesusy::Electroweak_constants::MW;
  input(MZ_pole) = flexiblesusy::Electroweak_constants::MZ;
  input(GFermi) = flexiblesusy::Electroweak_constants::gfermi;
  input(mc_mc) = MCHARM;
  input(MU_2GeV) = MUP;
  input(MD_2GeV) = MDOWN;
  input(MS_2GeV) = MSTRANGE;
  set_scale(flexiblesusy::Electroweak_constants::MZ);
  set_loops(3);
  set_thresholds(1);
}

Eigen::ArrayXd QedQcd::get() const
{
   Eigen::ArrayXd y(11);
   y(0) = a(0);
   y(1) = a(1);
   for (int i = 0; i < mf.size(); i++)
      y(i + 2) = mf(i);
   return y;
}

void QedQcd::set(const Eigen::ArrayXd& y)
{
   a(0) = y(0);
   a(1) = y(1);
   for (int i = 0; i < mf.size(); i++)
      mf(i) = y(i + 2);
}

Eigen::ArrayXd QedQcd::beta() const
{
   Eigen::ArrayXd dydx(11);
   dydx(0) = qedBeta();
   dydx(1) = qcdBeta();
   const auto y = massBeta();
   for (int i = 0; i < y.size(); i++)
      dydx(i + 2) = y(i);
   return dydx;
}

void QedQcd::runto_safe(double scale, double eps)
{
   try {
      run_to(scale, eps);
   } catch (...) {
      throw flexiblesusy::NonPerturbativeRunningQedQcdError(
         std::string("Non-perturbative running to Q = ")
         + flexiblesusy::ToString(scale)
         + " during determination of the SM(5) parameters.");
   }
}

//  Active flavours at energy mu
int QedQcd::flavours(double mu) const {
  int k = 0;
  // if (mu > mf(mTop - 1)) k++;
  if (mu > mf(mCharm - 1)) k++;
  if (mu > mf(mUp - 1)) k++;
  if (mu > mf(mDown - 1)) k++;
  if (mu > mf(mBottom - 1)) k++;
  if (mu > mf(mStrange - 1)) k++;
  return k;
}

ostream & operator <<(ostream &left, const QedQcd &m) {
  left << "mU: " << m.displayMass(mUp)
       << "  mC: " << m.displayMass(mCharm)
       << "  mt: " << m.displayMass(mTop)
       << "  mt^pole: " << m.displayPoleMt()
       << endl;
  left << "mD: " << m.displayMass(mDown)
       << "  mS: " << m.displayMass(mStrange)
       << "  mB: " << m.displayMass(mBottom)
       << "  mb(mb):  " << m.displayMbMb()
       << endl;
  left << "mE: " << m.displayMass(mElectron)
       << "  mM: " << m.displayMass(mMuon)
       <<  "  mT: " << m.displayMass(mTau)
       << "  mb^pole: " << m.displayPoleMb()
       << endl;
  left << "aE: " << 1.0 / m.displayAlpha(ALPHA)
       << "  aS: " << m.displayAlpha(ALPHAS)
       << "   Q: " << m.get_scale()
       << "  mT^pole: " << m.displayPoleMtau()
       << endl;
  left << "loops: " << m.get_loops()
       << "        thresholds: " << m.get_thresholds() << endl;

  return left;
}

//  returns qed beta function at energy mu < mtop
double QedQcd::qedBeta() const {
  double x;
  x = 24.0 / 9.0;
  if (get_scale() > mf(mCharm - 1)) x += 8.0 / 9.0;
  // if (get_scale() > mf(mTop - 1)) x += 8.0 / 9.0;
  if (get_scale() > mf(mBottom - 1)) x += 2.0 / 9.0;
  if (get_scale() > mf(mTau - 1)) x += 2.0 / 3.0;
  if (get_scale() > MW) x += -7.0 / 2.0;

  return (x * sqr(a(ALPHA - 1)) / PI);
}

//  next routine calculates beta function to 3 loops in qcd for The Standard
//  Model. Note that if quark masses are running, the number of active quarks
//  will take this into account. Returns beta
double QedQcd::qcdBeta() const {
  static const double INVPI = 1.0 / PI;
  const int quarkFlavours = flavours(get_scale());
  double qb0, qb1, qb2;
  qb0 = (11.0e0 - (2.0e0 / 3.0e0 * quarkFlavours)) / 4.0;
  qb1 = (102.0e0 - (38.0e0 * quarkFlavours) / 3.0e0) / 16.0;
  qb2 = (2.857e3 * 0.5 - (5.033e3 * quarkFlavours) / 18.0  +
         (3.25e2 * sqr(quarkFlavours) ) / 5.4e1) / 64;

  double qa0 = 0., qa1 = 0., qa2 = 0.;

  if (get_loops() > 0) qa0 = qb0 * INVPI;
  if (get_loops() > 1) qa1 = qb1 * sqr(INVPI);
  if (get_loops() > 2) qa2 = qb2 * sqr(INVPI) * INVPI;

  // add contributions of the one, two and three loop constributions resp.
  double beta;
  beta = -2.0 * sqr(displayAlpha(ALPHAS)) *
    (qa0 + qa1 * displayAlpha(ALPHAS) + qa2 *
     sqr(displayAlpha(ALPHAS)));
  return beta;
}

//(See comments for above function). returns a vector x(1..9) of fermion mass
//beta functions -- been checked!
Eigen::ArrayXd QedQcd::massBeta() const {
  static const double INVPI = 1.0 / PI, ZETA3 = 1.202056903159594;
  Eigen::ArrayXd x(9);

  // qcd bits: 1,2,3 loop resp.
  double qg1 = 0., qg2 = 0., qg3 = 0.;
  const int quarkFlavours = flavours(get_scale());
  if (get_loops() > 0) qg1 = INVPI;
  if (get_loops() > 1)
    qg2 = (202.0 / 3.0 - (20.0e0 * quarkFlavours) / 9.0) * sqr(INVPI) / 16.0;
  if (get_loops() > 2)
    qg3 = (1.249e3 - ((2.216e3 * quarkFlavours) / 27.0e0 +
                      1.6e2 * ZETA3 * quarkFlavours / 3.0e0) -
           140.0e0 * quarkFlavours * quarkFlavours / 81.0e0) * sqr(INVPI) *
      INVPI / 64.0;

  const double qcd = -2.0 * a(ALPHAS - 1) * (
     qg1  + qg2 * a(ALPHAS - 1) + qg3 * sqr(a(ALPHAS - 1)));
  const double qed = -a(ALPHA - 1) * INVPI / 2;

  for (int i = 0; i < 3; i++)   // up quarks
    x(i) = (qcd + 4.0 * qed / 3.0) * mf(i);
  for (int i = 3; i < 6; i++)   // down quarks
    x(i) = (qcd + qed / 3.0) * mf(i);
  for (int i = 6; i < 9; i++)   // leptons
    x(i) = 3.0 * qed * mf(i);

  // switch off relevant beta functions
  if (get_thresholds() > 0)
    for(int i = 0; i < x.size(); i++) {
      if (get_scale() < mf(i))
         x(i) = 0.0;
    }
  // nowadays, u,d,s masses defined at 2 GeV: don't run them below that
  if (get_scale() < 2.0)
     x(mUp - 1) = x(mDown - 1) = x(mStrange - 1) = 0.0;

  return x;
}

void QedQcd::runGauge(double x1, double x2)
{
  const double tol = 1.0e-5;
  Eigen::ArrayXd y(2);
  y(0) = displayAlpha(ALPHA);
  y(1) = displayAlpha(ALPHAS);

  flexiblesusy::Beta_function::Derivs derivs = [this] (double x, const Eigen::ArrayXd& y) {
     return gaugeDerivs(x, y);
  };

  call_rk(x1, x2, y, derivs, tol);

  setAlpha(ALPHA, y(0));
  setAlpha(ALPHAS, y(1));
}

// Done at pole mb: extracts running mb(polemb)
double QedQcd::extractRunningMb(double alphasMb) {
  double mbPole = displayPoleMb();

  if (get_scale() != mbPole) {
    ostringstream ii;
    ii << "QedQcd::extractRunningMb called at scale "
         << get_scale() << " instead of mbpole\n";
    throw flexiblesusy::SetupError(ii.str());
  }

  // Following is the MSbar correction from QCD, hep-ph/9912391 and ZPC48 673
  // (1990)
  double delta = 0.;
  if (get_loops() > 0) delta = delta + 4.0 / 3.0 * alphasMb / PI;
  if (get_loops() > 1)
    delta = delta + sqr(alphasMb / PI) *
      (10.1667 + (displayMass(mUp) + displayMass(mDown) +
                  displayMass(mCharm) + displayMass(mStrange)) / mbPole);
  if (get_loops() > 2)
    delta = delta + 101.45424 * alphasMb / PI * sqr(alphasMb / PI);

  double mbmb = mbPole * (1.0 - delta);

  return mbmb;
}

// Supposed to be done at mb(mb) -- MSbar, calculates pole mass
double QedQcd::extractPoleMb(double alphasMb) {

  if (get_scale() != displayMass(mBottom)) {
    ostringstream ii;
    ii << "QedQcd::extractPoleMb called at scale " << get_scale() <<
      " instead of mb(mb)\n";
    throw flexiblesusy::SetupError(ii.str());
  }

  // Following is the MSbar correction from QCD, hep-ph/9912391
  double delta = 0.0;
  if (get_loops() > 0) delta = delta + 4.0 / 3.0 * alphasMb / PI;
  if (get_loops() > 1) delta = delta + sqr(alphasMb / PI) *
    (9.2778 + (displayMass(mUp) + displayMass(mDown) + displayMass(mCharm) +
               displayMass(mStrange)) / mbPole);
  if (get_loops() > 2)
    delta = delta + 94.4182 * alphasMb / PI * sqr(alphasMb / PI);

  double mbPole = displayMass(mBottom) * (1.0 + delta);

  return mbPole;
}

// Calculates the running mass from the pole mass:
void QedQcd::calcRunningMb()
{
  const double tol = 1.0e-5;
  // Save initial object
  const auto saving(get());
  const double saveMu = get_scale();

  // Set arbitrarily low bottom mass to make sure it's included in the RGEs
  setMass(mBottom, 0.);
  run_to(displayPoleMb(), tol);
  const double mbAtPoleMb = extractRunningMb(displayAlpha(ALPHAS));
  setMass(mBottom, mbAtPoleMb);
  // Now, by running down to 1 GeV, you'll be left with mb(mb) since it will
  // decouple at this scale.
  run_to(1.0, tol);
  const double mbmb = displayMass(mBottom);

  // restore initial object
  set(saving);
  set_scale(saveMu);
  setMass(mBottom, mbmb);
}

// Calculates the pole mass from the running mass, which should be defined at
// mb
void QedQcd::calcPoleMb()
{
  const double alphasMZ = displayAlpha(ALPHAS);
  const double alphaMZ = displayAlpha(ALPHA);
  const double saveMu = get_scale();

  runGauge(get_scale(), displayMass(mBottom));
  const double poleMb = extractPoleMb(displayAlpha(ALPHAS));
  setPoleMb(poleMb);

  // Reset to erase numerical integration errors.
  setAlpha(ALPHAS, alphasMZ);
  setAlpha(ALPHA, alphaMZ);
  set_scale(saveMu);
}

// Takes QedQcd object created at MZ and spits it out at mt
void QedQcd::toMt()
{
  const double tol = 1.0e-5;

  setMass(mTop, getRunMtFromMz(displayPoleMt(), displayAlpha(ALPHAS), displayPoleMZ()));
  calcPoleMb();

  const double alphasMZ = displayAlpha(ALPHAS);
  const double alphaMZ = displayAlpha(ALPHA);
  const double mz = displayPoleMZ();

  runGauge(mz, 1.0);
  // Run whole lot up to pole top mass
  const double mt = this->displayPoleMt();
  run(1.0, mz, tol);

  // Reset alphas to erase numerical integration errors.
  setAlpha(ALPHAS, alphasMZ);
  setAlpha(ALPHA, alphaMZ);
  run(mz, mt, tol);
}

// Takes QedQcd object created at MZ and spits it out at MZ
void QedQcd::toMz()
{
  const double mt = input(MT_pole), as = a(ALPHAS - 1);
  setMass(mTop, getRunMtFromMz(mt, as, displayPoleMZ()));
  calcPoleMb();

  const double tol = 1.0e-5;
  const double alphasMZ = displayAlpha(ALPHAS);
  const double alphaMZ = displayAlpha(ALPHA);
  const double mz = displayPoleMZ();

  runGauge(mz, 1.0);
  run(1.0, mz, tol);
  // Reset alphas to erase numerical integration errors.
  setAlpha(ALPHAS, alphasMZ);
  setAlpha(ALPHA, alphaMZ);
}

/**
 * Calculates all running parameters in the SM w/o top quark at Q.
 * This function can be called multiple times, leading to the same
 * result (in contrast to toMz()).
 *
 * @param scale target renormalization scale
 * @param precision_goal precision goal
 * @param max_iterations maximum number of iterations
 */
void QedQcd::to(double scale, double precision_goal, unsigned max_iterations) {
   unsigned it = 0;
   bool converged = false;
   auto qedqcd_old(get()), qedqcd_new(get());
   const double running_precision = 0.1 * precision_goal;

   while (!converged && it < max_iterations) {
      // set alpha_i(MZ)
      runto_safe(displayPoleMZ(), running_precision);
      setAlpha(ALPHA, input(alpha_em_MSbar_at_MZ));
      setAlpha(ALPHAS, input(alpha_s_MSbar_at_MZ));

      // set mb(mb)
      runto_safe(displayMbMb(), running_precision);
      setMass(mBottom, displayMbMb());
      setPoleMb(extractPoleMb(displayAlpha(ALPHAS)));

      // set mc(mc)
      runto_safe(displayMcMc(), running_precision);
      setMass(mCharm, displayMcMc());

      // set mu, md, ms at 2 GeV
      runto_safe(2.0, running_precision);
      setMass(mUp, displayMu2GeV());
      setMass(mDown, displayMd2GeV());
      setMass(mStrange, displayMs2GeV());

      // set me, mm, ml at 2 GeV
      setMass(mElectron, displayPoleMel());
      setMass(mMuon, displayPoleMmuon());
      setMass(mTau, displayPoleMtau());

      // check convergence
      runto_safe(scale, running_precision);
      qedqcd_new = get();

      converged = flexiblesusy::MaxRelDiff(qedqcd_old, qedqcd_new) < precision_goal;

      qedqcd_old = qedqcd_new;

      it++;
   }

   // set alpha_i(MZ) on last time
   runto_safe(displayPoleMZ(), precision_goal);
   setAlpha(ALPHA, input(alpha_em_MSbar_at_MZ));
   setAlpha(ALPHAS, input(alpha_s_MSbar_at_MZ));

   runto_safe(scale, precision_goal);

   if (!converged && max_iterations > 0) {
      std::string msg =
         "Iteration to determine SM(5) parameters did not"
         " converge after " + std::to_string(max_iterations) +
         " iterations (precision goal: " + std::to_string(precision_goal)
         + ").";
      throw flexiblesusy::NoConvergenceError(max_iterations, msg);
   }
}

// This will calculate the three gauge couplings of the Standard Model at the
// scale m2.
// It's a simple one-loop calculation only and no
// thresholds are assumed. Range of validity is electroweak to top scale.
// alpha1 is in the GUT normalisation. sinth = sin^2 thetaW(Q) in MSbar
// scheme
Eigen::ArrayXd QedQcd::getGaugeMu(double m2, double sinth) const {
  using std::log;
  static const double INVPI = 1.0 / PI;
  Eigen::ArrayXd temp(3);

  const double aem = displayAlpha(ALPHA), m1 = get_scale();
  // Set alpha1,2 at scale m1 from data:
  double a1 = 5.0 * aem / (3.0 * (1.0 - sinth));
  double a2 = aem / sinth;

  const double mtpole = displayPoleMt();
  QedQcd oneset(*this);

  if (m1 < mtpole) {
    // Renormalise a1,a2 to threshold scale assuming topless SM with one
    // light Higgs doublet
    const double thresh = minimum(m2, mtpole);
    a1 = 1.0 / ( 1.0 / a1 + 4.0 * INVPI * 1.07e2 * log(m1 / thresh) / 2.4e2 );
    a2 = 1.0 / ( 1.0 / a2 - 4.0 * INVPI * 2.50e1 * log(m1 / thresh) / 4.8e1 );

    temp(0) = a1;
    temp(1) = a2;

    // calculate alphas(m2)
    if (m2 >= 1.0) {
       oneset.run_to(thresh);
    } else {
       oneset.run_to(1.0);
    }
    // Set alphas(m) to be what's already calculated.
    temp(2) = oneset.displayAlpha(ALPHAS);

    if (m2 > mtpole) {
      if (get_thresholds() > 0) {
        const double mtrun = oneset.displayMass(mTop);
        const double alphas_5f = oneset.displayAlpha(ALPHAS);
        const double alphas_sm = alphas_5f / (1.0 + INVPI * alphas_5f *
                                              log(mtrun / mtpole) / 3.0);
        oneset.setAlpha(ALPHAS, alphas_sm);
      }
      temp = oneset.runSMGauge(m2, temp);
    }
  } else {
    // Above the top threshold use SM RGEs only
    temp(0) = a1;
    temp(1) = a2;
    temp(2) = oneset.displayAlpha(ALPHAS);
    temp = oneset.runSMGauge(m2, temp);
  }

  return temp;
}

// Given the values of the SM gauge couplings alpha_i, i = 1, 2, 3, at
// the current scale, run to the scale end using SM RGEs.
// Range of validity is for scales greater than or equal to the
// top quark pole mass.
Eigen::ArrayXd QedQcd::runSMGauge(double end, const Eigen::ArrayXd& alphas)
{
  const double tol = 1.0e-5;
  const double start = get_scale();
  auto y = alphas;
  auto qedqcd(*this);

  flexiblesusy::Beta_function::Derivs derivs = [&qedqcd] (double x, const Eigen::ArrayXd& y) {
     return qedqcd.smGaugeDerivs(x, y);
  };

  call_rk(start, end, y, derivs, tol);

  return y;
}

void QedQcd::set_input(const Eigen::ArrayXd& pars)
{
   input = pars;
}

Eigen::ArrayXd QedQcd::display_input() const
{
   return input;
}

std::array<std::string, NUMBER_OF_LOW_ENERGY_INPUT_PARAMETERS> QedQcd::display_input_parameter_names()
{
   return QedQcd_input_parmeter_names;
}

bool operator ==(const QedQcd& a, const QedQcd& b)
{
   const double eps = 1e-10;

   return
      std::fabs(a.get_scale() - b.get_scale()) < eps &&
      std::fabs(a.get_loops() - b.get_loops()) < eps &&
      std::fabs(a.get_thresholds() - b.get_thresholds()) < eps &&
      std::fabs(a.displayAlpha(ALPHA) - b.displayAlpha(ALPHA)) < eps &&
      std::fabs(a.displayAlpha(ALPHAS) - b.displayAlpha(ALPHAS)) < eps &&
      std::fabs(a.displayMass(mUp) - b.displayMass(mUp)) < eps &&
      std::fabs(a.displayMass(mCharm) - b.displayMass(mCharm)) < eps &&
      std::fabs(a.displayMass(mTop) - b.displayMass(mTop)) < eps &&
      std::fabs(a.displayMass(mDown) - b.displayMass(mDown)) < eps &&
      std::fabs(a.displayMass(mStrange) - b.displayMass(mStrange)) < eps &&
      std::fabs(a.displayMass(mBottom) - b.displayMass(mBottom)) < eps &&
      std::fabs(a.displayMass(mElectron) - b.displayMass(mElectron)) < eps &&
      std::fabs(a.displayMass(mMuon) - b.displayMass(mMuon)) < eps &&
      std::fabs(a.displayMass(mTau) - b.displayMass(mTau)) < eps &&
      std::fabs(a.displayNeutrinoPoleMass(1) - b.displayNeutrinoPoleMass(1)) < eps &&
      std::fabs(a.displayNeutrinoPoleMass(2) - b.displayNeutrinoPoleMass(2)) < eps &&
      std::fabs(a.displayNeutrinoPoleMass(3) - b.displayNeutrinoPoleMass(3)) < eps &&
      std::fabs(a.displayPoleMt() - b.displayPoleMt()) < eps &&
      std::fabs(a.displayPoleMb() - b.displayPoleMb()) < eps &&
      std::fabs(a.displayPoleMtau() - b.displayPoleMtau()) < eps &&
      std::fabs(a.displayPoleMW() - b.displayPoleMW()) < eps &&
      std::fabs(a.displayPoleMZ() - b.displayPoleMZ()) < eps &&
      std::fabs(a.displayFermiConstant() - b.displayFermiConstant()) < eps;
}

} // namespace softsusy
