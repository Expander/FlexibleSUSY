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

#include "weinberg_angle.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "numerics.h"

namespace flexiblesusy {

namespace weinberg_angle {

using namespace softsusy;

Weinberg_angle::Weinberg_angle()
   : number_of_iterations(20)
   , precision_goal(1.0e-8)
   , alpha_em_drbar(Electroweak_constants::aem)
   , fermi_contant(Electroweak_constants::gfermi)
   , self_energy_z_at_mz(0.)
   , self_energy_z_at_0(0.)
   , self_energy_w_at_mw(0.)
   , rho_hat(0.)
{
}

Weinberg_angle::~Weinberg_angle()
{
}

void Weinberg_angle::set_number_of_iterations(unsigned n)
{
   number_of_iterations = n;
}

void Weinberg_angle::set_precision_goal(double p)
{
   precision_goal = p;
}

void Weinberg_angle::set_alpha_em_drbar(double a)
{
   alpha_em_drbar = a;
}

void Weinberg_angle::set_fermi_contant(double gfermi)
{
   fermi_contant = gfermi;
}

void Weinberg_angle::set_self_energy_z_at_mz(double s)
{
   self_energy_z_at_mz = s;
}

void Weinberg_angle::set_self_energy_z_at_0(double s)
{
   self_energy_z_at_0 = s;
}

void Weinberg_angle::set_self_energy_w_at_mw(double s)
{
   self_energy_w_at_mw = s;
}

double Weinberg_angle::get_rho_hat() const
{
   return rho_hat;
}

double Weinberg_angle::calculate() const
{
   return 0.;
}

double Weinberg_angle::calculate_delta_r(
   double scale,
   double rho,
   double sinThetaW,
   double mw_pole,
   double mz_pole,
   double alphaDRbar,
   double gY,                 // displayGaugeCoupling(1) * sqrt(0.6)
   double g2,                 // displayGaugeCoupling(2)
   double hmu,                // = displayYukawaElement(YE, 2, 2)
   double mselL,              // tree.me(1, 1)
   double msmuL,              // tree.me(1, 2)
   double msnue,              // tree.msnu(1)
   double msnumu,             // tree.msnu(2)
   const Eigen::ArrayXd& mneut, // tree.mnBpmz
   const Eigen::MatrixXcd& n,   // tree.nBpmz
   const Eigen::ArrayXd& mch,   // tree.mchBpmz
   const Eigen::MatrixXcd& u,   // tree.uBpmz
   const Eigen::MatrixXcd& v,    // tree.vBpmz
   double pizztMZ,
   double piwwt0,
   double mt,
   double gfermi,
   double g3,                 // displayGaugeCoupling(3)
   double tanBeta,
   double mh,
   double hmix12
) const
{
   const double outcos = Cos(ArcSin(sinThetaW));

   // 2-loop SM contribution
   const double sinb = Sin(ArcTan(tanBeta));

   const double xt = 3.0 * gfermi * Sqr(mt) / (8.0 * Sqr(Pi) * root2);

   const double dvb = calculate_delta_vb(
       scale,
       rho,
       sinThetaW,
       mw_pole,
       mz_pole,
       alphaDRbar,
       gY,
       g2,
       hmu,
       mselL,
       msmuL,
       msnue,
       msnumu,
       mneut,
       n,
       mch,
       u,
       v);

   const double mz = mz_pole;

   const double deltaR = rho * piwwt0 / Sqr(mw_pole) -
      pizztMZ / Sqr(mz) + dvb;

   // Dominant 2-loop SM term
   const double deltaR2LoopSm = alphaDRbar * Sqr(g3) /
      (16.0 * Sqr(Pi) * Pi * Sqr(sinThetaW) * Sqr(outcos)) *
      (2.145 * Sqr(mt) / Sqr(mz) + 0.575 * log(mt / mz) - 0.224
       - 0.144 * Sqr(mz) / Sqr(mt)) -
      Sqr(xt) * Sqr(hmix12) / Sqr(sinb) *
      rho_2(mh / mt) * (1.0 - deltaR) * rho / 3.0;

   const double deltaR_full = deltaR + deltaR2LoopSm;

   return deltaR_full;
}

double Weinberg_angle::calculate_delta_vb(
   double scale,
   double rho,
   double sinThetaW,
   double mw_pole,
   double mz_pole,
   double alphaDRbar,
   double gY,                 // displayGaugeCoupling(1) * sqrt(0.6)
   double g2,                 // displayGaugeCoupling(2)
   double hmu,                // = displayYukawaElement(YE, 2, 2)
   double mselL,              // tree.me(1, 1)
   double msmuL,              // tree.me(1, 2)
   double msnue,              // tree.msnu(1)
   double msnumu,             // tree.msnu(2)
   const Eigen::ArrayXd& mneut, // tree.mnBpmz
   const Eigen::MatrixXcd& n,   // tree.nBpmz
   const Eigen::ArrayXd& mch,   // tree.mchBpmz
   const Eigen::MatrixXcd& u,   // tree.uBpmz
   const Eigen::MatrixXcd& v    // tree.vBpmz
) const
{
  const double g       = g2;
  const double gp      = gY;
  const double costh   = mw_pole / mz_pole;
  const double cw2     = Sqr(costh);
  const double sw2     = 1.0 - cw2;
  const double sinThetaW2 = Sqr(sinThetaW);
  const double outcos  = Sqrt(1.0 - sinThetaW2);
  const double q       = scale;

  const int dimN =  mneut.rows();

  const double deltaVbSm =
     rho * alphaDRbar / (4.0 * Pi * sinThetaW2) *
     (6.0 + log(cw2) / sw2 *
      (3.5 - 2.5 * sw2 - sinThetaW2 * (5.0 - 1.5 * cw2 / Sqr(outcos))));

  Eigen::VectorXd bPsi0NuNul(Eigen::VectorXd::Zero(dimN)),
     bPsicNuSell(Eigen::VectorXd::Zero(2));
  Eigen::VectorXd bPsi0ESell(Eigen::VectorXd::Zero(dimN)),
     aPsicESnul(Eigen::VectorXd::Zero(2));
  Eigen::VectorXcd bChi0NuNul(Eigen::VectorXcd::Zero(dimN)),
     bChicNuSell(Eigen::VectorXcd::Zero(2));
  Eigen::VectorXcd bChi0ESell(Eigen::VectorXcd::Zero(dimN)),
     aChicESnul(Eigen::VectorXcd::Zero(2));

  bPsicNuSell(0) = g;
  bPsi0NuNul(1) = root2 * g * 0.5;
  bPsi0NuNul(0) = -gp / root2;
  aPsicESnul(0) = g;
  bPsi0ESell(0) = -gp / root2;
  bPsi0ESell(1) = -g * root2 * 0.5;

  bChicNuSell = u * bPsicNuSell;
  bChi0ESell =  n * bPsi0ESell;
  bChi0NuNul = n * bPsi0NuNul;

  aChicESnul = v.conjugate() * aPsicESnul;

  double deltaZnue = 0.0, deltaZe = 0.0;
  for (int i = 0; i < dimN; i++) {
     if (i < 2) {
        deltaZnue +=
           - Sqr(Abs(bChicNuSell(i))) * b1(0.0, mch(i), mselL, q);
        deltaZe +=
           - Sqr(Abs(aChicESnul(i))) * b1(0.0, mch(i), msnue, q);
     }
     deltaZnue +=
        - Sqr(Abs(bChi0NuNul(i))) * b1(0.0, mneut(i), msnue, q);
     deltaZe +=
        - Sqr(Abs(bChi0ESell(i))) * b1(0.0, mneut(i), mselL, q);
  }

  Eigen::VectorXd bPsicNuSmul(Eigen::VectorXd::Zero(2));
  Eigen::VectorXd bPsi0MuSmul(Eigen::VectorXd::Zero(dimN)),
     aPsicMuSnul(Eigen::VectorXd::Zero(2));
  Eigen::VectorXcd bChicNuSmul(Eigen::VectorXcd::Zero(2));
  Eigen::VectorXcd bChi0MuSmul(Eigen::VectorXcd::Zero(dimN)),
     aChicMuSnul(Eigen::VectorXcd::Zero(2));

  bPsicNuSmul(0) = g;
  bPsicNuSmul(1) = -hmu;
  aPsicMuSnul(0) = g;
  aPsicMuSnul(1) = -hmu;
  bPsi0MuSmul(0) = -gp / root2;
  bPsi0MuSmul(1) = -g * root2 * 0.5;

  bChicNuSmul = u * bPsicNuSmul;
  bChi0MuSmul =  n * bPsi0MuSmul;
  bChi0NuNul = n * bPsi0NuNul;
  aChicMuSnul = v.conjugate() * aPsicMuSnul;

  double deltaZnumu = 0.0, deltaZmu = 0.0;
  for(int i = 0; i < dimN; i++) {
     if (i < 2) {
        deltaZnumu +=
           - Sqr(Abs(bChicNuSmul(i))) * b1(0.0, mch(i), msmuL, q);
        deltaZmu +=
           - Sqr(Abs(aChicMuSnul(i))) * b1(0.0, mch(i), msnumu, q);
     }
     deltaZnumu +=
        - Sqr(Abs(bChi0NuNul(i))) * b1(0.0, mneut(i), msnumu, q);
     deltaZmu +=
        - Sqr(Abs(bChi0MuSmul(i))) * b1(0.0, mneut(i), msmuL, q);
  }

  Eigen::MatrixXd aPsi0PsicW(Eigen::MatrixXd::Zero(dimN,2)),
     bPsi0PsicW(Eigen::MatrixXd::Zero(dimN,2)),
     fW(Eigen::MatrixXd::Zero(dimN,2)),
     gW(Eigen::MatrixXd::Zero(dimN,2));
  Eigen::MatrixXcd aChi0ChicW(Eigen::MatrixXcd::Zero(dimN,2)),
     bChi0ChicW(Eigen::MatrixXcd::Zero(dimN,2));

  aPsi0PsicW(1, 0) = - g;
  bPsi0PsicW(1, 0) = - g;
  aPsi0PsicW(3, 1) = g / root2;
  bPsi0PsicW(2, 1) = -g / root2;

  /// These ought to be in physpars
  aChi0ChicW = n.conjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.adjoint();

  std::complex<double> deltaVE;
  for(int i = 0; i < 2; i++) {
     for(int j = 0; j < dimN; j++) {
        deltaVE += bChicNuSell(i) * Conj(bChi0ESell(j)) *
           (- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0(mselL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
            bChi0ChicW(j, i) *
            (b0(0.0, mch(i), mneut(j), q) + Sqr(mselL) *
             c0(mselL, mch(i), mneut(j)) - 0.5));
        deltaVE += - aChicESnul(i) * bChi0NuNul(j) *
           (- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0(msnue, mch(i), mneut(j)) + 1.0 / (root2 * g) *
            aChi0ChicW(j, i) *
            (b0(0.0, mch(i), mneut(j), q) + Sqr(msnue) *
             c0(msnue, mch(i), mneut(j)) - 0.5));
        if (i == 0) {
           deltaVE +=
              0.5 * Conj(bChi0ESell(j)) * bChi0NuNul(j) *
              (b0(0.0, mselL, msnue, q) + Sqr(mneut(j)) *
               c0(mneut(j), mselL, msnue) + 0.5);
        }
     }
  }

  std::complex<double> deltaVMu;
  for(int i = 0; i < 2; i++) {
     for(int j = 0; j < dimN; j++) {
        deltaVMu += bChicNuSmul(i) * Conj(bChi0MuSmul(j)) *
           (- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0(msmuL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
            bChi0ChicW(j, i) *
            (b0(0.0, mch(i), mneut(j), q) + Sqr(msmuL) *
             c0(msmuL, mch(i), mneut(j)) - 0.5));
        deltaVMu += - aChicMuSnul(i) * bChi0NuNul(j) *
           (- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
            c0(msnumu, mch(i), mneut(j)) + 1.0 / (root2 * g) *
            aChi0ChicW(j, i) *
            (b0(0.0, mch(i), mneut(j), q) + Sqr(msnumu) *
             c0(msnumu, mch(i), mneut(j)) - 0.5));
        if (i == 0) {
           deltaVMu +=
              0.5 * Conj(bChi0MuSmul(j)) * bChi0NuNul(j) *
              (b0(0.0, msmuL, msnumu, q) + Sqr(mneut(j)) * c0(mneut(j), msmuL,
                                                              msnumu) + 0.5);
        }
     }
  }

  std::complex<double> a1;
  for(int i = 0; i < 2; i++) {
     for(int j = 0; j < dimN; j++) {
        a1 += 0.5 * aChicMuSnul(i) * Conj(bChicNuSell(i)) *
           bChi0NuNul(j) * bChi0ESell(j) * mch(i) * mneut(j) *
           d0(mselL, msnumu, mch(i), mneut(j));
        a1 += 0.5 * Conj(aChicESnul(i)) * bChicNuSmul(i) *
           Conj(bChi0NuNul(j)) * Conj(bChi0MuSmul(j)) * mch(i) * mneut(j) *
           d0(msmuL, msnue, mch(i), mneut(j));
        a1 += bChicNuSmul(i) * Conj(bChicNuSell(i)) *
           Conj(bChi0MuSmul(j)) * bChi0ESell(j) *
           d27(msmuL, mselL, mch(i), mneut(j));
        a1 += Conj(aChicMuSnul(i)) * aChicESnul(i) *
           bChi0NuNul(j) * Conj(bChi0NuNul(j)) *
           d27(msnumu, msnue, mch(i), mneut(j));
     }
  }

  const double deltaVbSusy =
     (-sinThetaW2 * Sqr(outcos) / (2.0 * Pi * alphaDRbar) * Sqr(mz_pole)
      * a1.real() + deltaVE.real() + deltaVMu.real() +
      0.5 * (deltaZe + deltaZnue + deltaZmu + deltaZnumu) ) * oneOver16PiSqr;

  const double deltaVb = deltaVbSusy + deltaVbSm;

  return deltaVb;
}

double Weinberg_angle::rho_2(double r) const
{
   if (r <= 1.9) {
      return 19.0 - 16.5 * r + 43.0 * Sqr(r) / 12.0 + 7.0 / 120.0 * Sqr(r) * r -
         Pi * Sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * Sqr(r) + Sqr(r) * r /
                         256.0) - Sqr(Pi) * (2.0 - 2.0 * r + 0.5 * Sqr(r)) -
         Log(r) * (3.0 * r - 0.5 * Sqr(r));
   } else {
      const double rm1 = 1.0 / r, rm2 = Sqr(rm1), rm3 = rm2 * rm1,
         rm4 = rm3 * rm1, rm5 = rm4 * rm1;
      return Sqr(Log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 - 168.0
                            * rm4 - 612.0 * rm5) -
         Log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
                   8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5)
         + Sqr(Pi) * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
                      56.0 * rm4 - 204.0 * rm5)
         + 49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
         341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
   }
}

} // namespace weinberg_angle

} // namespace flexiblesusy
