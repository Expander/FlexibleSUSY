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
#include "linalg.h"
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
   const DoubleVector& mneut, // tree.mnBpmz
   const ComplexMatrix& n,    // tree.nBpmz
   const DoubleVector& mch,   // tree.mchBpmz
   const ComplexMatrix& u,    // tree.uBpmz
   const ComplexMatrix& v     // tree.vBpmz
) const
{
  const double g       = g2;
  const double gp      = gY;
  const double costh   = mw_pole / mz_pole;
  const double cw2     = sqr(costh);
  const double sw2     = 1.0 - cw2;
  const double outcos  = sqrt(1.0 - sqr(sinThetaW));
  const double q       = scale;

  //PA: get the dimension of menut
  const int dimN =  mneut.displayEnd();

  const double deltaVbSm =
    rho * alphaDRbar / (4.0 * PI * sqr(sinThetaW)) *
    (6.0 + log(cw2) / sw2 *
     (3.5 - 2.5 * sw2 - sqr(sinThetaW) * (5.0 - 1.5 * cw2 / sqr(outcos))));

  DoubleVector bPsi0NuNul(dimN), bPsicNuSell(2);
  DoubleVector bPsi0ESell(dimN), aPsicESnul(2);
  ComplexVector bChi0NuNul(dimN), bChicNuSell(2);
  ComplexVector bChi0ESell(dimN), aChicESnul(2);

  bPsicNuSell(1) = g;
  bPsi0NuNul(2) = root2 * g * 0.5;
  bPsi0NuNul(1) = - gp / root2;
  aPsicESnul(1) = g;
  bPsi0ESell(1) = -gp / root2;
  bPsi0ESell(2) = -g * root2 * 0.5;

  bChicNuSell = u * bPsicNuSell;
  bChi0ESell =  n * bPsi0ESell;
  bChi0NuNul = n * bPsi0NuNul;

  aChicESnul = v.complexConjugate() * aPsicESnul;

  double deltaZnue = 0.0, deltaZe = 0.0;
  int i; for(i=1; i<=dimN; i++) {
   if (i < 3) {
      deltaZnue = deltaZnue -
	sqr(bChicNuSell(i).mod()) * b1(0.0, mch(i), mselL, q);
      deltaZe = deltaZe -
	sqr(aChicESnul(i).mod()) * b1(0.0, mch(i), msnue, q);
    }
    deltaZnue = deltaZnue -
      sqr(bChi0NuNul(i).mod()) * b1(0.0, mneut(i), msnue, q);
    deltaZe = deltaZe -
      sqr(bChi0ESell(i).mod()) * b1(0.0, mneut(i), mselL, q);
  }

  DoubleVector bPsicNuSmul(2);
  DoubleVector bPsi0MuSmul(dimN), aPsicMuSnul(2);
  ComplexVector bChicNuSmul(2);
  ComplexVector bChi0MuSmul(dimN), aChicMuSnul(2);

  bPsicNuSmul(1) = g;
  bPsicNuSmul(2) = -hmu;
  aPsicMuSnul(1) = g;
  aPsicMuSnul(2) = -hmu;
  bPsi0MuSmul(1) = -gp / root2;
  bPsi0MuSmul(2) = -g * root2 * 0.5;

  bChicNuSmul = u * bPsicNuSmul;
  bChi0MuSmul =  n * bPsi0MuSmul;
  bChi0NuNul = n * bPsi0NuNul;
  aChicMuSnul = v.complexConjugate() * aPsicMuSnul;

  double deltaZnumu = 0.0, deltaZmu = 0.0;
  for(i=1; i<=dimN; i++) {
    if (i < 3) {
      deltaZnumu = deltaZnumu -
	sqr(bChicNuSmul(i).mod()) * b1(0.0, mch(i), msmuL, q);
      deltaZmu = deltaZmu -
	sqr(aChicMuSnul(i).mod()) * b1(0.0, mch(i), msnumu, q);
    }
    deltaZnumu = deltaZnumu -
      sqr(bChi0NuNul(i).mod()) * b1(0.0, mneut(i), msnumu, q);
    deltaZmu = deltaZmu -
      sqr(bChi0MuSmul(i).mod()) * b1(0.0, mneut(i), msmuL, q);
  }

  DoubleMatrix aPsi0PsicW(dimN, 2), bPsi0PsicW(dimN, 2), fW(dimN, 2), gW(dimN, 2);
  ComplexMatrix aChi0ChicW(dimN, 2), bChi0ChicW(dimN, 2);

  aPsi0PsicW(2, 1) = - g;
  bPsi0PsicW(2, 1) = - g;
  aPsi0PsicW(4, 2) = g / root2;
  bPsi0PsicW(3, 2) = -g / root2;

  /// These ought to be in physpars
  aChi0ChicW = n.complexConjugate() * aPsi0PsicW * v.transpose();
  bChi0ChicW = n * bPsi0PsicW * u.hermitianConjugate();

  Complex deltaVE = 0.0;
  int j; for(i=1; i<=2; i++)
    for(j=1; j<=dimN; j++) {
      deltaVE = deltaVE + bChicNuSell(i) * bChi0ESell(j).conj() *
	(- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(mselL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 bChi0ChicW(j, i) *
	 (b0(0.0, mch(i), mneut(j), q) + sqr(mselL) *
	  c0(mselL, mch(i), mneut(j)) - 0.5));
      deltaVE = deltaVE - aChicESnul(i) * bChi0NuNul(j) *
	(- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(msnue, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 aChi0ChicW(j, i) *
	 (b0(0.0, mch(i), mneut(j), q) + sqr(msnue) *
	  c0(msnue, mch(i), mneut(j)) - 0.5));
      if (i == 1)
	deltaVE = deltaVE +
	  0.5 * bChi0ESell(j).conj() * bChi0NuNul(j) *
	  (b0(0.0, mselL, msnue, q) + sqr(mneut(j)) *
	   c0(mneut(j), mselL, msnue) + 0.5);
    }

  Complex deltaVMu = 0.0;
  for(i=1; i<=2; i++)
    for(j=1; j<=dimN; j++) {
      deltaVMu = deltaVMu + bChicNuSmul(i) * bChi0MuSmul(j).conj() *
	(- root2 / g * aChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(msmuL, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 bChi0ChicW(j, i) *
	 (b0(0.0, mch(i), mneut(j), q) + sqr(msmuL) *
	  c0(msmuL, mch(i), mneut(j)) - 0.5));
      deltaVMu = deltaVMu - aChicMuSnul(i) * bChi0NuNul(j) *
	(- root2 / g * bChi0ChicW(j, i) * mch(i) * mneut(j) *
	 c0(msnumu, mch(i), mneut(j)) + 1.0 / (root2 * g) *
	 aChi0ChicW(j, i) *
	 (b0(0.0, mch(i), mneut(j), q) + sqr(msnumu) *
	  c0(msnumu, mch(i), mneut(j)) - 0.5));
      if (i == 1)
	deltaVMu = deltaVMu +
	  0.5 * bChi0MuSmul(j).conj() * bChi0NuNul(j) *
	  (b0(0.0, msmuL, msnumu, q) + sqr(mneut(j)) * c0(mneut(j), msmuL,
							  msnumu) + 0.5);
    }

  Complex a1(0.0, 0.0);
  for(i=1; i<=2; i++)
    for(j=1; j<=dimN; j++) {
      a1 = a1 + 0.5 * aChicMuSnul(i) * bChicNuSell(i).conj() *
	bChi0NuNul(j) * bChi0ESell(j) * mch(i) * mneut(j) *
	d0(mselL, msnumu, mch(i), mneut(j));
      a1 = a1 + 0.5 * aChicESnul(i).conj() * bChicNuSmul(i) *
	bChi0NuNul(j).conj() * bChi0MuSmul(j).conj() * mch(i) * mneut(j) *
	d0(msmuL, msnue, mch(i), mneut(j));
      a1 = a1 + bChicNuSmul(i) * bChicNuSell(i).conj() *
	bChi0MuSmul(j).conj() * bChi0ESell(j) *
	d27(msmuL, mselL, mch(i), mneut(j));
      a1 = a1 + aChicMuSnul(i).conj() * aChicESnul(i) *
	bChi0NuNul(j) * bChi0NuNul(j).conj() *
	d27(msnumu, msnue, mch(i), mneut(j));
    }

  const double deltaVbSusy =
    (-sqr(sinThetaW) * sqr(outcos) / (2.0 * PI * alphaDRbar) * sqr(mz_pole)
     * a1.real() + deltaVE.real() + deltaVMu.real() +
     0.5 * (deltaZe + deltaZnue + deltaZmu + deltaZnumu) ) /
    (sqr(PI) * 16.0);

  const double deltaVb = deltaVbSusy + deltaVbSm;

  return deltaVb;
}

} // namespace weinberg_angle

} // namespace flexiblesusy
