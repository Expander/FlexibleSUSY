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

#include "CMSSM_two_scale_model.hpp"
#include "weinberg_angle_pointer.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "numerics2.hpp"
#include "config.h"
#include "numerics.h"
#include "error.hpp"

#define ROOT2 Electroweak_constants::root2

namespace flexiblesusy {

namespace weinberg_angle {
/**
 * Sets the maximum number of iterations to 20, the number of loops to 2,
 * the precision goal to 1.0e-8, and the model pointer to the one
 * which is handed over as parameter.
 *
 * @param model_ pointer to the model for which the calculation shall be done
 */
Weinberg_angle_pointer::Weinberg_angle_pointer(const CMSSM<Two_scale>* model_)
   : number_of_iterations(20)
   , number_of_loops(2)
   , precision_goal(1.0e-8)
   , model(model_)
{
}

Weinberg_angle_pointer::~Weinberg_angle_pointer()
{
}

void Weinberg_angle_pointer::set_number_of_iterations(unsigned n)
{
   number_of_iterations = n;
}

void Weinberg_angle_pointer::set_number_of_loops(unsigned n)
{
   number_of_loops = n;
}

void Weinberg_angle_pointer::set_precision_goal(double p)
{
   precision_goal = p;
}

void Weinberg_angle_pointer::set_model_pointer(const CMSSM<Two_scale>* model_)
{
   model = model_;
}

/**
 * Calculates the DR-bar weak mixing angle \f$\sin\hat{\theta}_W\f$ as
 * defined in Eq. (C.3) from hep-ph/9606211 given the Fermi constant,
 * the Z-boson pole mass and the DR-bar electromagnetic coupling as input.
 *
 * The function throws an exception of type NoConvergenceError if the
 * iterative procedure to determine the weak mixing angle does not converge.
 *
 * @param rho_start initial guess for the rho-hat-parameter
 * @param sin_start initial guess for the sine of the weak mixing angle
 *
 * @return sine of the DR-bar weak mixing angle
 */
double Weinberg_angle_pointer::calculate(double rho_start, double sin_start)
{
   const double gY         = model->get_g1() * Sqrt(0.6);
   const double g2         = model->get_g2();
   const double e_drbar    = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(e_drbar) / (4.0 * Pi);
   const double mz_pole    = Electroweak_constants::MZ;
   const double scale      = model->get_scale();
   const double gfermi     = Electroweak_constants::gfermi;

   if (!is_equal(scale, mz_pole)) {
      WARNING("Weinberg_angle_pointer::calculate() called at scale "
              << scale << " != MZ_pole(" << mz_pole << ")");
   }

   unsigned iteration = 0;
   bool not_converged = true;
   double rho_old = rho_start, sin_old = sin_start;
   double rho_new = rho_start, sin_new = sin_start;

   while (not_converged && iteration < number_of_iterations) {
      const double deltaR = calculate_delta_r(rho_old, sin_old);

      double sin2thetasqO4 = Pi * alphaDRbar /
         (ROOT2 * Sqr(mz_pole) * gfermi * (1.0 - deltaR));

      if (sin2thetasqO4 >= 0.25)
         sin2thetasqO4 = 0.25;

      if (sin2thetasqO4 < 0.0)
         sin2thetasqO4 = 0.0;

      const double sin2theta = Sqrt(4.0 * sin2thetasqO4);
      const double theta = 0.5 * ArcSin(sin2theta);

      sin_new = Sin(theta);

      const double deltaRho = calculate_delta_rho(rho_old, sin_new);

      if (Abs(deltaRho) < 1.0)
         rho_new = 1.0 / (1.0 - deltaRho);
      else
         rho_new = 1.0;

      const double precision
         = Abs(rho_old / rho_new - 1.0) + Abs(sin_old / sin_new - 1.0);

      VERBOSE_MSG("Iteration step " << iteration
                  << ": prec=" << precision
                  << " dr=" << deltaR
                  << " drho=" << deltaRho
                  << " rho_new=" << rho_new
                  << " sin_new=" << sin_new);

      not_converged = precision >= precision_goal;

      rho_old = rho_new;
      sin_old = sin_new;
      iteration++;
   }

   if (not_converged)
      throw NoConvergenceError(number_of_iterations);

   return sin_new;
}

/**
 * Calculates the \f$\Delta\hat{\rho}\f$ corrections as defined in
 * Eqs. (C.4), (C.6) from hep-ph/9606211.
 *
 * @param rhohat rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\Delta\hat{\rho}\f$ as defined in (C.4) and (C.6) from hep-ph/9606211
 */
double Weinberg_angle_pointer::calculate_delta_rho(double rhohat, double sinThetaW)
{
   const double mz = Electroweak_constants::MZ;
   const double mw = Electroweak_constants::MW;
   const double mt = Electroweak_constants::PMTOP;
   const double mh = model->get_Mhh(0);
   const double g3 = model->get_g3();
   const double xt = 3.0 * Electroweak_constants::gfermi * Sqr(mt) *
                     ROOT2 * oneOver16PiSqr;

   const double sinb   = Sin(ArcTan(model->get_vu() / model->get_vd()));
   const double hmix_r = Sqr(model->get_ZH(0,1) / sinb);

   const double gY         = model->get_g1() * Sqrt(0.6);
   const double g2         = model->get_g2();
   const double e_drbar    = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(e_drbar) / (4.0 * Pi);

   const double mt_drbar    = model->get_MFu(2);
   const double pizztMZ     = Re(model->self_energy_VZ(mz));
   const double piwwtMW     = Re(model->self_energy_VWm(mw));
   double pizztMZ_corrected = pizztMZ;
   double piwwtMW_corrected = piwwtMW;
   if (model->get_thresholds() > 1) {
      pizztMZ_corrected =
         pizztMZ - calculate_self_energy_z_top(mz, mt_drbar)
                 + calculate_self_energy_z_top(mz, mt);
      piwwtMW_corrected =
         piwwtMW - calculate_self_energy_w_top(mw, mt_drbar)
                 + calculate_self_energy_w_top(mw, mt);
   }

   double deltaRho1Loop = 0.;
   if (number_of_loops > 0)
      deltaRho1Loop = pizztMZ_corrected / (rhohat * Sqr(mz)) -
         piwwtMW_corrected / Sqr(mw);

   double deltaRho2LoopSm = 0.;
   if (number_of_loops > 1) {
      deltaRho2LoopSm =
         alphaDRbar * Sqr(g3) * oneOver16PiSqr / (Pi * Sqr(sinThetaW)) *
         (-2.145 * Sqr(mt) / Sqr(mw) + 1.262 * log(mt / mz) - 2.24 -
          0.85 * Sqr(mz) / Sqr(mt)) +
         Sqr(xt) * hmix_r * rho_2(mh / mt) / 3.0;
   }

   const double deltaRho = deltaRho1Loop + deltaRho2LoopSm;

   return deltaRho;
}

/**
 * Calculates the \f$\Delta\hat{r}\f$ corrections as defined in
 * Eqs. (C.3), (C.5) from hep-ph/9606211.
 *
 * @param rhohat rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\Delta\hat{r}\f$ as defined in (C.3) and (C.5) from hep-ph/9606211
 */
double Weinberg_angle_pointer::calculate_delta_r(double rhohat, double sinThetaW)
{
   const double mz = Electroweak_constants::MZ;
   const double mw = Electroweak_constants::MW;
   const double mt = Electroweak_constants::PMTOP;
   const double mh = model->get_Mhh(0);
   const double g3 = model->get_g3();
   const double xt = 3.0 * Electroweak_constants::gfermi * Sqr(mt) *
                     ROOT2 * oneOver16PiSqr;
   const double outcos2 = 1.0 - Sqr(sinThetaW);

   const double sinb   = Sin(ArcTan(model->get_vu() / model->get_vd()));
   const double hmix_r = Sqr(model->get_ZH(0,1) / sinb);

   const double gY         = model->get_g1() * Sqrt(0.6);
   const double g2         = model->get_g2();
   const double e_drbar    = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(e_drbar) / (4.0 * Pi);

   const double mt_drbar    = model->get_MFu(2);
   const double pizztMZ     = Re(model->self_energy_VZ(mz));
   const double piwwt0      = Re(model->self_energy_VWm(0.));
   double pizztMZ_corrected = pizztMZ;
   double piwwt0_corrected  = piwwt0;
   if (model->get_thresholds() > 1) {
      pizztMZ_corrected =
         pizztMZ - calculate_self_energy_z_top(mz, mt_drbar)
                 + calculate_self_energy_z_top(mz, mt);
      piwwt0_corrected =
         piwwt0  - calculate_self_energy_w_top(0., mt_drbar)
                 + calculate_self_energy_w_top(0., mt);
   }

   double dvb = 0.;
   if (number_of_loops > 0)
      dvb = calculate_delta_vb(rhohat, sinThetaW);

   double deltaR1Loop = 0.;
   if (number_of_loops > 0)
      deltaR1Loop = rhohat * piwwt0_corrected / Sqr(mw) -
         pizztMZ_corrected / Sqr(mz) + dvb;

   double deltaR2LoopSm = 0.;
   if (number_of_loops > 1) {
      deltaR2LoopSm = alphaDRbar * Sqr(g3) * oneOver16PiSqr /
         (Pi * Sqr(sinThetaW) * outcos2) *
         (2.145 * Sqr(mt) / Sqr(mz) + 0.575 * log(mt / mz) - 0.224 -
          0.144 * Sqr(mz) / Sqr(mt)) -
         Sqr(xt) * hmix_r * rho_2(mh / mt) * (1.0 - deltaR1Loop) * rhohat / 3.0;
   }

   const double deltaR = deltaR1Loop + deltaR2LoopSm;

   return deltaR;
}

/**
 * Calculates the vertex, box and external wave-function
 * renormalizations \f$\delta_{\text{VB}}\f$ as given in
 * Eqs. (C.11)-(C.16), (C.20) from hep-ph/9606211.
 *
 * @param rhohat rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\delta_{\text{VB}}\f$ as defined in (C.11) from hep-ph/9606211
 */
double Weinberg_angle_pointer::calculate_delta_vb(double rhohat, double sinThetaW)
{
   const double deltaVbSm = calculate_delta_vb_sm(rhohat, sinThetaW);

   const double deltaVbSusy = calculate_delta_vb_susy(sinThetaW);

   const double deltaVb = deltaVbSm + deltaVbSusy;

   return deltaVb;
}

/**
 * Calculates the Standard Model vertex and box corrections
 * \f$\delta_{\text{VB}}^{\text{SM}}\f$ as given in Eq. (C.12) from
 * hep-ph/9606211 .
 *
 * @param rhohat rho-hat-parameter
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\delta_{\text{VB}}^{\text{SM}}\f$ as defined in (C.12)
 * from hep-ph/9606211
 */
double Weinberg_angle_pointer::calculate_delta_vb_sm(
   double rhohat, double sinThetaW)
{
   const double mz  = Electroweak_constants::MZ;
   const double mw  = Electroweak_constants::MW;
   const double cw2 = Sqr(mw / mz);
   const double sw2 = 1.0 - cw2;
   const double sinThetaW2 = Sqr(sinThetaW);
   const double outcos2    = 1.0 - sinThetaW2;

   const double gY         = model->get_g1() * Sqrt(0.6);
   const double g2         = model->get_g2();
   const double e_drbar    = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(e_drbar) / (4.0 * Pi);

   const double deltaVbSm = rhohat * alphaDRbar / (4.0 * Pi * sinThetaW2) *
      (6.0 + log(cw2) / sw2 *
       (3.5 - 2.5 * sw2 - sinThetaW2 * (5.0 - 1.5 * cw2 / outcos2)));

   return deltaVbSm;
}

/**
 * Calculates the SUSY vertex, box and external wave-function
 * renormalizations \f$\delta_{\text{VB}}^{\text{SUSY}}\f$ as given in
 * Eqs. (C.13)-(C.16), (C.20) from hep-ph/9606211 .
 *
 * @param sinThetaW sin(theta_W)
 *
 * @return \f$\delta_{\text{VB}}^{\text{SUSY}}\f$ as defined in (C.13)
 * from hep-ph/9606211
 */
double Weinberg_angle_pointer::calculate_delta_vb_susy(double sinThetaW)
{
   const double q          = model->get_scale();
   const double mz         = Electroweak_constants::MZ;
   const double gY         = model->get_g1() * Sqrt(0.6);
   const double g2         = model->get_g2();
   const double hmu        = Re(model->get_Ye(1,1));
   const double sinThetaW2 = Sqr(sinThetaW);
   const double outcos2    = 1.0 - sinThetaW2;

   const double e_drbar    = gY * g2 / Sqrt(Sqr(gY) + Sqr(g2));
   const double alphaDRbar = Sqr(e_drbar) / (4.0 * Pi);

   double msel  = 0.;
   double msmul = 0.;
   double msve  = 0.;
   double msvm  = 0.;
   const auto MSe = model->get_MSe();
   const auto ZE  = model->get_ZE();
   const auto MSv = model->get_MSv();
   const auto ZV  = model->get_ZV();
   for (int i = 0; i < decltype(MSe)::RowsAtCompileTime; i++) {
      msel  += AbsSqr(ZE(i,0))*MSe(i);
      msmul += AbsSqr(ZE(i,1))*MSe(i);
   }
   for (int i = 0; i < decltype(MSv)::RowsAtCompileTime; i++) {
      msve += AbsSqr(ZV(i,0))*MSv(i);
      msvm += AbsSqr(ZV(i,1))*MSv(i);
   }

   const Eigen::ArrayXd& mneut(model->get_MChi());
   const Eigen::ArrayXd& mch(model->get_MCha());
   const Eigen::MatrixXcd& n(model->get_ZN());
   const Eigen::MatrixXcd& u(model->get_UM());
   const Eigen::MatrixXcd& v(model->get_UP());

   const int dimN = mneut.rows();
   const int dimC = mch.rows();

   Eigen::VectorXd bPsi0NuSnul(Eigen::VectorXd::Zero(dimN)),
      bPsicNuSell(Eigen::VectorXd::Zero(dimC));
   Eigen::VectorXd bPsi0ESell(Eigen::VectorXd::Zero(dimN)),
      aPsicESnul(Eigen::VectorXd::Zero(dimC));
   Eigen::VectorXcd bChi0NuSnul, bChicNuSell;
   Eigen::VectorXcd bChi0ESell, aChicESnul;

   bPsi0NuSnul(0) = - gY / ROOT2;
   bPsi0NuSnul(1) =   g2 * ROOT2 * 0.5;
   bPsicNuSell(0) =   g2;
   bPsi0ESell(0)  = - gY / ROOT2;
   bPsi0ESell(1)  = - g2 * ROOT2 * 0.5;
   aPsicESnul(0)  =   g2;

   bChi0NuSnul = n * bPsi0NuSnul;
   bChicNuSell = u * bPsicNuSell;
   bChi0ESell  = n * bPsi0ESell;
   aChicESnul  = v.conjugate() * aPsicESnul;

   double deltaZnue = 0.0, deltaZe = 0.0;
   for (int i = 0; i < dimN; i++) {
      deltaZnue += - Sqr(Abs(bChi0NuSnul(i))) * b1(0.0, mneut(i), msve, q);
      deltaZe   += - Sqr(Abs(bChi0ESell(i)))  * b1(0.0, mneut(i), msel, q);
   }
   for (int i = 0; i < dimC; i++) {
      deltaZnue += - Sqr(Abs(bChicNuSell(i))) * b1(0.0, mch(i), msel, q);
      deltaZe   += - Sqr(Abs(aChicESnul(i)))  * b1(0.0, mch(i), msve, q);
   }

   Eigen::VectorXd bPsicNuSmul(Eigen::VectorXd::Zero(dimC));
   Eigen::VectorXd bPsi0MuSmul(Eigen::VectorXd::Zero(dimN)),
      aPsicMuSnul(Eigen::VectorXd::Zero(dimC));
   Eigen::VectorXcd bChicNuSmul;
   Eigen::VectorXcd bChi0MuSmul, aChicMuSnul;

   bPsicNuSmul(0) =   g2;
   bPsicNuSmul(1) = - hmu;
   bPsi0MuSmul(0) = - gY / ROOT2;
   bPsi0MuSmul(1) = - g2  * ROOT2 * 0.5;
   aPsicMuSnul(0) =   g2;
   aPsicMuSnul(1) = - hmu;

   bChicNuSmul = u * bPsicNuSmul;
   bChi0MuSmul = n * bPsi0MuSmul;
   aChicMuSnul = v.conjugate() * aPsicMuSnul;

   double deltaZnumu = 0.0, deltaZmu = 0.0;
   for (int i = 0; i < dimN; i++) {
      deltaZnumu += - Sqr(Abs(bChi0NuSnul(i))) * b1(0.0, mneut(i), msvm, q);
      deltaZmu   += - Sqr(Abs(bChi0MuSmul(i))) * b1(0.0, mneut(i), msmul, q);
   }
   for (int i = 0; i < dimC; i++) {
      deltaZnumu += - Sqr(Abs(bChicNuSmul(i))) * b1(0.0, mch(i), msmul, q);
      deltaZmu   += - Sqr(Abs(aChicMuSnul(i))) * b1(0.0, mch(i), msvm, q);
   }

   Eigen::MatrixXd aPsi0PsicW(Eigen::MatrixXd::Zero(dimN,dimC)),
      bPsi0PsicW(Eigen::MatrixXd::Zero(dimN,dimC));
   Eigen::MatrixXcd aChi0ChicW, bChi0ChicW;

   aPsi0PsicW(1, 0) = - g2;
   aPsi0PsicW(3, 1) =   g2 / ROOT2;
   bPsi0PsicW(1, 0) = - g2;
   bPsi0PsicW(2, 1) = - g2 / ROOT2;

   aChi0ChicW = n.conjugate() * aPsi0PsicW * v.transpose();
   bChi0ChicW = n * bPsi0PsicW * u.adjoint();

   const double b0_0_msel_msve  = b0(0.0, msel, msve, q);
   const double b0_0_msmul_msvm = b0(0.0, msmul, msvm, q);

   std::complex<double> deltaVE;
   for (int i = 0; i < dimC; i++) {
      for (int j = 0; j < dimN; j++) {
         const double c0_msel_mch_mneut = c0(msel, mch(i), mneut(j));
         const double c0_msve_mch_mneut = c0(msve, mch(i), mneut(j));
         const double b0_0_mch_mneut    = b0(0.0, mch(i), mneut(j), q);

         deltaVE += bChicNuSell(i) * Conj(bChi0ESell(j)) *
            (- ROOT2 / g2 * aChi0ChicW(j, i) * mch(i) * mneut(j) *
             c0_msel_mch_mneut + 1.0 / (ROOT2 * g2) * bChi0ChicW(j, i) *
             (b0_0_mch_mneut + Sqr(msel) * c0_msel_mch_mneut - 0.5));
         deltaVE += - aChicESnul(i) * bChi0NuSnul(j) *
            (- ROOT2 / g2 * bChi0ChicW(j, i) * mch(i) * mneut(j) *
             c0_msve_mch_mneut + 1.0 / (ROOT2 * g2) * aChi0ChicW(j, i) *
             (b0_0_mch_mneut + Sqr(msve) * c0_msve_mch_mneut - 0.5));
      }
   }

   for (int j = 0; j < dimN; j++) {
      deltaVE += 0.5 * Conj(bChi0ESell(j)) * bChi0NuSnul(j) *
         (b0_0_msel_msve + Sqr(mneut(j)) * c0(mneut(j), msel, msve) + 0.5);
   }

   std::complex<double> deltaVMu;
   for (int i = 0; i < dimC; i++) {
      for (int j = 0; j < dimN; j++) {
         const double c0_msmul_mch_mneut = c0(msmul, mch(i), mneut(j));
         const double c0_msvm_mch_mneut  = c0(msvm, mch(i), mneut(j));
         const double b0_0_mch_mneut     = b0(0.0, mch(i), mneut(j), q);

         deltaVMu += bChicNuSmul(i) * Conj(bChi0MuSmul(j)) *
            (- ROOT2 / g2 * aChi0ChicW(j, i) * mch(i) * mneut(j) *
             c0_msmul_mch_mneut + 1.0 / (ROOT2 * g2) * bChi0ChicW(j, i) *
             (b0_0_mch_mneut + Sqr(msmul) * c0_msmul_mch_mneut - 0.5));
         deltaVMu += - aChicMuSnul(i) * bChi0NuSnul(j) *
            (- ROOT2 / g2 * bChi0ChicW(j, i) * mch(i) * mneut(j) *
             c0_msvm_mch_mneut + 1.0 / (ROOT2 * g2) * aChi0ChicW(j, i) *
             (b0_0_mch_mneut + Sqr(msvm) * c0_msvm_mch_mneut - 0.5));
      }
   }

   for (int j = 0; j < dimN; j++) {
      deltaVMu += 0.5 * Conj(bChi0MuSmul(j)) * bChi0NuSnul(j) *
         (b0_0_msmul_msvm + Sqr(mneut(j)) * c0(mneut(j), msmul, msvm) + 0.5);
   }

   std::complex<double> a1;
   for(int i = 0; i < dimC; i++) {
      for(int j = 0; j < dimN; j++) {
         a1 += 0.5 * aChicMuSnul(i) * Conj(bChicNuSell(i)) *
            bChi0NuSnul(j) * bChi0ESell(j) * mch(i) * mneut(j) *
            d0(msel, msvm, mch(i), mneut(j));
         a1 += 0.5 * Conj(aChicESnul(i)) * bChicNuSmul(i) *
            Conj(bChi0NuSnul(j)) * Conj(bChi0MuSmul(j)) * mch(i) * mneut(j) *
            d0(msmul, msve, mch(i), mneut(j));
         a1 += bChicNuSmul(i) * Conj(bChicNuSell(i)) *
            Conj(bChi0MuSmul(j)) * bChi0ESell(j) *
            d27(msmul, msel, mch(i), mneut(j));
         a1 += Conj(aChicMuSnul(i)) * aChicESnul(i) *
            bChi0NuSnul(j) * Conj(bChi0NuSnul(j)) *
            d27(msvm, msve, mch(i), mneut(j));
      }
   }

   const double deltaVbSusy = oneOver16PiSqr *
      (- sinThetaW2 * outcos2 / (2.0 * Pi * alphaDRbar) * Sqr(mz) * a1.real() +
       deltaVE.real() + deltaVMu.real() +
       0.5 * (deltaZe + deltaZnue + deltaZmu + deltaZnumu));

   return deltaVbSusy;
}

/**
 * Calculates \f$\rho^{(2)}(r)\f$ as given in Eqs. (C.7)-(C.8) from
 * hep-ph/9606211.
 *
 * @param r ratio of Higgs mass over top quark mass
 *
 * @return \f$\rho^{(2)}(r)\f$
 */
double Weinberg_angle_pointer::rho_2(double r)
{
   const double Pi2 = Pi * Pi;

   if (r <= 1.9) {
      const double r2 = Sqr(r);
      return 19.0 - 16.5 * r + 43.0 / 12.0 * r2 + 7.0 / 120.0 * r2 * r -
         Pi * Sqrt(r) * (4.0 - 1.5 * r + 3.0 / 32.0 * r2 + r2 * r / 256.0) -
         Pi2 * (2.0 - 2.0 * r + 0.5 * r2) - Log(r) * (3.0 * r - 0.5 * r2);
   } else {
      const double rm1 = 1.0 / r, rm2 = Sqr(rm1), rm3 = rm2 * rm1,
         rm4 = rm3 * rm1, rm5 = rm4 * rm1;
      return Sqr(Log(r)) * (1.5 - 9.0 * rm1 - 15.0 * rm2 - 48.0 * rm3 -
                            168.0 * rm4 - 612.0 * rm5) -
         Log(r) * (13.5 + 4.0 * rm1 - 125.0 / 4.0 * rm2 - 558.0 / 5.0 * rm3 -
                   8307.0 / 20.0 * rm4 - 109321.0 / 70.0 * rm5) +
         Pi2 * (1.0 - 4.0 * rm1 - 5.0 * rm2 - 16.0 * rm3 -
                56.0 * rm4 - 204.0 * rm5) +
         49.0 / 4.0 + 2.0 / 3.0 * rm1 + 1613.0 / 48.0 * rm2 + 87.57 * rm3 +
         341959.0 / 1200.0 * rm4 + 9737663.0 / 9800.0 * rm5;
   }
}

/**
 * Calculates 1-loop top-quark contribution to Z boson self-energy.
 *
 * @param p momentum
 * @param mt top-quark mass
 *
 * @return 1-loop top-quark contribution to Z boson self-energy
 */
double Weinberg_angle_pointer::calculate_self_energy_z_top(double p, double mt)
{
   const double q   = model->get_scale();
   const double Nc  = 3.0;
   const double gY  = model->get_g1() * Sqrt(0.6);
   const double g2  = model->get_g2();
   const double gY2 = Sqr(gY);
   const double g22 = Sqr(g2);
   const double sw2 = gY2 / (gY2 + g22);
   const double cw2 = 1.0 - sw2;
   const double guL = 0.5 - 2.0 * sw2 / 3.0;
   const double guR = 2.0 * sw2 / 3.0;

   const double self_energy_z_top =
      Nc * Sqr(g2) / cw2 * oneOver16PiSqr *
      (hfn(p, mt, mt, q) * (Sqr(guL) + Sqr(guR)) -
       4.0 * guL * guR * Sqr(mt) * b0(p, mt, mt, q));

   return self_energy_z_top;
}

/**
 * Calculates 1-loop top-quark contribution to W boson self-energy.
 *
 * @param p momentum
 * @param mt top-quark mass
 *
 * @return 1-loop top-quark contribution to W boson self-energy
 */
double Weinberg_angle_pointer::calculate_self_energy_w_top(double p, double mt)
{
   const double q  = model->get_scale();
   const double mb = model->get_MFd(2);
   const double Nc = 3.0;
   const double g2 = model->get_g2();

   const double self_energy_w_top =
      0.5 * Nc * hfn(p, mt, mb, q) * Sqr(g2) * oneOver16PiSqr;

   return self_energy_w_top;
}

} // namespace weinberg_angle

} // namespace flexiblesusy
