#include <fstream>
#include <iostream>
#include <iomanip>

#define private public

#include "CMSSM_mass_eigenstates.hpp"
#include "CMSSM_weinberg_angle.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"
#include "scan.hpp"

using namespace flexiblesusy;

void ensure_tree_level_ewsb(CMSSM_mass_eigenstates& model)
{
// ensure that the EWSB eqs. are satisfied (Drees p.222)
   const double vu   = model.get_vu();
   const double vd   = model.get_vd();
   const double gY   = model.get_g1() * Sqrt(0.6);
   const double g2   = model.get_g2();
   const double Mu   = model.get_Mu();
   const double BMu  = model.get_BMu();
   const double mHd2 = BMu*vu/vd - (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   const double mHu2 = BMu*vd/vu + (Sqr(gY) + Sqr(g2))*(Sqr(vd) - Sqr(vu))/8. - Sqr(Mu);
   model.set_mHd2(mHd2);
   model.set_mHu2(mHu2);
}

// template version of this function already included
// in test_CMSSMNoFV.hpp and test_CMSSM_two_loop_spectrum.cpp
void setup_CMSSM_const(CMSSM_mass_eigenstates& model, const CMSSM_input_parameters& input)
{
   const double ALPHASMZ = Electroweak_constants::alpha3;
   const double ALPHAMZ  = Electroweak_constants::aem;
   const double sinthWsq = Electroweak_constants::sinThetaW2;
   const double alpha1   = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2   = ALPHAMZ / sinthWsq;
   const double g1       = Sqrt(4 * Pi * alpha1);
   const double g2       = Sqrt(4 * Pi * alpha2);
   const double g3       = Sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta  = input.TanBeta;
   const double sinBeta  = sin(atan(tanBeta));
   const double cosBeta  = cos(atan(tanBeta));
   const double M12      = input.m12;
   const double m0       = input.m0;
   const double a0       = input.Azero;
   const double root2    = Sqrt(2.0);
   const double vev      = Electroweak_constants::vev;
   const double vu       = vev * sinBeta;
   const double vd       = vev * cosBeta;
   const double susyMu   = input.SignMu * 120.0;
   const double BMu      = Sqr(2.0 * susyMu);
   const double scale    = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = Electroweak_constants::mtoprun * root2 / (vev * sinBeta);
   Yd(2,2) = Electroweak_constants::mbrun   * root2 / (vev * cosBeta);
   Ye(2,2) = Electroweak_constants::mtau    * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   model.set_input_parameters(input);
   model.set_scale(scale);
   model.set_loops(1);
   model.set_thresholds(3);
   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(M12);
   model.set_MassG(M12);
   model.set_MassWB(M12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(Sqr(m0));
   model.set_TYu(a0 * Yu);
   model.set_TYd(a0 * Yd);
   model.set_TYe(a0 * Ye);
   model.set_Mu(susyMu);
   model.set_BMu(BMu);
   model.set_vu(vu);
   model.set_vd(vd);

   ensure_tree_level_ewsb(model);
   model.calculate_DRbar_masses();
}

void setup_CMSSM_const_non_3rd_gen(CMSSM_mass_eigenstates& model,
                                   const CMSSM_input_parameters& input)
{
   setup_CMSSM_const(model, input);

   const double ymu = 0.1;

   Eigen::Matrix<double,3,3> Ye = model.get_Ye();
   Ye(1,1) = ymu;
   model.set_Ye(Ye);
}

int main()
{
   CMSSM_mass_eigenstates model;
   CMSSM_input_parameters input;
   input.TanBeta = 10.;
   input.SignMu = 1;
   input.Azero = 0.;

   CMSSM_weinberg_angle::Sm_parameters sm_parameters;
   sm_parameters.fermi_constant = Electroweak_constants::gfermi;
   sm_parameters.mw_pole = Electroweak_constants::MW;
   sm_parameters.mz_pole = Electroweak_constants::MZ;
   sm_parameters.mt_pole = Electroweak_constants::PMTOP;

   CMSSM_weinberg_angle wein(&model, sm_parameters);

   const std::vector<double> range(float_range_log(100., exp((20*log(10000)-log(100))/19), 20));
   ofstream outf("/home/bachi/Programme/FlexibleSUSY/markus_tests/sinThetaW_variation.dat");

   std::cout << "Start writing to file" << '\n';
   outf << "# "
        << std::setw(12) << std::left << "m0=m12/GeV" << ' '
        << std::setw(12) << std::left << "var dvbsusy" << ' '
        << std::setw(12) << std::left << "var 2loopsm"
        << '\n';

   for (std::vector<double>::const_iterator it = range.begin(),
           end = range.end(); it != end; ++it) {
      input.m0 =  *it;
      input.m12 = *it;
      setup_CMSSM_const_non_3rd_gen(model, input);

      double sin_theta_1 = wein.calculate();
      wein.disable_dvbsusy();
      double sin_theta_2 = wein.calculate();
      wein.enable_dvbsusy();
      wein.set_number_of_loops(1);
      double sin_theta_3 = wein.calculate();
      wein.set_number_of_loops(2);

      outf << "  "
           << std::setw(12) << std::left << *it << ' '
           << std::setw(12) << std::left << 1 - sin_theta_2/sin_theta_1 << ' '
           << std::setw(12) << std::left << 1 - sin_theta_3/sin_theta_1
           << '\n';
   }

   std::cout << "End writing to file" << '\n';

   return 0;
}
