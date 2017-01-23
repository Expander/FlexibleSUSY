
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_effective_couplings

#include <boost/test/unit_test.hpp>

#include "effective_couplings.hpp"
#include "SM_mass_eigenstates.hpp"
#include "SM_effective_couplings.hpp"
#include "SM_two_scale_ewsb_solver.hpp"
#include "physical_input.hpp"
#include "standard_model.hpp"
#include "wrappers.hpp"

#include <complex>

using namespace flexiblesusy;
using namespace effective_couplings;

class Standard_model_tester {
public:
   Standard_model_tester(const SM_mass_eigenstates& model,
                         softsusy::QedQcd qedqcd_,
                         Physical_input physical_inputs_)
      : model(model), qedqcd(qedqcd_), physical_inputs(physical_inputs_)
      , rg_improve(true), include_qcd_corrections(true) {}
   ~Standard_model_tester() {}
   void do_rg_improve(bool flag) { rg_improve = flag; }
   void do_include_qcd_corrections(bool flag) { include_qcd_corrections = flag; }
   void set_model(const SM_mass_eigenstates& model_) { model = model_; }

   std::complex<double> get_eff_CphhVPVP();
   std::complex<double> get_eff_CphhVGVG();

private:
   SM_mass_eigenstates model;
   softsusy::QedQcd qedqcd;
   Physical_input physical_inputs;
   bool rg_improve;
   bool include_qcd_corrections;

   double number_of_active_flavours(double m) const;
   double qcd_scaling_factor(double m) const;
   void run_SM_gauge_couplings_to(double m);
};

double Standard_model_tester::number_of_active_flavours(double m) const
{
   if (m < qedqcd.displayMbMb()) {
      return 4.0;
   } else if (m < qedqcd.displayPoleMt()) {
      return 5.0;
   } else {
      return 6.0;
   }
}

double Standard_model_tester::qcd_scaling_factor(double m) const
{
   const double Nf = number_of_active_flavours(m);
   const double mtpole = qedqcd.displayPoleMt();
   const double l = Log(Sqr(m) / Sqr(mtpole));
   const double g3 = model.get_g3();

   const double nlo_qcd = 0.025330295910584444*(23.75 - 1.1666666666666667*Nf)*
      Sqr(g3);
   const double nnlo_qcd = 0.000641623890917771*Power(g3,4)*(370.1956513893174
      + 2.375*l + (-47.18640261449638 + 0.6666666666666666*l)*Nf +
      0.9017702481178881*Sqr(Nf));
   const double nnnlo_qcd = 0.000016252523020247696*Power(g3,6)*(467.683620788
      + 122.440972222*l + 10.9409722222*Sqr(l));

   return Sqrt(1.0 + nlo_qcd + nnlo_qcd + nnnlo_qcd);
}

void Standard_model_tester::run_SM_gauge_couplings_to(double m)
{
   using namespace standard_model;

   Standard_model sm;

   sm.set_loops(2);
   sm.set_thresholds(2);
   sm.set_physical_input(physical_inputs);

   sm.initialise_from_input(qedqcd);
   sm.run_to(m);

   model.set_g3(sm.get_g3());
}

// @note in the current mixed scheme, the LO expressions
// are not correctly reproduced, i.e. the ratios coupling * vev / mass
// do not reduce to 1 - this should be fixed
std::complex<double> Standard_model_tester::get_eff_CphhVPVP()
{
   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_pars(model.get());

   const double Mhh = model.get_physical().Mhh;

   run_SM_gauge_couplings_to(0.5 * Mhh);

   Eigen::Array<double,3,1> MFu(model.get_physical().MFu);
   const Eigen::Array<double,3,1> MFd(model.get_physical().MFd);
   const Eigen::Array<double,3,1> MFe(model.get_physical().MFe);
   const double MVWp = model.get_physical().MVWp;

   double qcd_fermion = 1.0;
   if (include_qcd_corrections) {
      qcd_fermion = 1.0 - Sqr(model.get_g3()) / (4.0 * Sqr(Pi));
   }

   std::complex<double> result = AS1(0.25 * Sqr(Mhh) / Sqr(MVWp));

   for (int i = 0; i < 3; ++i) {
      if (MFu(i) > Mhh) {
         result += 4.0 * qcd_fermion * AS12(0.25 * Sqr(Mhh) / Sqr(MFu(i))) / 3.0;
      } else {
        result += 4.0 * AS12(0.25 * Sqr(Mhh) / Sqr(MFu(i))) / 3.0;
      }
      if (MFd(i) > Mhh) {
         result += qcd_fermion * AS12(0.25 * Sqr(Mhh) / Sqr(MFd(i))) / 3.0;
      } else {
         result += AS12(0.25 * Sqr(Mhh) / Sqr(MFd(i))) / 3.0;
      }
      result += AS12(0.25 * Sqr(Mhh) / Sqr(MFe(i)));
   }

   result *= physical_inputs.get(Physical_input::alpha_em_0)
      * Sqrt(qedqcd.displayFermiConstant()) / (Power(2.0, 0.75) * Pi);

   model.set_scale(scale);
   model.set(saved_pars);

   return result;
}

std::complex<double> Standard_model_tester::get_eff_CphhVGVG()
{
   const double scale = model.get_scale();
   const Eigen::ArrayXd saved_pars(model.get());

   const double Mhh = model.get_physical().Mhh;

   run_SM_gauge_couplings_to(Mhh);
   model.calculate_DRbar_masses();

   Eigen::Array<double,3,1> MFu(model.get_physical().MFu);
   MFu(2) = qedqcd.displayPoleMt();
   const Eigen::Array<double,3,1> MFd(model.get_physical().MFd);

   double qcd_fermion = 1.0;
   if (include_qcd_corrections) {
      qcd_fermion = 1.0 + 2.0 * Sqr(model.get_g3()) / (3.0 * Sqr(Pi));
   }

   std::complex<double> result;

   for (int i = 0; i < 3; ++i) {
      result += qcd_fermion * AS12(0.25 * Sqr(Mhh) / Sqr(MFd(i)));
      result += qcd_fermion * AS12(0.25 * Sqr(Mhh) / Sqr(MFu(i)));
   }
   result *= std::complex<double>(0.75,0);

   if (include_qcd_corrections) {
      result *= qcd_scaling_factor(Mhh);
   }

   const double alpha_s = 0.07957747154594767*Sqr(model.get_g3());
   result *= alpha_s * Sqrt(qedqcd.displayFermiConstant())
      * Power(2.0, 0.25) / (3.0 * Pi);

   model.set(saved_pars);
   model.set_scale(scale);

   return result;
}

void set_test_model_parameters(SM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd)
{
   model.do_calculate_sm_pole_masses(true);


   const double v = 1.0 / Sqrt(qedqcd.displayFermiConstant() * Sqrt(2.0));

   const double g1 = Sqrt(5.0 / 3.0) * 3.58533449e-1;
   const double g2 = 6.47712483e-1;
   const double g3 = 1.16182994;

   Eigen::Matrix<double,3,3> Yu;
   Eigen::Matrix<double,3,3> Yd;
   Eigen::Matrix<double,3,3> Ye;

   Yu << 7.55665947e-6, 0.0, 0.0,
         0.0, 3.45513221e-3, 0.0,
         0.0, 0.0, Sqrt(2.0) * qedqcd.displayPoleMt() / v;

   Yd << 1.50336005e-5, 0.0, 0.0,
         0.0, 3.29156700e-4, 0.0,
         0.0, 0.0, 1.56338961e-2;

   Ye << 2.87235703e-6, 0.0, 0.0,
         0.0, 5.93912276e-4, 0.0,
         0.0, 0.0, 1.00093274e-2;

   const double Lambdax = 0.252805415;
   const double mu2 = 8555.44881e3;

   model.set_g1(g1);
   model.set_g2(g2);
   model.set_g3(g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_Lambdax(Lambdax);
   model.set_mu2(mu2);
   model.set_v(v);
   
   model.set_scale(173.34);

   // necessary to match LO expressions above
   model.set_pole_mass_loop_order(0);
   model.set_ewsb_loop_order(0);
   model.calculate_DRbar_masses();
   model.solve_ewsb_tree_level();
   model.calculate_spectrum();
}

BOOST_AUTO_TEST_CASE( test_LO_effective_couplings )
{
   softsusy::QedQcd qedqcd;
   Physical_input physical_inputs;

   SM_mass_eigenstates model;
   SM_ewsb_solver<Two_scale> ewsb_solver;
   model.set_ewsb_solver(&ewsb_solver);
   set_test_model_parameters(model, qedqcd);

   Standard_model_tester sm_tester(model, qedqcd, physical_inputs);
   sm_tester.do_rg_improve(false);
   sm_tester.do_include_qcd_corrections(false);

   const std::complex<double> expected_cphhVPVP = sm_tester.get_eff_CphhVPVP();
   const std::complex<double> expected_cphhVGVG = sm_tester.get_eff_CphhVGVG();

   SM_effective_couplings eff_cp(model, qedqcd, physical_inputs);
   eff_cp.do_run_couplings(false);
   eff_cp.do_include_qcd_corrections(false);
   eff_cp.calculate_effective_couplings();

   const std::complex<double> obtained_cphhVPVP = eff_cp.get_eff_CphhVPVP();
   const std::complex<double> obtained_cphhVGVG = eff_cp.get_eff_CphhVGVG();

   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cphhVPVP), Re(obtained_cphhVPVP), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cphhVPVP), Im(obtained_cphhVPVP), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Re(expected_cphhVGVG), Re(obtained_cphhVGVG), 1.0e-10);
   BOOST_CHECK_CLOSE_FRACTION(-Im(expected_cphhVGVG), Im(obtained_cphhVGVG), 1.0e-10);
}
