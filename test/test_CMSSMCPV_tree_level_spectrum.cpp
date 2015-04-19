
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_CMSSMCPV_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test.h"
#include "test_CMSSMCPV.hpp"
#include "CMSSM_two_scale_model.hpp"
#include "CMSSMCPV_two_scale_model.hpp"

#define COMPARE_MASS(p,dev) TEST_CLOSE(a.get_##p(), b.get_##p(), dev);

using namespace flexiblesusy;

void test_spectrum_equality(const CMSSMCPV<Two_scale>& a, const CMSSM<Two_scale>& b)
{
   COMPARE_MASS(MHpm, 1e-10);
   COMPARE_MASS(MChi, 1e-10);
   COMPARE_MASS(MCha, 1e-10);
   COMPARE_MASS(MGlu, 1e-10);
   COMPARE_MASS(MVP , 1e-10);
   COMPARE_MASS(MVZ , 1e-10);
   COMPARE_MASS(MVWm, 1e-10);
   COMPARE_MASS(MVG , 1e-10);
   COMPARE_MASS(MSu , 1e-10);
   COMPARE_MASS(MSd , 1e-10);
   COMPARE_MASS(MSe , 1e-10);
   COMPARE_MASS(MSv , 1e-10);
   COMPARE_MASS(MFu , 1e-10);
   COMPARE_MASS(MFd , 1e-10);
   COMPARE_MASS(MFe , 1e-10);
   COMPARE_MASS(MFv , 1e-10);

   // Higgses
   const Eigen::Array<double,2,1> Mhh(b.get_Mhh());
   const Eigen::Array<double,2,1> MAh(b.get_MAh());
   const Eigen::Array<double,4,1> Mhh_full(a.get_Mhh());
   Eigen::Array<double,4,1> Mhh_mixed;
   Mhh_mixed << Mhh(0), Mhh(1), MAh(0), MAh(1);

   std::sort(Mhh_mixed.data(), Mhh_mixed.data() + Mhh_mixed.size());

   TEST_EQUALITY(Mhh_full, Mhh_mixed);
}

void copy_parameters(const CMSSMCPV<Two_scale>& a, CMSSM<Two_scale>& b)
{
   b.set_scale (a.get_scale ());
   b.set_loops (a.get_loops ());
   b.set_g1    (a.get_g1    ());
   b.set_g2    (a.get_g2    ());
   b.set_g3    (a.get_g3    ());
   b.set_Yu    (a.get_Yu    ().real());
   b.set_Yd    (a.get_Yd    ().real());
   b.set_Ye    (a.get_Ye    ().real());
   b.set_MassB (a.get_MassB ().real());
   b.set_MassG (a.get_MassG ().real());
   b.set_MassWB(a.get_MassWB().real());
   b.set_mq2   (a.get_mq2   ().real());
   b.set_ml2   (a.get_ml2   ().real());
   b.set_md2   (a.get_md2   ().real());
   b.set_mu2   (a.get_mu2   ().real());
   b.set_me2   (a.get_me2   ().real());
   b.set_mHd2  (a.get_mHd2  ());
   b.set_mHu2  (a.get_mHu2  ());
   b.set_TYu   (a.get_TYu   ().real());
   b.set_TYd   (a.get_TYd   ().real());
   b.set_TYe   (a.get_TYe   ().real());
   b.set_Mu    (a.get_Mu    ().real());
   b.set_BMu   (a.get_BMu   ().real());
   b.set_vu    (a.get_vu    ());
   b.set_vd    (a.get_vd    ());
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_tree_level_spectrum )
{
   CMSSMCPV_input_parameters input;
   input.TanBeta = 10.;
   input.m0 = 125.;
   input.m12 = 200.;
   input.PhaseMu = 1.;
   input.Azero = 0.;
   input.etaInput = 0;

   CMSSMCPV<Two_scale> m1(input);
   CMSSM<Two_scale> m2;

   setup_CMSSMCPV(m1, input);
   m1.solve_ewsb_tree_level();
   copy_parameters(m1, m2);
   m1.calculate_DRbar_masses();
   m2.calculate_DRbar_masses();

   test_spectrum_equality(m1, m2);

   BOOST_REQUIRE(gErrors == 0);
   if (gErrors) {
      BOOST_FAIL("spectra are not equal");
      gErrors = 0;
   }
}

void check_goldstone_masses(const CMSSMCPV_input_parameters& input)
{
   CMSSMCPV<Two_scale> m;
   setup_CMSSMCPV(m, input);

   m.solve_ewsb_tree_level();
   m.calculate_DRbar_masses();
   m.reorder_DRbar_masses();

   BOOST_CHECK_CLOSE(m.get_Mhh(0) , m.get_MVZ() , 1e-10);
   BOOST_CHECK_CLOSE(m.get_MHpm(0), m.get_MVWm(), 1e-10);
}

BOOST_AUTO_TEST_CASE( test_CMSSMCPV_goldstone_boson_masses )
{
   for (unsigned m = 0; m < 100; m++) {
      for (unsigned e = 0; e < 100; e++) {
         CMSSMCPV_input_parameters input;
         input.TanBeta = 10.;
         input.m0 = 125.;
         input.m12 = 200.;
         input.PhaseMu = std::complex<double>(std::polar(1., 2*Pi*m/100));
         input.Azero = 100.;
         input.etaInput = 2*Pi*e/100;

         check_goldstone_masses(input);
      }
   }
}
