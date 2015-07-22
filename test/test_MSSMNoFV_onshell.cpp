
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSMNoFV_onshell

#include <boost/test/unit_test.hpp>

#include "test_CMSSMNoFV.hpp"
#include "wrappers.hpp"
#include "MSSMNoFV_onshell.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;
using namespace gm2os;

void setup(MSSMNoFV_mass_eigenstates& model)
{
   const double TanBeta = 10;
   const double SignMu = 1;
   const double m0 = 125;
   const double m12 = 500;
   const double Azero = 0;
   const double sinBeta = sin(atan(TanBeta));
   const double cosBeta = cos(atan(TanBeta));
   const double vev = 245.;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double root2 = Sqrt(2.0);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   model.set_scale(1000);
   model.set_loops(2);
   model.set_g1(Electroweak_constants::g1);
   model.set_g2(Electroweak_constants::g2);
   model.set_g3(Electroweak_constants::g3);
   model.set_Yu(Yu);
   model.set_Yd(Yd);
   model.set_Ye(Ye);
   model.set_MassB(m12);
   model.set_MassG(m12);
   model.set_MassWB(m12);
   model.set_mq2(mm0);
   model.set_ml2(mm0);
   model.set_md2(mm0);
   model.set_mu2(mm0);
   model.set_me2(mm0);
   model.set_mHd2(Sqr(m0));
   model.set_mHu2(-Sqr(m0));
   model.set_TYu(Azero * Yu);
   model.set_TYd(Azero * Yd);
   model.set_TYe(Azero * Ye);
   model.set_vu(vu);
   model.set_vd(vd);

   model.solve_ewsb_tree_level();
   model.calculate_DRbar_masses();

   model.solve_ewsb();
   model.do_calculate_sm_pole_masses(true);
   model.calculate_spectrum();
}

BOOST_AUTO_TEST_CASE( test_DRbar_os_conversion )
{
   MSSMNoFV_mass_eigenstates model;
   setup(model);

   MSSMNoFV_onshell osmodel(model);

   BOOST_CHECK_GT(std::abs(osmodel.get_MCha(0) - osmodel.get_physical().MCha(0)), 0.8);
   BOOST_CHECK_GT(std::abs(osmodel.get_MCha(1) - osmodel.get_physical().MCha(1)), 5.);
   BOOST_CHECK_GT(std::abs(osmodel.get_MChi(0) - osmodel.get_physical().MCha(0)), 5.);

   osmodel.convert_to_onshell();

   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MCha(0), osmodel.get_physical().MCha(0), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MCha(1), osmodel.get_physical().MCha(1), 1e-5);
   BOOST_CHECK_CLOSE_FRACTION(osmodel.get_MChi(0), osmodel.get_physical().MChi(0), 1e-4);

}
