
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_low_scale_constraint

#include <boost/test/unit_test.hpp>
#include "test.h"
#include <functional>
#include <Eigen/Dense>

#define private public

#include "MSSM_model.hpp"
#include "MSSM_low_scale_constraint.hpp"
#include "wrappers.hpp"

MSSM setup_MSSM()
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = 10;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = 100.0;
   const double m0 = 250.0;
   const double a0 = 50.0;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = 120.0;
   const double BMu = Sqr(2.0 * susyMu);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   MSSM m;
   m.set_scale(91);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);
   m.set_MassB(M12);
   m.set_MassG(M12);
   m.set_MassWB(M12);
   m.set_mq2(mm0);
   m.set_ml2(mm0);
   m.set_md2(mm0);
   m.set_mu2(mm0);
   m.set_me2(mm0);
   m.set_mHd2(Sqr(m0));
   m.set_mHu2(Sqr(m0));
   m.set_TYu(a0 * Yu);
   m.set_TYd(a0 * Yd);
   m.set_TYe(a0 * Ye);
   m.set_Mu(susyMu);
   m.set_BMu(BMu);
   m.set_vu(vu);
   m.set_vd(vd);

   return m;
}

DoubleVector calculate_gauge_couplings(MSSM model, MSSM_low_scale_constraint constraint, double scale)
{
   model.set_scale(scale);
   constraint.set_model(&model);
   constraint.calculate_DRbar_gauge_couplings();

   DoubleVector g(3);
   g(1) = model.get_g1();
   g(2) = model.get_g2();
   g(3) = model.get_g3();

   return g;
}

BOOST_AUTO_TEST_CASE( test_threshold_corrections )
{
   MSSM m(setup_MSSM());
   m.solve_ewsb(0);
   m.calculate_DRbar_parameters();

   MSSM_input_parameters input;
   MSSM_low_scale_constraint constraint(input);

   const double Q1 = constraint.get_scale();
   const double Q2 = 2. * Q1;
   const double oneOver16PiSqr = 1./(16. * M_PI * M_PI);
   const double gut_normalization = 3./5.;
   DoubleVector g_old(3);
   g_old(1) = m.get_g1();
   g_old(2) = m.get_g2();
   g_old(3) = m.get_g3();
   DoubleVector prefactor(3);
   for (int i = 1; i <= 3; i++)
      prefactor(i) = 1. / (oneOver16PiSqr * Power(g_old(i),3));

   const DoubleVector g_Q1(calculate_gauge_couplings(m, constraint, Q1));
   const DoubleVector g_Q2(calculate_gauge_couplings(m, constraint, Q2));

   const DoubleVector beta_numeric((g_Q1 - g_Q2) * prefactor * (1. / log(Q1/Q2)));
   DoubleVector beta_SM(3);
   beta_SM(1) = 41./6. * gut_normalization;
   beta_SM(2) = -19./6.;
   beta_SM(3) = -7.;
   DoubleVector beta_MSSM(3);
   beta_MSSM(1) = 11. * gut_normalization;
   beta_MSSM(2) = 1.;
   beta_MSSM(3) = -3.;

   // BOOST_CHECK_CLOSE_FRACTION(beta_numeric(1), beta_MSSM(1) - beta_SM(1), 0.5);
   // BOOST_CHECK_CLOSE_FRACTION(beta_numeric(2), beta_MSSM(2) - beta_SM(2), 0.5);
   BOOST_CHECK_CLOSE_FRACTION(beta_numeric(3), beta_MSSM(3) - beta_SM(3), 0.011);
}
