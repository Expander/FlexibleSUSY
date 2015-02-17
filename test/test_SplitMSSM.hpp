
#ifndef TEST_SplitMSSM_H
#define TEST_SplitMSSM_H

#include "wrappers.hpp"
#include "ew_input.hpp"
#include "SplitMSSM_two_scale_model.hpp"

using namespace flexiblesusy;

void setup_SplitMSSM_const(SplitMSSM<Two_scale>& m, const SplitMSSM_input_parameters& input)
{
   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double lambda = input.LambdaInput;
   const double kappa = 0.01;
   const double M1 = input.M1Input;
   const double M2 = input.M2Input;
   const double M3 = input.M3Input;
   const double g1u = input.g1uInput;
   const double g1d = input.g1dInput;
   const double g2u = input.g2uInput;
   const double g2d = input.g2dInput;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double Mu = input.MuInput;
   const double scale = Electroweak_constants::MZ;
   const double mu2 = Sqr(130.);

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / vev;
   Yd(2,2) = 2.9     * root2 / vev;
   Ye(2,2) = 1.77699 * root2 / vev;

   m.set_scale(scale);
   m.set_loops(1);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(Yu);
   m.set_Yd(Yd);
   m.set_Ye(Ye);

   m.set_Lambdax(lambda);
   m.set_g1u(g1u);
   m.set_g1d(g1d);
   m.set_g2u(g2u);
   m.set_g2d(g2d);
   m.set_MassB(M1);
   m.set_MassG(M3);
   m.set_MassWB(M2);
   m.set_Mu(Mu);
   m.set_mu2(mu2);
   m.set_v(vev);
}

void ensure_n_loop_ewsb(SplitMSSM<Two_scale>& m, int loop_level)
{
   const double precision = m.get_ewsb_iteration_precision();
   m.set_ewsb_loop_order(loop_level);
   m.solve_ewsb();

   if (loop_level == 1)
      BOOST_CHECK_SMALL(m.get_ewsb_eq_hh_1() - m.tadpole_hh().real(), precision);
}

#endif
