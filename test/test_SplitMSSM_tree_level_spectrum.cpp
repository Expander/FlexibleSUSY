
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SplitMSSM_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SplitMSSM.hpp"
#include "test_FullMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "SplitMSSM_two_scale_model.hpp"
#include "FullMSSM_two_scale_model.hpp"

using namespace flexiblesusy;

/// calculates Higgs mass in convention of arXiv:0705.1496
double calc_mh_tree(const SplitMSSM<Two_scale>& m)
{
   const double lambda = m.get_Lambdax();
   const double v = m.get_v();
   const double mh_tree = 2*lambda*Sqr(v);

   return Sqrt(mh_tree);
}

/// calculates chargino mass matrix in convention of arXiv:0705.1496
Eigen::Matrix<double,2,2> calc_mass_matrix_Cha(const SplitMSSM<Two_scale>& m)
{
   const double M2 = m.get_MassWB();
   const double g2u = m.get_g2u();
   const double g2d = m.get_g2d();
   const double Mu = m.get_Mu();
   const double v = m.get_v();

   Eigen::Matrix<double,2,2> mcha_tree;
   mcha_tree <<
      M2   , g2u*v,
      g2d*v, Mu;

   return mcha_tree;
}

/// calculates neutralino mass matrix in convention of arXiv:0705.1496
Eigen::Matrix<double,4,4> calc_mass_matrix_Chi(const SplitMSSM<Two_scale>& m)
{
   const double M1 = m.get_MassB();
   const double M2 = m.get_MassWB();
   const double g1u = m.get_g1u();
   const double g1d = m.get_g1d();
   const double g2u = m.get_g2u();
   const double g2d = m.get_g2d();
   const double Mu = m.get_Mu();
   const double v = m.get_v();
   const double root2 = sqrt(2.0);

   Eigen::Matrix<double,4,4> mchi_tree;
   mchi_tree <<
      M1, 0 , -g1d*v/root2,  g1u*v/root2,
      0 , M2,  g2d*v/root2, -g2u*v/root2,
      0 , 0 ,  0          , -Mu         ,
      0 , 0 ,  0          ,  0          ;

   Symmetrize(mchi_tree);

   return mchi_tree;
}

BOOST_AUTO_TEST_CASE( test_SplitMSSM_tree_level_masses_convention )
{
   SplitMSSM_input_parameters input;
   input.LambdaInput = 0.1;
   input.g1uInput = 0.2;
   input.g1dInput = 0.3;
   input.g2uInput = 0.4;
   input.g2dInput = 0.5;
   input.M1Input = 100.;
   input.M2Input = 200.;
   input.M3Input = 300.;
   input.MuInput = 400.;
   SplitMSSM<Two_scale> m;
   setup_SplitMSSM_const(m, input);

   m.calculate_DRbar_masses();

   // Higgs mass
   const double Mhh_em = calc_mh_tree(m);
   const double Mhh = m.get_Mhh();

   BOOST_CHECK_CLOSE(Mhh_em, Mhh, 1.0e-10);

   // Gluino mass
   const double MGlu_em = m.get_MassG();
   const double MGlu = m.get_MGlu();

   BOOST_CHECK_CLOSE(MGlu_em, MGlu, 1.0e-10);

   // chargino mass matrix
   const Eigen::Matrix<double,2,2> M_Cha_em(calc_mass_matrix_Cha(m));
   const Eigen::Matrix<double,2,2> M_Cha(m.get_mass_matrix_Cha());

   for (int i = 0; i < 2; i++) {
      for (int k = 0; k < 2; k++) {
         BOOST_CHECK_CLOSE(M_Cha_em(i,k), M_Cha(i,k), 1.0e-10);
      }
   }

   // neutralino mass matrix
   const Eigen::Matrix<double,4,4> M_Chi_em(calc_mass_matrix_Chi(m));
   const Eigen::Matrix<double,4,4> M_Chi(m.get_mass_matrix_Chi());

   for (int i = 0; i < 4; i++) {
      for (int k = 0; k < 4; k++) {
         BOOST_CHECK_CLOSE(M_Chi_em(i,k), M_Chi(i,k), 1.0e-10);
      }
   }

   // top mass matrix
   // const Eigen::Matrix<double,3,3> M_top_em(calc_mass_matrix_top(m));
   // const Eigen::Matrix<double,3,3> M_top(m.get_mass_matrix_top());

   // for (int i = 0; i < 3; i++) {
   //    for (int k = 0; k < 3; k++) {
   //       BOOST_CHECK_EQUAL(M_top_em(i,k), M_top(i,k));
   //    }
   // }
}

void match_tree_level(const FullMSSM<Two_scale>& mssm,
                      SplitMSSM<Two_scale>& split)
{
   const double g1 = mssm.get_g1();
   const double g2 = mssm.get_g2();
   const double g3 = mssm.get_g3();
   const double GUT = sqrt(0.6);
   const double gY = GUT * g1;
   const double vu = mssm.get_vu();
   const double vd = mssm.get_vd();
   const double tan_beta = vu / vd;
   const double beta = ArcTan(tan_beta);
   const double sin_beta = Sin(beta);
   const double cos_beta = Cos(beta);
   const double cos_2beta = Cos(2*beta);

   const double lambda = 0.25 * (Sqr(g2) + Sqr(gY)) * Sqr(cos_2beta);
   const double vev = Sqrt(Sqr(vu) + Sqr(vd));
   const double g1u = gY * sin_beta;
   const double g1d = gY * cos_beta;
   const double g2u = g2 * sin_beta;
   const double g2d = g2 * cos_beta;

   const Eigen::Matrix<double,2,2> ZH(mssm.get_ZH());

   const double lambda_vertex =
      1./4. * (Sqr(gY) + Sqr(g2)) * Sqr(Sqr(ZH(0,0)) - Sqr(ZH(0,1)));

   BOOST_MESSAGE("lambda in terms of cos(2*beta)  = " << lambda);
   BOOST_MESSAGE("lambda in terms of cos(2*alpha) = " << lambda_vertex);

   split.set_scale(mssm.get_scale());
   split.set_loops(mssm.get_loops());
   split.set_g1(gY);
   split.set_g2(g2);
   split.set_g3(g3);
   split.set_Yu(mssm.get_Yu() * sin_beta);
   split.set_Yd(mssm.get_Yd() * cos_beta);
   split.set_Ye(mssm.get_Ye() * cos_beta);
   split.set_MassB(mssm.get_MassB());
   split.set_MassG(mssm.get_MassG());
   split.set_MassWB(mssm.get_MassWB());
   split.set_Mu(mssm.get_Mu());

   split.set_Lambdax(lambda);
   split.set_v(vev);
   split.set_g1u(g1u);
   split.set_g1d(g1d);
   split.set_g2u(g2u);
   split.set_g2d(g2d);
   // split.set_mu2(100); // unfixed
}

BOOST_AUTO_TEST_CASE( test_SplitMSSM_FullMSSM_tree_level_masses_convention )
{
   FullMSSM<Two_scale> mssm;
   FullMSSM_input_parameters input_full;

   input_full.m0 = 125.;
   input_full.m12 = 500.;
   input_full.TanBeta = 10.;
   input_full.SignMu = 1;
   input_full.Azero = 0.;

   setup_FullMSSM(mssm, input_full);
   mssm.solve_ewsb_tree_level();
   mssm.calculate_DRbar_masses();

   SplitMSSM<Two_scale> split;

   match_tree_level(mssm, split);
   split.solve_ewsb_tree_level(); // fixes mu2
   split.calculate_DRbar_masses();

   BOOST_CHECK_CLOSE(mssm.get_MVZ(), split.get_MVZ(), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MVWm(), split.get_MVWp(), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MAh(0), split.get_MAh(), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MHpm(0), split.get_MHp(), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFu(0), split.get_MFu(0), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFu(1), split.get_MFu(1), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFu(2), split.get_MFu(2), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFd(0), split.get_MFd(0), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFd(1), split.get_MFd(1), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFd(2), split.get_MFd(2), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFe(0), split.get_MFe(0), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFe(1), split.get_MFe(1), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MFe(2), split.get_MFe(2), 1.0e-10);

   BOOST_CHECK_CLOSE(mssm.get_MGlu(), split.get_MGlu(), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MCha(0), split.get_MCha(0), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MCha(1), split.get_MCha(1), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MChi(0), split.get_MChi(0), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MChi(1), split.get_MChi(1), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MChi(2), split.get_MChi(2), 1.0e-10);
   BOOST_CHECK_CLOSE(mssm.get_MChi(3), split.get_MChi(3), 1.0e-10);

   BOOST_CHECK_CLOSE(mssm.get_Mhh(0), split.get_Mhh(), 1.0e-10);
}
