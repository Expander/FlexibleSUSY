
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SplitMSSM_tree_level_spectrum

#include <boost/test/unit_test.hpp>

#include "test_SplitMSSM.hpp"
#include "wrappers.hpp"
#include "conversion.hpp"
#include "SplitMSSM_two_scale_model.hpp"

using namespace flexiblesusy;

/// calculates Higgs mass in convention of arXiv:0705.1496
double calc_mh_tree(const SplitMSSM<Two_scale>& m)
{
   const double lambda = m.get_Lambdax();
   const double v = m.get_v();
   const double mh_tree = 2*lambda*Sqr(v);

   return Sqrt(mh_tree);
}

/// calculates Higgs mass in convention of arXiv:0705.1496
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

/// calculates Higgs mass in convention of arXiv:0705.1496
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

BOOST_AUTO_TEST_CASE( test_SplitMSSM_tree_level_masses )
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
