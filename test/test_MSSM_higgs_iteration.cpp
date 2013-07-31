
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_higgs_iteration

#include <boost/test/unit_test.hpp>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>

#include "error.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"
#include "two_scale_solver.hpp"
#include "two_scale_running_precision.hpp"
#include "MSSM_model.hpp"
#include "MSSM_input_parameters.hpp"
#include "MSSM_high_scale_constraint.hpp"
#include "MSSM_susy_scale_constraint.hpp"
#include "MSSM_low_scale_constraint.hpp"
#include "MSSM_convergence_tester.hpp"
#include "MSSM_initial_guesser.hpp"
#include "test_MSSM.hpp"

#define SM(p) Electroweak_constants::p

double mHmZchi2(const gsl_vector* x, void* params)
{
   MSSM* model = static_cast<MSSM*>(params);

   const double vd = gsl_vector_get(x, 0);
   const double vu = gsl_vector_get(x, 1);

   model->set_vd(vd);
   model->set_vu(vu);

   model->calculate_DRbar_parameters();
   model->calculate_Mhh_pole_1loop();
   model->calculate_MVZ_pole_1loop();

   const double mH = model->get_physical().Mhh(1);
   const double mZ = model->get_physical().MVZ;

   return Sqr(SM(MZ) - mZ) + Sqr(SM(MH) - mH);
}

int minimize(MSSM* model)
{
   const gsl_multimin_fminimizer_type *T =
      gsl_multimin_fminimizer_nmsimplex2;
   gsl_multimin_fminimizer *s = NULL;
   gsl_vector *ss, *x;
   gsl_multimin_function minex_func;

   const size_t max_iter = 100;
   const size_t dimension = 2;
   const double precision = 1.0e-2;

   size_t iter = 0;
   int status;

   // Starting point
   x = gsl_vector_alloc(2);
   gsl_vector_set(x, 0, SM(vev));
   gsl_vector_set(x, 1, SM(vev));

   // Set initial step sizes to 1
   ss = gsl_vector_alloc(dimension);
   gsl_vector_set_all(ss, 1.0);

   // Initialize method and iterate
   minex_func.n = dimension;
   minex_func.f = mHmZchi2;
   minex_func.params = model;

   s = gsl_multimin_fminimizer_alloc(T, dimension);
   gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

   do
   {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
         break;

      const double size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size, precision);

      INFO("it " << iter << ", x = (" << gsl_vector_get(s->x, 0) << ","
           << gsl_vector_get(s->x, 1) << "), f() = " << s->fval
           << ", size = " << size);
   } while (status == GSL_CONTINUE && iter < max_iter);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
   BOOST_CHECK_LT(iter, max_iter);

   if (status == GSL_SUCCESS && iter < max_iter) {
      printf("converged at\n");
      INFO("it " << iter << ", x = (" << gsl_vector_get(s->x, 0) << ","
           << gsl_vector_get(s->x, 1) << "), f() = " << s->fval);
   }

   gsl_vector_free(x);
   gsl_vector_free(ss);
   gsl_multimin_fminimizer_free (s);

   return status;
}

BOOST_AUTO_TEST_CASE( test_MSSM_higgs_iteration )
{
   MSSM_input_parameters input;
   MSSM model;

   const double ALPHASMZ = 0.1176;
   const double ALPHAMZ = 1.0 / 127.918;
   const double sinthWsq = 0.23122;
   const double alpha1 = 5 * ALPHAMZ / (3 * (1 - sinthWsq));
   const double alpha2 = ALPHAMZ / sinthWsq;
   const double g1 = sqrt(4 * Pi * alpha1);
   const double g2 = sqrt(4 * Pi * alpha2);
   const double g3 = sqrt(4 * Pi * ALPHASMZ);
   const double tanBeta = input.TanBeta;
   const double sinBeta = sin(atan(tanBeta));
   const double cosBeta = cos(atan(tanBeta));
   const double M12 = input.m12;
   const double m0 = input.m0;
   const double a0 = input.Azero;
   const double root2 = sqrt(2.0);
   const double vev = 246.0;
   const double vu = vev * sinBeta;
   const double vd = vev * cosBeta;
   const double susyMu = input.SignMu * 120.0;
   const double BMu = Sqr(2.0 * susyMu);
   const double scale = Electroweak_constants::MZ;

   Eigen::Matrix<double,3,3> Yu(Eigen::Matrix<double,3,3>::Zero()),
      Yd(Eigen::Matrix<double,3,3>::Zero()),
      Ye(Eigen::Matrix<double,3,3>::Zero()),
      mm0(Eigen::Matrix<double,3,3>::Zero());
   Yu(2,2) = 165.0   * root2 / (vev * sinBeta);
   Yd(2,2) = 2.9     * root2 / (vev * cosBeta);
   Ye(2,2) = 1.77699 * root2 / (vev * cosBeta);
   mm0 = Sqr(m0) * Eigen::Matrix<double,3,3>::Identity();

   model.set_scale(scale);
   model.set_loops(1);
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

   const int status = minimize(&model);

   BOOST_CHECK_EQUAL(status, GSL_SUCCESS);
}
