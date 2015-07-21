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

#include "gm2_1loop.hpp"
#include "MSSMNoFV_onshell.hpp"
#include "ffunctions.hpp"
#include "linalg2.hpp"
#include "logger.hpp"
#include "wrappers.hpp"

#include <complex>

namespace flexiblesusy {
namespace gm2 {

// approximations

double amuWHnu(const MSSMNoFV_onshell& model) {
   double tan_beta = model.get_TB();
   double M2 = model.get_MassWB();
   double MUE = model.get_Mu();
   double MSv_2 = model.get_MSvmL();

   return ( sqr(model.get_g2()) * 2. * oneOver16PiSqr
            * (sqr(model.get_MM()) * M2 * MUE * tan_beta)
            / sqr(sqr(MSv_2))
            * Fa(sqr(M2 / MSv_2), sqr(MUE / MSv_2)) );
}

double amuWHmuL(const MSSMNoFV_onshell& model) {
   double tan_beta = model.get_TB();
   double M2 = model.get_MassWB();
   double MUE = model.get_Mu();
   double MSL_2 = sqrt(model.get_ml2()(1, 1));

   return ( - sqr(model.get_g2()) * oneOver16PiSqr
            * (sqr(model.get_MM()) * M2 * MUE * tan_beta)
            / sqr(sqr(MSL_2))
            * Fb(sqr(M2 / MSL_2), sqr(MUE / MSL_2)) );
}

double amuBHmuL(const MSSMNoFV_onshell& model) {
   double tan_beta = model.get_TB();
   double M1 = model.get_MassB();
   double MUE = model.get_Mu();
   double MSL_2 = sqrt(model.get_ml2()(1, 1));
   double gY = model.get_gY();

   return ( sqr(gY) * oneOver16PiSqr
            * (sqr(model.get_MM()) * M1 * MUE * tan_beta)
            / sqr(sqr(MSL_2))
            * Fb(sqr(M1 / MSL_2), sqr(MUE / MSL_2)) );
}

double amuBHmuR(const MSSMNoFV_onshell& model) {
   double tan_beta = model.get_TB();
   double M1 = model.get_MassB();
   double MUE = model.get_Mu();
   double MSE_2 = sqrt(model.get_me2()(1, 1));
   double gY = model.get_gY();

   return ( - sqr(gY) * 2. * oneOver16PiSqr
            * (sqr(model.get_MM()) * M1 * MUE * tan_beta)
            / sqr(sqr(MSE_2))
            * Fb(sqr(M1 / MSE_2), sqr(MUE / MSE_2)) );
}

double amuBmuLmuR(const MSSMNoFV_onshell& model) {
   double tan_beta = model.get_TB();
   double M1 = model.get_MassB();
   double MUE = model.get_Mu();
   double MSL_2 = sqrt(model.get_ml2()(1, 1));
   double MSE_2 = sqrt(model.get_me2()(1, 1));
   double gY = model.get_gY();

   return ( sqr(gY) * 2. * oneOver16PiSqr
            * (sqr(model.get_MM()) * MUE * tan_beta)
            / (M1 * sqr(M1))
            * Fb(sqr(MSL_2 / M1), sqr(MSE_2 / M1)) );
}

double amu1Lapprox(const MSSMNoFV_onshell& model) {

   return ( amuWHnu(model) + amuWHmuL(model) + amuBHmuL(model)
            + amuBHmuR(model) + amuBmuLmuR(model) );
}

// complete computation

Eigen::Matrix<std::complex<double>,4,2> n_L(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<std::complex<double>,4,2> result;
   double gY(model.get_gY());
   double g2(model.get_g2());
   double ymu(model.get_Ye()(1, 1));
   Eigen::Matrix<std::complex<double>,4,4> ZN(model.get_ZN());
   Eigen::Array<double,2,1> m_smu(model.get_MSmu());
   Eigen::Matrix<double,2,2> u_smu(model.get_USmu());
   for(int i=0; i <4; ++i) {
      for(int m=0; m <2; ++m) {
         result(i, m) = ( (1. / sqrt(2.) * ( gY * ZN(i, 0)
                           + g2 * ZN(i, 1) ) * Conj(u_smu(m, 0)))
                         - ymu * ZN(i, 2) * Conj(u_smu(m, 1)) );
      }
   }

   return result;
}

Eigen::Matrix<std::complex<double>,4,2> n_R(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<std::complex<double>,4,2> result;
   double gY(model.get_gY());
   double ymu(model.get_Ye()(1, 1));
   Eigen::Matrix<std::complex<double>,4,4> ZN(model.get_ZN());
   Eigen::Array<double,2,1> m_smu(model.get_MSmu());
   Eigen::Matrix<double,2,2> u_smu(model.get_USmu());
   for(int i=0; i <4; ++i) {
      for(int m=0; m <2; ++m) {
         result(i, m) = ( sqrt(2.) * gY * ZN(i, 0) * u_smu(m, 1) 
                       + ymu * ZN(i, 2) * u_smu(m, 0) );
      }
   }

   return result;
}

Eigen::Array<std::complex<double>,2,1> c_L(const MSSMNoFV_onshell& model) {
   Eigen::Array<std::complex<double>,2,1> result;
   double g2 = model.get_g2();
   Eigen::Matrix<std::complex<double>,2,2> UP(model.get_UP());
   for(int k=0; k<2; ++k) {
      result(k) = - g2 * UP(k, 0);
   }

   return result;
}

Eigen::Array<std::complex<double>,2,1> c_R(const MSSMNoFV_onshell& model) {
   Eigen::Array<std::complex<double>,2,1> result;
   double ymu = model.get_Ye()(1, 1);
   Eigen::Matrix<std::complex<double>,2,2> UM(model.get_UM());
   for(int k=0; k<2; ++k) {
      result(k) = ymu * UM(k, 1);
   }

   return result;
}

Eigen::Array<double,2,1> AAC(const MSSMNoFV_onshell& model) {
   Eigen::Array<double,2,1> result;
   Eigen::Array<std::complex<double>,2,1> c__L(c_L(model));
   Eigen::Array<std::complex<double>,2,1> c__R(c_R(model));
   for(int k=0; k<2; ++k) {
      result(k) = norm(c__L(k)) + norm(c__R(k));
   }

   return result;
}

Eigen::Matrix<double,4,2> AAN(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<double,4,2> result;
   Eigen::Matrix<std::complex<double>,4,2> n__L(n_L(model));
   Eigen::Matrix<std::complex<double>,4,2> n__R(n_R(model));
   for(int i=0; i<4; ++i) {
      for(int m=0; m<2; ++m) {
         result(i, m) = norm(n__L(i, m)) + norm(n__R(i, m));
      }
   }

   return result;
}

Eigen::Array<double,2,1> BBC(const MSSMNoFV_onshell& model) {
   Eigen::Array<double,2,1> result;
   Eigen::Array<std::complex<double>,2,1> c__L(c_L(model));
   Eigen::Array<std::complex<double>,2,1> c__R(c_R(model));
   for(int k=0; k<2; ++k) {
      result(k) = real(c__L(k) * c__R(k)) / model.get_MM();
   }

   return result;
}

Eigen::Matrix<double,4,2> BBN(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<double,4,2> result;
   Eigen::Matrix<std::complex<double>,4,2> n__L = n_L(model);
   Eigen::Matrix<std::complex<double>,4,2> n__R = n_R(model);
   for(int i=0; i<4; ++i) {
      for(int m=0; m<2; ++m) {
         result(i, m) = real(n__L(i, m) * n__R(i, m)) / model.get_MM();
      }
   }

   return result;
}

Eigen::Matrix<double,4,2> x_im(const MSSMNoFV_onshell& model) {
   Eigen::Matrix<double,4,2> result;
   Eigen::Array<double,2,1> m_smu(model.get_MSmu());
   Eigen::Matrix<double,2,2> u_smu(model.get_USmu());
   Eigen::Array<double,4,1> MChi(model.get_MChi());
   for(int i=0; i <4; ++i) {
      for(int m=0; m <2; ++m) {
         result(i, m) = sqr(MChi(i) / m_smu(m));
      }
   }

   return result;
}

Eigen::Array<double,2,1> x_k(const MSSMNoFV_onshell& model) {
   Eigen::Array<double,2,1> result;
   for(int k=0; k<2; ++k) {
      result(k) = sqr(model.get_MCha()(k) / model.get_MSvmL()); // !!!
   }

   return result;
}

double amuChi0(const MSSMNoFV_onshell& model) {
   double result = 0.;
   Eigen::Array<double,2,1> m_smu(model.get_MSmu());
   Eigen::Matrix<double,2,2> u_smu(model.get_USmu());
   Eigen::Matrix<double,4,2> AAN_(AAN(model));
   Eigen::Matrix<double,4,2> BBN_(BBN(model));
   Eigen::Matrix<double,4,2> x__im(x_im(model));
   Eigen::Array<double,4,1> MChi(model.get_MChi());
   for(int i=0; i<4; ++i) {
      for(int m=0; m<2; ++m) {
         result += ( - AAN_(i, m) * F1N(x__im(i, m)) / (12. * sqr(m_smu(m)))
                     + MChi(i) * BBN_(i, m) * F2N(x__im(i, m))
                      / (3. * sqr(m_smu(m))) );
      }
   }

   return result * sqr(model.get_MM()) * oneOver16PiSqr;
}

double amuChipm(const MSSMNoFV_onshell& model) {
   double result = 0.;
   Eigen::Array<double,2,1> x__k(x_k(model));
   const double MSvm(model.get_MSvmL());
   Eigen::Array<double,2,1> AAC_(AAC(model));
   Eigen::Array<double,2,1> BBC_(BBC(model));
   Eigen::Array<double,2,1> MCha(model.get_MCha());
   for(int k=0; k<2; ++k) {
      result += ( AAC_(k) * F1C(x__k(k)) / (12. * sqr(MSvm))
                 + 2. * MCha(k) * BBC_(k) * F2C(x__k(k))
                  / (3. * sqr(MSvm)) );
   }

   return result * sqr(model.get_MM()) * oneOver16PiSqr;
}

double calculate_gm2_1loop(const MSSMNoFV_onshell& model) {

   return amuChi0(model) + amuChipm(model);
}

} // namespace gm2
} // namespace flexiblesusy
