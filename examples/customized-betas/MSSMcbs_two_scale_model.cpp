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

#include "MSSMcbs_two_scale_model.hpp"
#include "wrappers.hpp"
#include <Eigen/Core>

namespace flexiblesusy {

using namespace CMSSM_info;

MSSMcbs<Two_scale>::MSSMcbs(const CMSSM_input_parameters& input_)
   : CMSSM<Two_scale>(input_)
{
}

MSSMcbs<Two_scale>::~MSSMcbs()
{
}

Eigen::ArrayXd MSSMcbs<Two_scale>::beta() const
{
   return calc_beta().get();
}

CMSSM_soft_parameters MSSMcbs<Two_scale>::calc_beta() const
{
   CMSSM_soft_parameters betas(CMSSM<Two_scale>::calc_beta());
   if (get_loops() <= 2) return betas;

   // 3-loop g3 beta function from http://www.liv.ac.uk/~dij/betas/

   Eigen::Matrix<double,3,3> Yt = Yu.transpose();
   Eigen::Matrix<double,3,3> Yb = Yd.transpose();
   Eigen::Matrix<double,3,3> Yl = Ye.transpose();
   Eigen::Matrix<double,3,3> Ytc = Yt.adjoint();
   Eigen::Matrix<double,3,3> Ybc = Yb.adjoint();
   Eigen::Matrix<double,3,3> Ylc = Yl.adjoint();
   double a1 = Sqr(g1);
   double a2 = Sqr(g2);
   double a3 = Sqr(g3);

   double bg33 =
      + Sqr(a1)*a3 * ( - 1702/75.0 )
      + a1*a2*a3 * ( - 3/5.0 )
      + a1*Sqr(a3) * ( 22/15.0 )
      + a1*a3 * ( - 44/15.0*(Ytc*Yt).trace() - 32/15.0*(Ybc*Yb).trace() )
      + Sqr(a2)*a3 * ( - 27 )
      + a2*Sqr(a3) * ( 6 )
      + a2*a3 * ( - 12*(Ytc*Yt).trace() - 12*(Ybc*Yb).trace() )
      + Power(a3,3) * ( 347/3.0 )
      + Sqr(a3) * ( - 104/3.0*(Ytc*Yt).trace() - 104/3.0*(Ybc*Yb).trace() )
      + a3 * ( 18*Sqr((Ytc*Yt).trace()) + 12*(Ytc*Yt*Ytc*Yt).trace()
	     + 8*(Ytc*Yt*Ybc*Yb).trace()
	     + 18*Sqr((Ybc*Yb).trace()) + 6*(Ybc*Yb).trace()*(Ylc*Yl).trace()
	     + 12*(Ybc*Yb*Ybc*Yb).trace() );

   betas.set_g3(betas.get_g3() + Power(oneOver16PiSqr,3) * g3 * bg33);

   return betas;
}

} // namespace flexiblesusy
