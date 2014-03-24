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

#include <cmath>
#include <Eigen/LU>
#include "numerics.hpp"

#include "two_scale_constraint.hpp" // for cast_model
#include "MSSM_MSSMRHN_two_scale_matching.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSMRHN_two_scale_model.hpp"

using namespace std;
using namespace flexiblesusy;

MSSM_MSSMRHN_matching<Two_scale>::MSSM_MSSMRHN_matching() :
    Matching<Two_scale>(),
    fixed_scale(0),
    lower(0), upper(0),
    inputPars()
{
    // no reasonable guess without input
    scale = initial_scale_guess = 0;
}

MSSM_MSSMRHN_matching<Two_scale>::MSSM_MSSMRHN_matching(const MSSMRHN_input_parameters& inputPars_) :
    Matching<Two_scale>(),
    lower(0), upper(0),
    inputPars(inputPars_)
{
    make_initial_scale_guess();
}

void MSSM_MSSMRHN_matching<Two_scale>::match_low_to_high_scale_model()
{
    upper->set_Yd(lower->get_Yd());
    upper->set_Ye(lower->get_Ye());
    upper->set_Yu(lower->get_Yu());
    upper->set_Mu(lower->get_Mu());
    upper->set_g1(lower->get_g1());
    upper->set_g2(lower->get_g2());
    upper->set_g3(lower->get_g3());
    upper->set_vd(lower->get_vd());
    upper->set_vu(lower->get_vu());

    upper->set_TYd(lower->get_TYd());
    upper->set_TYe(lower->get_TYe());
    upper->set_TYu(lower->get_TYu());
    upper->set_BMu(lower->get_BMu());
    upper->set_mq2(lower->get_mq2());
    upper->set_ml2(lower->get_ml2());
    upper->set_mHd2(lower->get_mHd2());
    upper->set_mHu2(lower->get_mHu2());
    upper->set_md2(lower->get_md2());
    upper->set_mu2(lower->get_mu2());
    upper->set_me2(lower->get_me2());
    upper->set_MassB(lower->get_MassB());
    upper->set_MassWB(lower->get_MassWB());
    upper->set_MassG(lower->get_MassG());

    upper->set_scale(lower->get_scale());
}

void MSSM_MSSMRHN_matching<Two_scale>::match_high_to_low_scale_model()
{
    update_scale();

    lower->set_Yd(upper->get_Yd());
    lower->set_Ye(upper->get_Ye());
    lower->set_Yu(upper->get_Yu());
    lower->set_Mu(upper->get_Mu());
    lower->set_g1(upper->get_g1());
    lower->set_g2(upper->get_g2());
    lower->set_g3(upper->get_g3());
    lower->set_vd(upper->get_vd());
    lower->set_vu(upper->get_vu());

    lower->set_TYd(upper->get_TYd());
    lower->set_TYe(upper->get_TYe());
    lower->set_TYu(upper->get_TYu());
    lower->set_BMu(upper->get_BMu());
    lower->set_mq2(upper->get_mq2());
    lower->set_ml2(upper->get_ml2());
    lower->set_mHd2(upper->get_mHd2());
    lower->set_mHu2(upper->get_mHu2());
    lower->set_md2(upper->get_md2());
    lower->set_mu2(upper->get_mu2());
    lower->set_me2(upper->get_me2());
    lower->set_MassB(upper->get_MassB());
    lower->set_MassWB(upper->get_MassWB());
    lower->set_MassG(upper->get_MassG());

    lower->set_scale(upper->get_scale());
}

double MSSM_MSSMRHN_matching<Two_scale>::get_scale() const
{
    return scale;
}

double MSSM_MSSMRHN_matching<Two_scale>::get_initial_scale_guess() const
{
   return initial_scale_guess;
}

void MSSM_MSSMRHN_matching<Two_scale>::set_models(Two_scale_model *lower_, Two_scale_model *upper_)
{
    lower = cast_model<MSSM<Two_scale> >(lower_);
    upper = cast_model<MSSMRHN<Two_scale> >(upper_);
}

void MSSM_MSSMRHN_matching<Two_scale>::set_upper_input_parameters(const MSSMRHN_input_parameters& inputPars_)
{
   inputPars = inputPars_;
   make_initial_scale_guess();
}

void MSSM_MSSMRHN_matching<Two_scale>::make_initial_scale_guess()
{
    double RHN_scale = pow(abs(inputPars.MvInput.determinant()), 1.0/3);
    scale = initial_scale_guess = RHN_scale;
}

void MSSM_MSSMRHN_matching<Two_scale>::reset()
{
   scale = initial_scale_guess;
   fixed_scale = 0.;
   lower = NULL;
   upper = NULL;
}

void MSSM_MSSMRHN_matching<Two_scale>::update_scale()
{
    if (!is_zero(fixed_scale)) {
	scale = fixed_scale;
	return;
    }

    double RHN_scale = pow(abs(upper->get_Mv().determinant()), 1.0/3);
    scale = RHN_scale;
}
