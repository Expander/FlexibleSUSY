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

#include "MSSM_MSSMRHN_two_scale_initial_guesser.hpp"
#include "MSSM_two_scale_model.hpp"
#include "MSSMRHN_two_scale_model.hpp"
#include "lowe.h"
#include "error.hpp"
#include "ew_input.hpp"
#include "wrappers.hpp"

#include <Eigen/Core>
#include <cassert>

namespace flexiblesusy {

#define INPUTPARAMETER(p) input_pars.p
#define MODEL1PARAMETER(p) model_1->get_##p()
#define MODEL2PARAMETER(p) model_2->get_##p()
#define SM(p) Electroweak_constants::p
// #define MODEL model

MSSM_MSSMRHN_initial_guesser<Two_scale>::MSSM_MSSMRHN_initial_guesser(
   MSSM<Two_scale>* model_1_, MSSMRHN<Two_scale>* model_2_,
   const MSSMRHN_input_parameters& input_pars_,
   const QedQcd& oneset_,
   const MSSM_low_scale_constraint<Two_scale>& low_constraint_,
   const MSSM_susy_scale_constraint<Two_scale>& susy_constraint_,
   const MSSMRHN_high_scale_constraint<Two_scale>& high_constraint_,
   const MSSM_MSSMRHN_matching<Two_scale>& matching_
)
   : Initial_guesser<Two_scale>()
   , model_1(model_1_), model_2(model_2_)
   , input_pars(input_pars_)
   , oneset(oneset_)
   , low_constraint(low_constraint_)
   , susy_constraint(susy_constraint_)
   , high_constraint(high_constraint_)
   , matching(matching_)
{
   assert(model_1 && model_2 && "MSSM_MSSMRHN_initial_guesser: Error: pointers to models must not be zero");
}

MSSM_MSSMRHN_initial_guesser<Two_scale>::~MSSM_MSSMRHN_initial_guesser()
{
}

void MSSM_MSSMRHN_initial_guesser<Two_scale>::guess()
{
   guess_susy_parameters();
   guess_soft_parameters();
}

void MSSM_MSSMRHN_initial_guesser<Two_scale>::guess_susy_parameters()
{
   QedQcd leAtMt(oneset);
   const double MZ = Electroweak_constants::MZ;
   const double MW = Electroweak_constants::MW;
   const double sinThetaW2 = 1.0 - Sqr(MW / MZ);
   const double mtpole = leAtMt.displayPoleMt();

   // guess gauge couplings at mt
   const DoubleVector alpha_sm(leAtMt.getGaugeMu(mtpole, sinThetaW2));

   model_1->set_g1(sqrt(4.0 * M_PI * alpha_sm(1)));
   model_1->set_g2(sqrt(4.0 * M_PI * alpha_sm(2)));
   model_1->set_g3(sqrt(4.0 * M_PI * alpha_sm(3)));
   model_1->set_scale(mtpole);
   model_1->set_loops(2);

   // apply user-defined initial guess at the low scale
   const auto TanBeta = INPUTPARAMETER(TanBeta);

   model_1->set_vd(SM(vev)/Sqrt(1 + Sqr(TanBeta)));
   model_1->set_vu((TanBeta*SM(vev))/Sqrt(1 + Sqr(TanBeta)));

   Eigen::Matrix<double,3,3> topDRbar(Eigen::Matrix<double,3,3>::Zero()),
      bottomDRbar(Eigen::Matrix<double,3,3>::Zero()),
      electronDRbar(Eigen::Matrix<double,3,3>::Zero());
   topDRbar(0,0)      = leAtMt.displayMass(mUp);
   topDRbar(1,1)      = leAtMt.displayMass(mCharm);
   topDRbar(2,2)      = leAtMt.displayMass(mTop) - 30.0;
   bottomDRbar(0,0)   = leAtMt.displayMass(mDown);
   bottomDRbar(1,1)   = leAtMt.displayMass(mStrange);
   bottomDRbar(2,2)   = leAtMt.displayMass(mBottom);
   electronDRbar(0,0) = leAtMt.displayMass(mElectron);
   electronDRbar(1,1) = leAtMt.displayMass(mMuon);
   electronDRbar(2,2) = leAtMt.displayMass(mTau);

   Eigen::Matrix<double,3,3> new_Yu, new_Yd, new_Ye;

   const auto vd = MODEL1PARAMETER(vd);
   const auto vu = MODEL1PARAMETER(vu);
   new_Yu = ((1.4142135623730951*topDRbar)/vu).transpose();
   new_Yd = ((1.4142135623730951*bottomDRbar)/vd).transpose();
   new_Ye = ((1.4142135623730951*electronDRbar)/vd).transpose();

   model_1->set_Yu(new_Yu);
   model_1->set_Yd(new_Yd);
   model_1->set_Ye(new_Ye);
}

void MSSM_MSSMRHN_initial_guesser<Two_scale>::guess_soft_parameters()
{
   const double low_scale_guess = low_constraint.get_initial_scale_guess();
   const double high_scale_guess = high_constraint.get_initial_scale_guess();
   const double matching_scale_guess = matching.get_initial_scale_guess();

   model_1->run_to(matching_scale_guess);
   matching.set_models(model_1, model_2);
   matching.match_low_to_high_scale_model();

   model_2->run_to(high_scale_guess);

   // apply high-scale constraint
   high_constraint.set_model(model_2);
   high_constraint.apply();

   // apply user-defined initial guess at the high scale
   model_2->set_Mu(1.);
   model_2->set_BMu(0.);

   model_2->run_to(matching_scale_guess);
   matching.match_high_to_low_scale_model();

   model_1->run_to(low_scale_guess);

   // apply EWSB constraint
   model_1->solve_ewsb_tree_level();

   // calculate tree-level spectrum
   model_1->calculate_DRbar_parameters();
   model_1->set_thresholds(3);
   model_1->set_loops(2);
}

} // namespace flexiblesusy
