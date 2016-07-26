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


#include "standard_model_two_scale_low_scale_constraint.hpp"
#include "standard_model_two_scale_model.hpp"
#include "wrappers.hpp"
#include "logger.hpp"
#include "ew_input.hpp"
#include "gsl_utils.hpp"
#include "minimizer.hpp"
#include "root_finder.hpp"
#include "threshold_loop_functions.hpp"
#include "weinberg_angle.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace flexiblesusy {
namespace standard_model {

#define MODELPARAMETER(p) model->get_##p()
#define PHASE(p) model->get_##p()
#define BETAPARAMETER(p) beta_functions.get_##p()
#define BETA(p) beta_##p
#define LowEnergyConstant(p) Electroweak_constants::p
#define MZPole qedqcd.displayPoleMZ()
#define STANDARDDEVIATION(p) Electroweak_constants::Error_##p
#define Pole(p) model->get_physical().p
#define SCALE model->get_scale()
#define MODEL model
#define MODELCLASSNAME StandardModel<Two_scale>
#define CKM ckm
#define PMNS pmns
#define THETAW theta_w
#define THRESHOLD static_cast<int>(model->get_thresholds())
#define ALPHA_EM_DRBAR alpha_em_drbar
#define CALCULATE_DRBAR_MASSES() model->calculate_DRbar_masses()

Standard_model_low_scale_constraint<Two_scale>::Standard_model_low_scale_constraint()
   : Constraint<Two_scale>()
   , scale(0.)
   , model(0)
   , qedqcd()
   , ckm()
   , pmns()
   , threshold_corrections_loop_order(1)
{
   ckm << 1., 0., 0.,
          0., 1., 0.,
          0., 0., 1.;

   pmns << 1., 0., 0.,
           0., 1., 0.,
           0., 0., 1.;
}

Standard_model_low_scale_constraint<Two_scale>::Standard_model_low_scale_constraint(
   StandardModel<Two_scale>* model_, const softsusy::QedQcd& qedqcd_)
   : Constraint<Two_scale>()
   , model(model_)
   , qedqcd(qedqcd_)
{
   initialize();
}

Standard_model_low_scale_constraint<Two_scale>::~Standard_model_low_scale_constraint()
{
}

void Standard_model_low_scale_constraint<Two_scale>::apply()
{
   assert(model && "Error: Standard_model_low_scale_constraint::apply():"
          " model pointer must not be zero");

   qedqcd.runto(scale, 1.0e-5);
   model->set_low_energy_data(qedqcd);
   model->calculate_DRbar_masses();

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mz_pole  = qedqcd.displayPoleMZ();

   double delta_alpha_em = 0.;
   double delta_alpha_s  = 0.;

   if (model->get_thresholds()) {
      delta_alpha_em = model->calculate_delta_alpha_em(alpha_em);
      delta_alpha_s  = model->calculate_delta_alpha_s(alpha_s);
   }

   const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
   const double alpha_s_drbar  = alpha_s  / (1.0 - delta_alpha_s);
   const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);
   const double mZ = model->get_thresholds() ?
      model->calculate_MVZ_DRbar(mz_pole) : mz_pole;
   const double theta_w = model->calculate_theta_w(alpha_em_drbar);

   const auto g1 = MODELPARAMETER(g1);
   const auto g2 = MODELPARAMETER(g2);

   MODEL->set_v(Re((2*mZ)/Sqrt(0.6*Sqr(g1) + Sqr(g2))));
   model->calculate_Yu_DRbar();
   model->calculate_Yd_DRbar();
   model->calculate_Ye_DRbar();

   model->set_g1(1.2909944487358056*e_drbar*Sec(theta_w));
   model->set_g2(e_drbar*Csc(theta_w));
   model->set_g3(3.5449077018110318*Sqrt(alpha_s_drbar));

   model->recalculate_mw_pole();
}

const Eigen::Matrix<std::complex<double>,3,3>& Standard_model_low_scale_constraint<Two_scale>::get_ckm()
{
   return ckm;
}

const Eigen::Matrix<std::complex<double>,3,3>& Standard_model_low_scale_constraint<Two_scale>::get_pmns()
{
   return pmns;
}

double Standard_model_low_scale_constraint<Two_scale>::get_scale() const
{
   return scale;
}

void Standard_model_low_scale_constraint<Two_scale>::set_model(Two_scale_model* model_)
{
   model = cast_model<StandardModel<Two_scale>*>(model_);
}

void Standard_model_low_scale_constraint<Two_scale>::set_sm_parameters(
   const softsusy::QedQcd& qedqcd_)
{
   qedqcd = qedqcd_;
}

const softsusy::QedQcd& Standard_model_low_scale_constraint<Two_scale>::get_sm_parameters() const
{
   return qedqcd;
}

void Standard_model_low_scale_constraint<Two_scale>::clear()
{
   scale = 0.;
   model = NULL;
   qedqcd = softsusy::QedQcd();
}

void Standard_model_low_scale_constraint<Two_scale>::initialize()
{
   assert(model && "Standard_model_low_scale_constraint<Two_scale>::"
          "initialize(): model pointer is zero.");

   scale = MZPole;
   ckm = qedqcd.get_complex_ckm();
   pmns = qedqcd.get_complex_pmns();
}

} // namespace standard_model
} // namespace flexiblesusy
