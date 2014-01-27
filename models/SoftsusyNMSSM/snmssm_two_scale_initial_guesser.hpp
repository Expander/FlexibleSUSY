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

#ifndef SoftsusyNMSSM_TWO_SCALE_INITIAL_GUESSER_H
#define SoftsusyNMSSM_TWO_SCALE_INITIAL_GUESSER_H

#include "two_scale_initial_guesser.hpp"
#include "snmssm_parameter_point.hpp"
#include "lowe.h"
#include "linalg.h"

namespace flexiblesusy {

template<class T> class SNmssm;
class Two_scale;

class SoftsusyNMSSM_low_scale_constraint;
class SoftsusyNMSSM_susy_scale_constraint;
class SoftsusyNMSSM_sugra_constraint;

class SoftsusyNMSSM_initial_guesser : public Initial_guesser<Two_scale> {
public:
   SoftsusyNMSSM_initial_guesser(SNmssm<Two_scale>*, const SoftsusyNMSSM_parameter_point&,
                          const SoftsusyNMSSM_low_scale_constraint&,
                          const SoftsusyNMSSM_susy_scale_constraint&,
                          const SoftsusyNMSSM_sugra_constraint&);
   virtual ~SoftsusyNMSSM_initial_guesser();
   virtual void guess();
   void set_QedQcd(const QedQcd& qedqcd) { oneset = qedqcd; }

private:
   SNmssm<Two_scale>* nmssm;  ///< SNmssm model
   QedQcd oneset;             ///< low-energy parameters
   SoftsusyNMSSM_parameter_point pp; ///< SNmssm parameter point
   bool ewsbBCscale;          ///< EWSB at susy scale
};

} // namespace flexiblesusy

#endif
