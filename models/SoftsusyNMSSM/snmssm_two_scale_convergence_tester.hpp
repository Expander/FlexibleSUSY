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

#ifndef SoftsusyNMSSM_TWO_SCALE_CONVERGENCE_TESTER_H
#define SoftsusyNMSSM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "two_scale_convergence_tester_skeleton.hpp"
#include "snmssm_two_scale.hpp"

namespace flexiblesusy {

class SNmssm_convergence_tester : public Convergence_tester_skeleton<SNmssm<Two_scale> > {
public:
   SNmssm_convergence_tester(SNmssm<Two_scale>*, double);
   virtual ~SNmssm_convergence_tester();

protected:
   virtual double max_rel_diff() const;

private:
   double sumTol(const SNmssm<Two_scale>&, const SNmssm<Two_scale>&) const;
};

}

#endif
