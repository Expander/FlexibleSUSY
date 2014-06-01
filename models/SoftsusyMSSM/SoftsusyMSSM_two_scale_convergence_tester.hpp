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

#ifndef SoftsusyMSSM_TWO_SCALE_CONVERGENCE_TESTER_H
#define SoftsusyMSSM_TWO_SCALE_CONVERGENCE_TESTER_H

#include "two_scale_convergence_tester_drbar.hpp"
#include "SoftsusyMSSM_two_scale.hpp"

namespace flexiblesusy {

class SoftsusyMSSM_convergence_tester : public Convergence_tester_DRbar<SoftsusyMSSM<Two_scale> > {
public:
   SoftsusyMSSM_convergence_tester(SoftsusyMSSM<Two_scale>*, double);
   virtual ~SoftsusyMSSM_convergence_tester();

protected:
   virtual double max_rel_diff() const;

private:
   double sumTol(const SoftsusyMSSM<Two_scale>&, const SoftsusyMSSM<Two_scale>&) const;
};

}

#endif
