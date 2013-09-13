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

#include "mssm_two_scale_convergence_tester.hpp"

namespace flexiblesusy {

Mssm_convergence_tester::Mssm_convergence_tester(Mssm<Two_scale>* mssm_, double accuracy_goal_)
   : Convergence_tester_skeleton<Mssm<Two_scale> >(mssm_, accuracy_goal_)
{
}

Mssm_convergence_tester::~Mssm_convergence_tester()
{
}

double Mssm_convergence_tester::max_rel_diff() const
{
   return sumTol(*get_model(), *get_last_iteration_model());
}

double Mssm_convergence_tester::sumTol(const Mssm<Two_scale>& in, const Mssm<Two_scale>& out) const
{
  drBarPars inforLoops(in.displayDrBarPars()),
    outforLoops(out.displayDrBarPars());

  DoubleVector sT(32);
  int k = 1;

  double sTin  = fabs(inforLoops.mh0(1)); double sTout = fabs(outforLoops.mh0(1));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mA0(1)); sTout = fabs(outforLoops.mA0(1));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mh0(2)); sTout = fabs(outforLoops.mh0(2));
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mHpm); sTout = fabs(outforLoops.mHpm);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  int i; for (i=1; i<=3; i++) {
    sTin  = fabs(inforLoops.msnu(i));
    sTout = fabs(outforLoops.msnu(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }
  for (i=1; i<=2; i++) {
    sTin = fabs(inforLoops.mch(i));
    sTout = fabs(outforLoops.mch(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }
  for (i=1; i<=4; i++) {
    sTin = fabs(inforLoops.mneut(i));
    sTout = fabs(outforLoops.mneut(i));
    sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
    k++;
  }
  sTin = fabs(inforLoops.mGluino);
  sTout = fabs(outforLoops.mGluino);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
  k++;
  int j; for (j=1; j<=3; j++)
    for(i=1; i<=2; i++) {
      sTin = fabs(inforLoops.mu(i, j));
      sTout = fabs(outforLoops.mu(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k++;
      sTin = fabs(inforLoops.md(i, j));
      sTout = fabs(outforLoops.md(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k++;
      sTin = fabs(inforLoops.me(i, j));
      sTout = fabs(outforLoops.me(i, j));
      sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout));
      k++;
    }
  return sT.max();
}

}
