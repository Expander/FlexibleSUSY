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
#include "logger.hpp"

#include <limits>

Mssm_convergence_tester::Mssm_convergence_tester(Mssm<Two_scale>* mssm_, double accuracy_goal_, unsigned long max_iterations_)
   : Convergence_tester<Two_scale>()
   , mssm(mssm_)
   , last_iteration()
   , it_count(0)
   , maximum_iterations(max_iterations_)
   , accuracy_goal(accuracy_goal_)
{
}

Mssm_convergence_tester::~Mssm_convergence_tester()
{
}

bool Mssm_convergence_tester::accuracy_goal_reached()
{
   bool precision_reached;
   if (it_count == 0) {
      // this is the first run => no comparison possible => assume
      // that accuracy goal has not been reached
      precision_reached = false;
   } else {
      if (scale_has_changed() && rel_scale_difference() > accuracy_goal) {
         WARNING("scale has changed by " << scale_difference()
                 << " GeV (" << rel_scale_difference()
                 << "%), parameter comparison might fail");
      }

      const double fracDiff = sumTol(*mssm, last_iteration);
      precision_reached = fracDiff < accuracy_goal;
      VERBOSE_MSG("Mssm_convergence_tester: max. rel. difference: " << fracDiff
                  << " (goal: " << accuracy_goal << ")");
   }

   // save old model parameters
   last_iteration = *mssm;
   ++it_count;

   return precision_reached;
}

unsigned int Mssm_convergence_tester::max_iterations() const
{
   return maximum_iterations;
}

bool Mssm_convergence_tester::scale_has_changed() const
{
   return !is_zero(scale_difference());
}

double Mssm_convergence_tester::scale_difference() const
{
   return mssm->getScale() - last_iteration.getScale();
}

double Mssm_convergence_tester::rel_scale_difference() const
{
   const double diff = scale_difference();
   if (!is_zero(last_iteration.getScale()))
      return diff / last_iteration.getScale();
   return std::numeric_limits<double>::infinity();
}

double Mssm_convergence_tester::sumTol(const Mssm<Two_scale>& in, const Mssm<Two_scale>& out)
{
  drBarPars inforLoops(in.displayDrBarPars()),
    outforLoops(out.displayDrBarPars());

  DoubleVector sT(32);
  int k = 1;

  double sTin  = fabs(inforLoops.mh0); double sTout = fabs(outforLoops.mh0);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mA0); sTout = fabs(outforLoops.mA0);
  sT(k) = fabs(1.0 - minimum(sTin, sTout) / maximum(sTin, sTout)); k++;
  sTin  = fabs(inforLoops.mH0); sTout = fabs(outforLoops.mH0);
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
