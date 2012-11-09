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

#ifndef TWO_SCALE_SOLVER_H
#define TWO_SCALE_SOLVER_H

#include "rg_flow.hpp"
#include "two_scale_model.hpp"
#include <vector>

class Two_scale;
class Two_scale_model;
class Two_scale_matching;

template<>
class RGFlow<Two_scale> {
public:
   RGFlow(const std::vector<Two_scale_model*>&);

   void addMatchingCondition(const Two_scale_matching*);
   void run_up();
   void run_down();
   void solve();

private:
   std::vector<Two_scale_model*> rge;
   std::vector<const Two_scale_matching*> matching;
   unsigned int maxIterations;

   bool accuracyGoalReached() const;
};

typedef RGFlow<Two_scale> Two_scale_solver;

inline RGFlow<Two_scale>::RGFlow(const std::vector<Two_scale_model*>& rge_)
   : rge(rge_)
   , matching()
   , maxIterations(10)
{
}

inline void RGFlow<Two_scale>::solve()
{
   // initial run
   run_up();
   run_down();

   for (unsigned int iter = 0;
        iter < maxIterations && !accuracyGoalReached();
        ++iter) {
      run_up();
      run_down();
   }
}

inline void RGFlow<Two_scale>::run_up()
{
   for (std::vector<Two_scale_model*>::iterator model = rge.begin(), end = rge.end();
        model != end; ++model) {
      (*model)->run_up();
   }
}

inline void RGFlow<Two_scale>::run_down()
{
   for (std::vector<Two_scale_model*>::iterator model = rge.begin(), end = rge.end();
        model != end; ++model) {
      (*model)->run_down();
   }
}

inline void RGFlow<Two_scale>::addMatchingCondition(const Two_scale_matching* mc)
{
   matching.push_back(mc);
}

inline bool RGFlow<Two_scale>::accuracyGoalReached() const
{
   return true;
}

#endif
