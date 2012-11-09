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
#include <vector>

class RGE;
class Two_scale;

template<>
class RGFlow<Two_scale> {
public:
   RGFlow(const std::vector<RGE*>&);
   void solve();

private:
   std::vector<RGE*> rge;
};

typedef RGFlow<Two_scale> Two_scale_solver;

inline RGFlow<Two_scale>::RGFlow(const std::vector<RGE*>& rge_)
   : rge(rge_)
{
}

inline void RGFlow<Two_scale>::solve()
{
}

#endif
