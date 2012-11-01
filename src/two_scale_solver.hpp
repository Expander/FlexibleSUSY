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

#include "rge_flow.hpp"
#include <cassert>

class RGE;
class Two_scale;

template<>
class RGEFlow<Two_scale> {
public:
   RGEFlow(RGE*);
   void solve();

private:
   RGE* rge;
};

typedef RGEFlow<Two_scale> Two_scale_solver;

RGEFlow<Two_scale>::RGEFlow(RGE* rge_)
   : rge(rge_)
{
   assert(rge_ && "rge must not be 0");
}

void RGEFlow<Two_scale>::solve()
{
}

#endif
