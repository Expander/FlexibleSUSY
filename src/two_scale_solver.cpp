
#include "two_scale_solver.hpp"
#include "linalg.h"
#include "rge.h"

#include <cassert>

Two_scale_solver::Two_scale_solver(RGE* rge_)
   : rge(rge_)
{
   assert(rge_ && "rge must not be 0");
}

Two_scale_solver::~Two_scale_solver()
{
}

DoubleVector Two_scale_solver::solve(double scale)
{
   DoubleVector v(rge->howMany());
   return v;
}
