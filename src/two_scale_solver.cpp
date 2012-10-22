
#include "two_scale_solver.hpp"
#include "linalg.h"
#include "rge.h"

Two_scale_solver::Two_scale_solver(RGE* rge_)
   : rge(rge_)
{
}

Two_scale_solver::~Two_scale_solver()
{
}

DoubleVector Two_scale_solver::solve(double scale)
{
   DoubleVector v(10);
   return v;
}
