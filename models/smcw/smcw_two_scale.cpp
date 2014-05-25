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

#include "smcw_two_scale.hpp"

#include <cassert>

namespace flexiblesusy {

StandardModelCW<Two_scale>::StandardModelCW()
   : StandardModel<Two_scale>()
   , g4(0.0)
   , lambda(0.0)
   , vs(0.0)
{
   setPars(numStandardModelCWPars);
}

StandardModelCW<Two_scale>::StandardModelCW(const StandardModelCW<Two_scale>& smcw)
   : StandardModel<Two_scale>(smcw)
   , g4(smcw.g4)
   , lambda(smcw.lambda)
   , vs(smcw.vs)
{
   setPars(numStandardModelCWPars);
}

StandardModelCW<Two_scale>::StandardModelCW(const StandardModel<Two_scale>& sm, double g4_, double lambda_, double vs_)
   : StandardModel<Two_scale>(sm)
   , g4(g4_)
   , lambda(lambda_)
   , vs(vs_)
{
   setPars(numStandardModelCWPars);
}

StandardModelCW<Two_scale>::StandardModelCW(const DoubleMatrix& SMu,
                                            const DoubleMatrix& SMd,
                                            const DoubleMatrix& SMe,
                                            const DoubleVector& g_,
                                            double lambda_, double vs_)
   : StandardModel<Two_scale>(SMu, SMd, SMe, g_)
   , g4(g_(4))
   , lambda(lambda_)
   , vs(vs_)
{
   setPars(numStandardModelCWPars);
}

StandardModelCW<Two_scale>::~StandardModelCW()
{
}

const StandardModelCW<Two_scale>& StandardModelCW<Two_scale>::operator=(const StandardModelCW<Two_scale>& smcw)
{
   if (this == &smcw) return *this;
   StandardModel<Two_scale>::operator=(smcw);
   g4 = smcw.g4;
   lambda = smcw.lambda;
   vs = smcw.vs;
   return *this;
}

void StandardModelCW<Two_scale>::setGaugeCoupling(int i, double f)
{
   assert(i <= 4 && "gauge coupling index > 4");
   assert(i >= 1 && "gauge coupling index < 1");
   if (i == 4)
      g4 = f;
   else
      StandardModel<Two_scale>::setGaugeCoupling(i, f);
}

void StandardModelCW<Two_scale>::setAllGauge(const DoubleVector& v)
{
   assert(v.displayStart() == 1 && v.displayEnd() >= 4);
   g4 = v(4);
   for (int i = 1; i <= 3; ++i)
      StandardModel<Two_scale>::setGaugeCoupling(i, v(i));
}

DoubleVector StandardModelCW<Two_scale>::displayGauge() const
{
   DoubleVector g(StandardModel<Two_scale>::displayGauge());
   g.setEnd(4);
   g(4) = g4;
   return g;
}

double StandardModelCW<Two_scale>::displayGaugeCoupling(int i) const
{
   assert(i >= 1 && "i < 1");
   assert(i <= 4 && "i > 4");
   if (i == 4)
      return g4;
   else
      return StandardModel<Two_scale>::displayGaugeCoupling(i);
}

const DoubleVector StandardModelCW<Two_scale>::display() const
{
   DoubleVector y(StandardModel<Two_scale>::display());
   y.setEnd(numStandardModelCWPars);
   y(numStandardModelCWPars - 2) = g4;
   y(numStandardModelCWPars - 1) = lambda;
   y(numStandardModelCWPars)     = vs;

   return y;
}

void StandardModelCW<Two_scale>::set(const DoubleVector& y)
{
   assert(y.displayStart() == 1 && y.displayEnd() >= numStandardModelCWPars);
   StandardModel<Two_scale>::set(y);
   g4     = y(numStandardModelCWPars - 2);
   lambda = y(numStandardModelCWPars - 1);
   vs     = y(numStandardModelCWPars);
}

std::ostream& operator <<(std::ostream& left, const StandardModelCW<Two_scale>& s)
{
   left << "SMCW parameters at Q: " << s.get_scale()
        << '\n'
        << " Y^U" << s.displayYukawaMatrix(StandardModelCW<Two_scale>::YU)
        << " Y^D" << s.displayYukawaMatrix(StandardModelCW<Two_scale>::YD)
        << " Y^E" << s.displayYukawaMatrix(StandardModelCW<Two_scale>::YE)
        << '\n'
        << " g1: " << s.displayGaugeCoupling(1)
        << " g2: " << s.displayGaugeCoupling(2)
        << " g3: " << s.displayGaugeCoupling(3)
        << " g4: " << s.displayGaugeCoupling(4)
        << '\n'
        << " lambda: " << s.displayLambda()
        << " vs: " << s.displayVs()
        << '\n'
        << " thresholds: " << s.displayThresholds()
        << " #loops: " << s.displayLoops() << '\n';
   return left;
}

// Outputs derivatives (DRbar scheme) in the form of ds
StandardModelCW<Two_scale> StandardModelCW<Two_scale>::calc_beta() const
{
   static const double oneO16Pisq = 1.0 / (16.0 * PI * PI);
   double dg4, dlambda;

   dg4 = oneO16Pisq * std::pow(g4, 3) / 3.0;
   dlambda = oneO16Pisq * (2.0 / 3.0) * (5.0 * lambda * lambda
                                         - 18.0 * g4 * g4 * lambda
                                         + 54.0 * std::pow(g4, 4));

   return StandardModelCW(StandardModel<Two_scale>::calc_beta(), dg4, dlambda, 0.0);
}

DoubleVector StandardModelCW<Two_scale>::beta() const
{
   return calc_beta().display();
}

double StandardModelCW<Two_scale>::calcZprimeMass() const
{
   return vs * g4;
}

} // namespace flexiblesusy
