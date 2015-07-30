// ====================================================================
// This file is part of GM2Calc.
//
// GM2Calc is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// GM2Calc is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with GM2Calc.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef GM2_FFUNCTIONS_H
#define GM2_FFUNCTIONS_H

namespace flexiblesusy {
namespace gm2calc {

double F1C(double);
double F2C(double);
double F3C(double);
double F4C(double);
double F1N(double);
double F2N(double);
double F3N(double);
double F4N(double);
double Fa(double, double);
double Fb(double, double);
double G3(double);
double G4(double);
double H2(double, double);
double Iabc(double, double, double);
/// \f$f_{PS}(z)\f$, Eq. (70) arXiv:hep-ph/0609168
double f_PS(double);
/// \f$f_S(z)\f$, Eq. (71) arXiv:hep-ph/0609168
double f_S(double);
/// \f$f_{\tilde{f}}(z)\f$, Eq. (72) arXiv:hep-ph/0609168
double f_sferm(double);

double cube(double);
double quad(double);
int sign(double);
double signed_abs_sqrt(double);
double signed_sqr(double);
double sqr(double);

} // namespace gm2calc
} // namespace flexiblesusy

#endif
