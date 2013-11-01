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

#ifndef pv_hpp
#define pv_hpp

namespace flexiblesusy {

namespace passarino_veltman {

double ReA0(double m2, double scl2);
double ReB0(double p2, double m21, double m22, double scl2);
double ReB1(double p2, double m21, double m22, double scl2);
double ReB00(double p2, double m21, double m22, double scl2);
double ReB22(double p2, double m21, double m22, double scl2);
double ReH0(double p2, double m21, double m22, double scl2);
double ReF0(double p2, double m21, double m22, double scl2);
double ReG0(double p2, double m21, double m22, double scl2);

}

}

#endif // pv_hpp
