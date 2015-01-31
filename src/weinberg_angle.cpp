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

#include "weinberg_angle.hpp"
#include "ew_input.hpp"

namespace flexiblesusy {

namespace weinberg_angle {

Weinberg_angle::Weinberg_angle()
   : number_of_iterations(20)
   , precision_goal(1.0e-8)
   , alpha_em_drbar(Electroweak_constants::aem)
   , fermi_contant(Electroweak_constants::gfermi)
   , self_energy_z_at_mz(0.)
   , self_energy_z_at_0(0.)
   , self_energy_w_at_mw(0.)
{
}

Weinberg_angle::~Weinberg_angle()
{
}

void Weinberg_angle::set_number_of_iterations(unsigned n)
{
   number_of_iterations = n;
}

void Weinberg_angle::set_precision_goal(double p)
{
   precision_goal = p;
}

void Weinberg_angle::set_alpha_em_drbar(double a)
{
   alpha_em_drbar = a;
}

void Weinberg_angle::set_fermi_contant(double gfermi)
{
   fermi_contant = gfermi;
}

void Weinberg_angle::set_self_energy_z_at_mz(double s)
{
   self_energy_z_at_mz = s;
}

void Weinberg_angle::set_self_energy_z_at_0(double s)
{
   self_energy_z_at_0 = s;
}

void Weinberg_angle::set_self_energy_w_at_mw(double s)
{
   self_energy_w_at_mw = s;
}

double Weinberg_angle::calculate() const
{
   return 0.;
}

} // namespace weinberg_angle

} // namespace flexiblesusy
