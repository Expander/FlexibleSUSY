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

#ifndef PROBLEMS_H
#define PROBLEMS_H

#include <cassert>

namespace flexiblesusy {

template <unsigned Number_of_particles>
class Problems {
public:
   Problems();
   ~Problems() {}

   void flag_tachyon(unsigned);
   void flag_no_ewsb()         { failed_ewsb = true; }
   void flag_no_convergence()  { failed_convergence = true; }
   void flag_no_perturbative() { non_perturbative = true; }

   void unflag_tachyon(unsigned);
   void unflag_no_ewsb()         { failed_ewsb = false; }
   void unflag_no_convergence()  { failed_convergence = false; }
   void unflag_no_perturbative() { non_perturbative = false; }

   bool is_tachyon(unsigned) const;
   bool have_tachyon() const;
   bool no_ewsb() const         { return failed_ewsb; }
   bool no_convergence() const  { return failed_convergence; }
   bool no_perturbative() const { return non_perturbative; }

private:
   bool tachyons[Number_of_particles];
   bool failed_ewsb, failed_convergence, non_perturbative;
};

template <unsigned Number_of_particles>
Problems<Number_of_particles>::Problems()
   : tachyons() // intializes all elements to zero (= false)
   , failed_ewsb(false)
   , failed_convergence(false)
   , non_perturbative(false)
{
}

template <unsigned Number_of_particles>
void Problems<Number_of_particles>::flag_tachyon(unsigned particle)
{
   assert(particle < Number_of_particles
          && "Error: particle index out of bounds");
   tachyons[particle] = true;
}

template <unsigned Number_of_particles>
void Problems<Number_of_particles>::unflag_tachyon(unsigned particle)
{
   assert(particle < Number_of_particles
          && "Error: particle index out of bounds");
   tachyons[particle] = false;
}

template <unsigned Number_of_particles>
bool Problems<Number_of_particles>::is_tachyon(unsigned particle) const
{
   assert(particle < Number_of_particles
          && "Error: particle index out of bounds");
   return tachyons[particle];
}

template <unsigned Number_of_particles>
bool Problems<Number_of_particles>::have_tachyon() const
{
   for (unsigned i = 0; i < Number_of_particles; ++i) {
      if (tachyons[i])
         return true;
   }
   return false;
}

} // namespace flexiblesusy

#endif
