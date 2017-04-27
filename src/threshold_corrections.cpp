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

#include "threshold_corrections.hpp"
#include "error.hpp"
#include <cmath>
#include <string>

namespace flexiblesusy {

namespace {

/// returns digit [0-9] in flags at position pos
int get_digit(Threshold_corrections::Flags_t flags, int pos)
{
   if (pos < 0) {
      throw OutOfBoundsError(
         "get_digit: position ( " + std::to_string(pos) + ") must be positive");
   }

   return static_cast<Threshold_corrections::Flags_t>(flags * std::pow(10.l, -pos)) % 10;
}

/// sets digit [0-9] in flags at position pos
void set_digit(Threshold_corrections::Flags_t& flags, int pos, int digit)
{
   if (pos < 0) {
      throw OutOfBoundsError(
         "set_digit: position ( " + std::to_string(pos) + ") must be positive");
   }

   if (digit < 0 || digit > 9) {
      throw OutOfBoundsError(
         "set_digit: digit ( " + std::to_string(digit) + ") must be within [0-9]");
   }

   const int old_digit = get_digit(flags, pos);

   flags += (digit - old_digit) * std::pow(10.l,pos);
}

} // anonymous namespace

Threshold_corrections::Threshold_corrections(Flags_t flags)
{
   set(flags);
}

void Threshold_corrections::set(Flags_t flags)
{
   alpha_em    = get_digit(flags, 0);
   sin_theta_w = get_digit(flags, 1);
   alpha_s     = get_digit(flags, 2);
   mw          = get_digit(flags, 3);
   mz          = get_digit(flags, 4);
   mt          = get_digit(flags, 5);
   mb          = get_digit(flags, 6);
   me          = get_digit(flags, 7);
   mm          = get_digit(flags, 8);
   mtau        = get_digit(flags, 9);
}

Threshold_corrections::Flags_t Threshold_corrections::get() const
{
   Flags_t flags = 0;

   set_digit(flags, 0, alpha_em);
   set_digit(flags, 1, sin_theta_w);
   set_digit(flags, 2, alpha_s);
   set_digit(flags, 3, mw);
   set_digit(flags, 4, mz);
   set_digit(flags, 5, mt);
   set_digit(flags, 6, mb);
   set_digit(flags, 7, me);
   set_digit(flags, 8, mm);
   set_digit(flags, 9, mtau);

   return flags;
}

} // namespace flexiblesusy
