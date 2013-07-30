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

#include "stopwatch.hpp"

namespace flexiblesusy {

Stopwatch::Stopwatch()
   : time(0)
{
}

Stopwatch::~Stopwatch()
{
}

void Stopwatch::start()
{
   time = clock();
}

void Stopwatch::stop()
{
   time = clock() - time;
}

double Stopwatch::get_clicks()
{
   return time;
}

double Stopwatch::get_time_in_seconds()
{
   return static_cast<double>(time)/CLOCKS_PER_SEC;
}

} // namespace flexiblesusy
