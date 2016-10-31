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

#ifndef IF_H
#define IF_H

#include <cstddef>

namespace flexiblesusy {

#define IF(cond,a,b) lazy_if([&](){ return (cond) ? (a) : (b); })

template<class Function>
auto lazy_if(Function f) -> decltype(f())
{
   return f();
}

} // namespace flexiblesusy

#endif // sum_hpp
