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

#ifndef ARRAY_VIEW_H
#define ARRAY_VIEW_H

#include "error.hpp"
#include <cstddef>
#include <string>

namespace flexiblesusy {

template <typename ElementType>
class Dynamic_array_view
{
public:
   using Element_t = ElementType;
   using Index_t = std::ptrdiff_t;
   using Pointer_t = Element_t*;
   using Iterator_t = Pointer_t;
   using Const_iterator_t = const Iterator_t;
   using Reference_t = Element_t&;

   constexpr Dynamic_array_view() noexcept = default;
   constexpr Dynamic_array_view(Pointer_t p, Index_t l) noexcept
      : ptr(p), len(l) {}
   constexpr Dynamic_array_view(Pointer_t first, Pointer_t last) noexcept
      : ptr(first), len(std::distance(first, last)) {}

   constexpr Pointer_t data() const noexcept { return ptr; }
   constexpr bool empty() const noexcept { return size() == 0; }
   constexpr Index_t size() const noexcept { return len; }

   Iterator_t begin() const noexcept { return ptr; }
   Iterator_t end() const noexcept { return ptr + len; }

   Const_iterator_t cbegin() const noexcept { return ptr; }
   Const_iterator_t cend() const noexcept { return ptr + len; }

   Reference_t operator[](Index_t idx) const {
      check_range(idx);
      return ptr[idx];
   }

private:
   Pointer_t ptr{nullptr};
   Index_t len{0};

   constexpr void check_range(Index_t idx) const {
      if (idx < 0 || idx >= len)
         throw OutOfBoundsError(
            "Dynamic_array_view index " + std::to_string(idx)
            + " out of range [0, " + std::to_string(len) + ").");
   }
};

template <typename ElementType>
constexpr Dynamic_array_view<ElementType> make_dynamic_array_view(ElementType* ptr, std::ptrdiff_t len)
{
   return Dynamic_array_view<ElementType>(ptr, len);
}

template <typename ElementType>
constexpr Dynamic_array_view<ElementType> make_dynamic_array_view(ElementType* first, ElementType* last)
{
   return Dynamic_array_view<ElementType>(first, last);
}

template <typename ElementType>
std::ostream& operator<<(std::ostream& ostr, const Dynamic_array_view<ElementType> av)
{
   if (av.empty())
      return ostr;

   ostr << "[";

   for (typename Dynamic_array_view<ElementType>::Index_t i = 0; i < av.size(); i++) {
      ostr << av[i];
      if (i < av.size() - 1)
         ostr << ", ";
   }

   ostr << "]";

   return ostr;
}

} // namespace flexiblesusy

#endif
