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

#ifndef RAII_H
#define RAII_H

namespace flexiblesusy {

/**
 * @class RAII_save
 * @brief Saves value of variable and restores it at destruction
 */
template <typename T>
class RAII_save {
public:
   RAII_save(T& var_) noexcept : var(var_), value(var_) {}
   ~RAII_save() { var = value; }

private:
   T& var;
   T value{};
};

template <typename T>
constexpr RAII_save<T> make_raii_save(T& var)
{
   return RAII_save<T>(var);
}

} // namespace flexiblesusy

#endif
