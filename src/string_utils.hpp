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

#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <complex>
#include <string>
#include <vector>

namespace flexiblesusy {

/// concatenate strings
std::string concat(const std::vector<std::string>&);

/// concatenate strings with separator
template <typename T>
std::string concat(const std::vector<std::string>& strings, const T& separator)
{
   std::string result;

   for (const auto& s: strings)
      result += s + separator;

   return result;
}

std::string to_string(char);
std::string to_string(unsigned char);
std::string to_string(unsigned short);
std::string to_string(unsigned int);
std::string to_string(unsigned long);
std::string to_string(unsigned long long);
std::string to_string(signed char);
std::string to_string(signed short);
std::string to_string(signed int);
std::string to_string(signed long);
std::string to_string(signed long long);
std::string to_string(double);
std::string to_string(const std::complex<double>&);

} // namespace flexiblesusy

#endif
