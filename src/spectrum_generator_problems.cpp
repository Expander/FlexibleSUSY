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

#include "spectrum_generator_problems.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <iostream>

namespace flexiblesusy {

Spectrum_generator_problems::Spectrum_generator_problems(const std::vector<Problems>& problems_)
   : problems(problems_)
{
}

Spectrum_generator_problems::Spectrum_generator_problems(std::vector<Problems>&& problems_)
   : problems(std::move(problems_))
{
}

void Spectrum_generator_problems::set_model_problems(const std::vector<Problems>& p)
{
   problems = p;
}

void Spectrum_generator_problems::set_model_problems(std::vector<Problems>&& p)
{
   problems = std::move(p);
}

void Spectrum_generator_problems::clear()
{
   failed_convergence = false;
}

bool Spectrum_generator_problems::have_problem() const
{
   return no_convergence() ||
      std::any_of(problems.cbegin(), problems.cend(),
                  [] (const Problems& p) { return p.have_problem(); });
}

bool Spectrum_generator_problems::have_warning() const
{
   return std::any_of(problems.cbegin(), problems.cend(),
                      [] (const Problems& p) { return p.have_warning(); });
}

std::vector<std::string> Spectrum_generator_problems::get_problem_strings() const
{
   std::vector<std::string> result;

   for (const auto& p: problems) {
      const auto strings = p.get_problem_strings();
      result.insert(result.end(), strings.cbegin(), strings.cend());
   }

   if (no_convergence())
      result.push_back("no convergence");

   return result;
}

std::vector<std::string> Spectrum_generator_problems::get_warning_strings() const
{
   std::vector<std::string> result;

   for (const auto& p: problems) {
      const auto strings = p.get_warning_strings();
      result.insert(result.end(), strings.cbegin(), strings.cend());
   }

   return result;
}

std::string Spectrum_generator_problems::get_problem_string() const
{
   return concat(get_problem_strings(), '\n');
}

std::string Spectrum_generator_problems::get_warning_string() const
{
   return concat(get_warning_strings(), '\n');
}

void Spectrum_generator_problems::print_problems(std::ostream& ostr) const
{
   if (!have_problem())
      return;

   ostr << get_problem_string();
}

void Spectrum_generator_problems::print_warnings(std::ostream& ostr) const
{
   if (!have_warning())
      return;

   ostr << get_warning_string();
}

const std::vector<Problems>& Spectrum_generator_problems::get_model_problems() const
{
   return problems;
}

std::vector<Problems>& Spectrum_generator_problems::get_model_problems()
{
   return problems;
}

void Spectrum_generator_problems::flag_no_convergence()
{
   failed_convergence = true;
}

void Spectrum_generator_problems::unflag_no_convergence()
{
   failed_convergence = false;
}

bool Spectrum_generator_problems::no_convergence() const
{
   return failed_convergence;
}

std::ostream& operator<<(std::ostream& ostr, const Spectrum_generator_problems& problems)
{
   problems.print_problems(ostr);
   problems.print_warnings(ostr);
   return ostr;
}

} // namespace flexiblesusy
