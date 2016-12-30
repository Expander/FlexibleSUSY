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

#include "logger.hpp"
#include "config.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <map>
#include <vector>

namespace flexiblesusy {

/**
 * @class Problems
 * @brief stores problem flags for the spectrum generator
 */
template <unsigned Number_of_particles, unsigned Number_of_parameters>
class Problems {
public:
   enum : unsigned {
      NUMBER_OF_PARTICLES = Number_of_particles,
      NUMBER_OF_PARAMETERS = Number_of_parameters
   };

   Problems(const std::array<std::string, Number_of_particles>&,
            const std::array<std::string, Number_of_parameters>&);

   void flag_bad_mass(unsigned particle, bool flag = true)
      { bad_masses.at(particle) = flag; }
   void flag_running_tachyon(unsigned particle, bool flag = true) {
      running_tachyons.at(particle) = flag;
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
      if (flag) WARNING("running " << particle_names[particle] << " tachyon");
#endif
   }
   void flag_pole_tachyon(unsigned particle, bool flag = true) {
      pole_tachyons.at(particle) = flag;
#if defined(ENABLE_VERBOSE) || defined(ENABLE_DEBUG)
      if (flag) WARNING("pole " << particle_names[particle] << " tachyon");
#endif
   }
   void flag_thrown(const std::string& msg = "") {
      thrown = true;
      exception_msg = msg;
   }
   void flag_no_ewsb()
      { failed_ewsb = true; }
   void flag_no_convergence()
      { failed_convergence = true; }
   void flag_no_perturbative()
      { non_perturbative = true; }
   void flag_no_pole_mass_convergence(unsigned particle)
      { failed_pole_mass_convergence.at(particle) = true; }
   void flag_non_perturbative_parameter(int parameter, double value, double scale, double threshold = 0.) {
      if (parameter < -1 || parameter >= static_cast<int>(Number_of_parameters))
         ERROR("Parameter index " << parameter << " out of range ["
               << -1 <<  ", " << (Number_of_parameters - 1) << "]");
      non_pert_pars[parameter] = NonPerturbativeValue(value, scale, threshold);
   }
   void flag_no_rho_convergence()
      { failed_rho_convergence = true; }

   void unflag_bad_mass(unsigned particle)
      { bad_masses.at(particle) = false; }
   void unflag_running_tachyon(unsigned particle)
      { running_tachyons.at(particle) = false; }
   void unflag_pole_tachyon(unsigned particle)
      { pole_tachyons.at(particle) = false; }
   void unflag_all_tachyons() {
      running_tachyons = std::array<bool, Number_of_particles>{};
      pole_tachyons    = std::array<bool, Number_of_particles>{};
   }
   void unflag_thrown()
      { thrown = false; exception_msg = ""; }
   void unflag_no_ewsb()
      { failed_ewsb = false; }
   void unflag_no_convergence()
      { failed_convergence = false; }
   void unflag_no_perturbative()
      { non_perturbative = false; }
   void unflag_no_pole_mass_convergence(unsigned particle)
      { failed_pole_mass_convergence.at(particle) = false; }
   void unflag_non_perturbative_parameter(int parameter)
      { non_pert_pars.erase(parameter); }
   void unflag_no_rho_convergence() { failed_rho_convergence = false; }

   bool is_bad_mass(unsigned particle) const
      { return bad_masses.at(particle); }
   bool is_running_tachyon(unsigned particle) const
      { return running_tachyons.at(particle); }
   bool is_pole_tachyon(unsigned particle) const
      { return pole_tachyons.at(particle); }
   bool have_bad_mass() const
      { return std::any_of(bad_masses.cbegin(), bad_masses.cend(), [](bool x){ return x; }); }
   bool have_running_tachyon() const {
      return std::any_of(running_tachyons.cbegin(), running_tachyons.cend(), [](bool x){ return x; });
   }
   bool have_pole_tachyon() const {
      return std::any_of(pole_tachyons.cbegin(), pole_tachyons.cend(), [](bool x){ return x; });
   }
   bool have_tachyon() const
      { return have_running_tachyon() || have_pole_tachyon(); }
   bool have_thrown() const
      { return thrown; }
   bool have_non_perturbative_parameter() const
      { return !non_pert_pars.empty(); }
   bool have_failed_pole_mass_convergence() const
      { return std::any_of(failed_pole_mass_convergence.cbegin(), failed_pole_mass_convergence.cend(), [](bool x){ return x; }); }
   bool no_ewsb() const            { return failed_ewsb; }
   bool no_convergence() const     { return failed_convergence; }
   bool no_perturbative() const    { return non_perturbative; }
   bool no_rho_convergence() const { return failed_rho_convergence; }

   void clear();                      ///< clear all problems
   bool have_problem() const;         ///< problems which yield invalid spectrum
   bool have_warning() const;         ///< warnings
   std::vector<std::string> get_problem_strings() const;
   std::vector<std::string> get_warning_strings() const;
   std::string get_problem_string() const { return concat(get_problem_strings(), '\n'); }
   std::string get_warning_string() const { return concat(get_warning_strings(), '\n'); }
   void print_problems(std::ostream& = std::cout) const;
   void print_warnings(std::ostream& = std::cout) const;

   std::array<bool, Number_of_particles> get_bad_masses() const
      { return bad_masses; }
   std::array<bool, Number_of_particles> get_running_tachyons() const
      { return running_tachyons; }
   std::array<bool, Number_of_particles> get_pole_tachyons() const
      { return pole_tachyons; }
   std::array<bool, Number_of_particles> get_failed_pole_mass_convergence() const
      { return failed_pole_mass_convergence; }

private:
   struct NonPerturbativeValue {
      NonPerturbativeValue() = default;
      NonPerturbativeValue(double value_, double scale_, double threshold_)
         : value(value_), scale(scale_), threshold(threshold_) {}
      double value{0.}, scale{0.}, threshold{0.};
   };

   std::array<bool, Number_of_particles> bad_masses; ///< imprecise mass eigenvalues
   std::array<bool, Number_of_particles> running_tachyons; ///< tachyonic particles (running mass)
   std::array<bool, Number_of_particles> pole_tachyons; ///< tachyonic particles (pole mass)
   std::array<bool, Number_of_particles> failed_pole_mass_convergence; ///< no convergence during pole mass calculation
   std::map<int, NonPerturbativeValue> non_pert_pars; ///< non-perturbative parmeters
   std::string exception_msg;          ///< exception message
   std::array<std::string, Number_of_particles> particle_names; ///< particle names
   std::array<std::string, Number_of_parameters> parameter_names; ///< parameter names
   bool thrown;                        ///< excepton thrown
   bool failed_ewsb;                   ///< no EWSB
   bool failed_convergence;            ///< no convergence
   bool non_perturbative;              ///< non-perturbative running
   bool failed_rho_convergence;        ///< rho-parameter not converged

   std::string get_parameter_name(int) const; ///< returns parameter name
   static std::string concat(const std::vector<std::string>&, char); ///< concatenate strings
};

template <unsigned Number_of_particles, unsigned Number_of_parameters>
Problems<Number_of_particles, Number_of_parameters>::Problems(
   const std::array<std::string, Number_of_particles>& particle_names_,
   const std::array<std::string, Number_of_parameters>& parameter_names_)
   : bad_masses() // intializes all elements to zero (= false)
   , running_tachyons()
   , pole_tachyons()
   , failed_pole_mass_convergence()
   , non_pert_pars()
   , exception_msg("")
   , particle_names(particle_names_)
   , parameter_names(parameter_names_)
   , thrown(false)
   , failed_ewsb(false)
   , failed_convergence(false)
   , non_perturbative(false)
   , failed_rho_convergence(false)
{
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
void Problems<Number_of_particles, Number_of_parameters>::clear()
{
   bad_masses = std::array<bool, Number_of_particles>{};
   running_tachyons = std::array<bool, Number_of_particles>{};
   pole_tachyons = std::array<bool, Number_of_particles>{};
   failed_pole_mass_convergence = std::array<bool, Number_of_particles>{};
   non_pert_pars.clear();
   exception_msg = "";
   failed_ewsb = false;
   failed_convergence = false;
   non_perturbative = false;
   failed_rho_convergence = false;
   thrown = false;
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
bool Problems<Number_of_particles, Number_of_parameters>::have_problem() const
{
   return have_tachyon() || failed_ewsb || failed_convergence
      || non_perturbative || failed_rho_convergence || thrown
      || have_failed_pole_mass_convergence()
      || have_non_perturbative_parameter();
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
bool Problems<Number_of_particles, Number_of_parameters>::have_warning() const
{
   return have_bad_mass();
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
std::string Problems<Number_of_particles, Number_of_parameters>::get_parameter_name(int idx) const
{
   if (idx == -1)
      return "Q";

   return parameter_names.at(idx);
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
std::vector<std::string> Problems<Number_of_particles, Number_of_parameters>::get_problem_strings() const
{
   std::vector<std::string> strings;

   for (unsigned i = 0; i < Number_of_particles; ++i) {
      if (running_tachyons[i])
         strings.push_back("running tachyon " + particle_names[i]);
   }
   for (unsigned i = 0; i < Number_of_particles; ++i) {
      if (pole_tachyons[i])
         strings.push_back("pole tachyon " + particle_names[i]);
   }
   if (failed_ewsb)
      strings.push_back("no ewsb");
   if (failed_convergence)
      strings.push_back("no convergence");
   if (non_perturbative)
      strings.push_back("non-perturbative");
   if (failed_rho_convergence)
      strings.push_back("no rho convergence");
   if (thrown)
      strings.push_back("exception thrown(" + exception_msg + ")");
   for (unsigned i = 0; i < Number_of_particles; ++i) {
      if (failed_pole_mass_convergence[i])
         strings.push_back("no M" + particle_names[i] + " pole convergence");
   }

   for (const auto& par: non_pert_pars) {
      const std::string par_name = get_parameter_name(par.first);
      std::string str("non-perturbative " + par_name);
      if (par.second.threshold > 0) {
         str += " [|" + par_name + "|(" +
                std::to_string(par.second.scale) + ") = " +
                std::to_string(par.second.value) +
                " > " + std::to_string(par.second.threshold) + "]";
      } else {
         str += " [" + par_name + "(" +
                std::to_string(par.second.scale) +
                ") = " + std::to_string(par.second.value) + "]";
      }
      strings.push_back(str);
   }

   return strings;
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
std::vector<std::string> Problems<Number_of_particles, Number_of_parameters>::get_warning_strings() const
{
   std::vector<std::string> strings;

   for (unsigned i = 0; i < Number_of_particles; ++i) {
      if (bad_masses[i])
         strings.push_back("Warning: imprecise M" + particle_names[i]);
   }

   return strings;
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
void Problems<Number_of_particles, Number_of_parameters>::print_problems(std::ostream& ostr) const
{
   if (!have_problem())
      return;

   ostr << get_problem_string();
}


template <unsigned Number_of_particles, unsigned Number_of_parameters>
void Problems<Number_of_particles, Number_of_parameters>::print_warnings(std::ostream& ostr) const
{
   if (!have_warning())
      return;

   ostr << get_warning_string();
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
std::string Problems<Number_of_particles, Number_of_parameters>::concat(
   const std::vector<std::string>& strings, char separator)
{
   std::string result;

   for (const auto& s: strings)
      result += s + separator;

   return result;
}

template <unsigned Number_of_particles, unsigned Number_of_parameters>
std::ostream& operator<<(std::ostream& ostr, const Problems<Number_of_particles, Number_of_parameters>& problems)
{
   problems.print_problems(ostr);
   problems.print_warnings(ostr);
   return ostr;
}

} // namespace flexiblesusy

#endif
