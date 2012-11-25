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

#ifndef RUNNING_COUPLING_H
#define RUNNING_COUPLING_H

#include <vector>
#include <cmath>
#include <algorithm>

#include "linalg.h"
#include "logger.hpp"

/**
 * @class Running_coupling
 * @brief stores the gauge couplings at different scales
 *
 * Usage:
 *
 *    Running_coupling rc;
 *    rc.run(sm, 100, 1.e12, 50, true); // sm is a two scale model
 *    rc.write_to_file("running_coupling.dat");
 */
class Running_coupling {
public:
   typedef std::pair<double, DoubleVector> TTouple;///< touple of scale and couplings

   Running_coupling();
   virtual ~Running_coupling();

   template <class T>
   void run(T, double, double, unsigned int, bool include_endpoint = false);
   TTouple get_max_scale() const;
   void reset();
   void write_to_file(const string&) const;

protected:
   typedef std::vector<TTouple> TData; ///< container for the scales and couplings
   struct TDataComp {
      bool operator() (const TData::value_type& i,const TData::value_type& j) {
         return i.first < j.first;
      }
   };

   /// write a comment line
   void write_comment_line(char, std::ofstream&, std::size_t, int) const;

private:
   TData couplings;
};

/**
 * Add running couplings between scale q1 and q2.
 *
 * @param rge class with RGEs and parameters
 * @param q1 scale to start at
 * @param q2 end scale
 * @param number_of_steps number of steps
 * @param include_endpoint include the endpoint q2 in the running
 *        (false by default)
 */
template <class T>
void Running_coupling::run(T rge, double q1, double q2,
                           unsigned int number_of_steps, bool include_endpoint)
{
   if (q1 <= 0.0 || q2 <= 0.0) {
      ERROR("negative scales are not allowed: q1=" << q1 << ", q2=" << q2);
      return;
   }

   if (number_of_steps < 1)
      number_of_steps = 1;

   // if the endpoint should be included, the scale loop must run from
   // (n == 0) to (n == number_of_steps); otherwise it runs from (n == 0) to (n
   // == number_of_steps - 1)
   const unsigned int endpoint_offset = include_endpoint ? 1 : 0;

   // run from q1 to q2
   for (unsigned int n = 0; n < number_of_steps + endpoint_offset; ++n) {
      const double scale = exp(log(q1) + n * (log(q2) - log(q1)) / number_of_steps);
      rge.run_to(scale);
      couplings.push_back(TData::value_type(scale, rge.displayGauge()));
   }

   std::sort(couplings.begin(), couplings.end(), TDataComp());
}

#endif
