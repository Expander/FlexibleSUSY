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

#ifndef TWO_SCALE_SOLVER_H
#define TWO_SCALE_SOLVER_H

#include "rg_flow.hpp"
#include "two_scale_model.hpp"

#include <vector>
#include <string>

class Two_scale;
class Two_scale_model;
class Two_scale_matching;

template<>
class RGFlow<Two_scale> {
public:
   class Error {
   public:
      Error(const std::string& message_) : message(message) {}
      std::string what() const { return message; }
   private:
      std::string message;
   };

   RGFlow();

   void add_matching_condition(const Two_scale_matching*);
   void add_model(Two_scale_model*);
   void run_up();
   void run_down();
   void solve();

private:
   std::vector<Two_scale_model*> models;
   std::vector<const Two_scale_matching*> matching;
   unsigned int maxIterations;

   bool accuracy_goal_reached() const;
   void check_setup() const;
};

typedef RGFlow<Two_scale> Two_scale_solver;

inline RGFlow<Two_scale>::RGFlow()
   : models()
   , matching()
   , maxIterations(10)
{
}

inline void RGFlow<Two_scale>::solve()
{
   check_setup();

   // initial run
   run_up();
   run_down();

   for (unsigned int iter = 0;
        iter < maxIterations && !accuracy_goal_reached();
        ++iter) {
      run_up();
      run_down();
   }
}

inline void RGFlow<Two_scale>::check_setup() const
{
   if (models.size() + 1 != matching.size()) {
      std::string message;
      message = "RGFlow<Two_scale>::Error: number of models does not match"
         "number of matching conditions - 1";
      throw Error(message);
   }
}

inline void RGFlow<Two_scale>::run_up()
{
   for (std::vector<Two_scale_model*>::iterator model = models.begin(),
           end = models.end();
        model != end; ++model) {
      (*model)->run_up();
   }
}

inline void RGFlow<Two_scale>::run_down()
{
   for (std::vector<Two_scale_model*>::iterator model = models.begin(),
           end = models.end();
        model != end; ++model) {
      (*model)->run_down();
   }
}

inline void RGFlow<Two_scale>::add_matching_condition(const Two_scale_matching* mc)
{
   matching.push_back(mc);
}

inline void RGFlow<Two_scale>::add_model(Two_scale_model* m)
{
   models.push_back(m);
}

inline bool RGFlow<Two_scale>::accuracy_goal_reached() const
{
   return true;
}

#endif
