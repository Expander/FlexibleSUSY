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
#include "two_scale_matching.hpp"

#include <vector>
#include <string>

class Two_scale;

template<>
class RGFlow<Two_scale> {
public:
   class Error {
   public:
      Error(const std::string& message_) : message(message_) {}
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
   std::vector<const Two_scale_matching*> matching_condition;
   unsigned int maxIterations;

   bool accuracy_goal_reached() const;
   void check_setup() const;
};

typedef RGFlow<Two_scale> Two_scale_solver;

inline RGFlow<Two_scale>::RGFlow()
   : models()
   , matching_condition()
   , maxIterations(10)
{
}

inline void RGFlow<Two_scale>::solve()
{
   check_setup();

   // initial run
   run_up();
   run_down();

   unsigned int iter = maxIterations;
   while (iter-- && !accuracy_goal_reached()) {
      run_up();
      run_down();
   }
}

inline void RGFlow<Two_scale>::check_setup() const
{
   if (models.size() != matching_condition.size() + 1) {
      std::string message;
      message = "RGFlow<Two_scale>::Error: number of models does not match "
         "number of matching conditions + 1";
      throw Error(message);
   }

   for (std::vector<Two_scale_model*>::const_iterator model = models.begin(),
           end = models.end(); model != end; ++model) {
      if (!*model) {
         std::string message;
         message = "RGFlow<Two_scale>::Error: model pointer is NULL";
         throw Error(message);
      }
   }

   for (std::vector<const Two_scale_matching*>::const_iterator
           mc = matching_condition.begin(),
           end = matching_condition.end(); mc != end; ++mc) {
      if (!*mc) {
         std::string message;
         message = "RGFlow<Two_scale>::Error: matching condition "
            "pointer is NULL";
         throw Error(message);
      }
   }
}

inline void RGFlow<Two_scale>::run_up()
{
   const std::size_t number_of_models = models.size();
   std::size_t i = 0;
   while (i != number_of_models) {
      // init model parameters from low-scale model
      if (i > 0) {
         const Two_scale_model* low_scale_model = models[i - 1];
         const Two_scale_matching* mc = matching_condition[i - 1];
         models[i]->setParameters(
            mc->calcHighFromLowScaleParameters(
               low_scale_model->getParameters()));
      }
      models[i++]->run_up();
   }
}

inline void RGFlow<Two_scale>::run_down()
{
   const std::size_t number_of_models = models.size();
   std::size_t i = number_of_models;
   while (i--) {
      // init model parameters from high-scale model
      if (i < number_of_models - 1) {
         const Two_scale_model* high_scale_model = models[i + 1];
         const Two_scale_matching* mc = matching_condition[i];
         models[i]->setParameters(
            mc->calcLowFromHighScaleParameters(
               high_scale_model->getParameters()));
      }
      models[i]->run_down();
   }
}

inline void RGFlow<Two_scale>::add_matching_condition(const Two_scale_matching* mc)
{
   matching_condition.push_back(mc);
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
