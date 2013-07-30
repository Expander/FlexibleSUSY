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

#include "error.hpp"

#include <sstream>

namespace flexiblesusy {

TachyonError::TachyonError(const Two_scale_model* model_,
                           const std::string& particle_name_,
                           int particle_index_)
   : model(model_)
   , particle_name(particle_name_)
   , particle_index(particle_index_)
{
}

std::string TachyonError::what() const
{
   std::stringstream message;
   message << "TachyonError: tachyonic particle detected: "
           << particle_name << "(" << particle_index << ")";
   return message.str();
}

NoEWSBError::NoEWSBError(const Two_scale_model* model_,
                         double requested_precision_)
   : model(model_)
   , requested_precision(requested_precision_)
{
}

std::string NoEWSBError::what() const
{
   std::stringstream message;
   message << "NoEWSBError: Could not find a solution to the EWSB"
      " equations with requested precision " << requested_precision;
   return message.str();
}

}
