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

#ifndef MSSMCBS_TWO_SCALE_H
#define MSSMCBS_TWO_SCALE_H

#include "MSSMcbs_model.hpp"
#include "CMSSM_two_scale_model.hpp"

namespace flexiblesusy {

template<>
class MSSMcbs<Two_scale> : public CMSSM<Two_scale> {
public:
   explicit MSSMcbs(const CMSSM_input_parameters& input_ = CMSSM_input_parameters());
   virtual ~MSSMcbs();

   virtual Eigen::ArrayXd beta() const;
   CMSSM_soft_parameters calc_beta() const;
};

} // namespace flexiblesusy

#endif
