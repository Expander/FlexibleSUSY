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

#ifndef MSSMD5O_MSSMRHN_TWO_SCALE_MATCHING_H
#define MSSMD5O_MSSMRHN_TWO_SCALE_MATCHING_H

#include "MSSMD5O_MSSMRHN_matching.hpp"
#include "MSSMD5O_input_parameters.hpp"
#include "two_scale_matching.hpp"

namespace flexiblesusy {

template<class T>
class MSSMD5O;

template<class T>
class MSSMRHN;

class Two_scale;

template<>
class MSSMD5O_MSSMRHN_matching<Two_scale> : public Matching<Two_scale> {
public:
    MSSMD5O_MSSMRHN_matching();
    MSSMD5O_MSSMRHN_matching(const MSSMD5O_input_parameters&);
    void match_low_to_high_scale_model();
    void match_high_to_low_scale_model();
    double get_scale() const;
    void set_models(Two_scale_model *lower, Two_scale_model *upper);
    double get_initial_scale_guess() const;
    void set_lower_input_parameters(const MSSMD5O_input_parameters&);
    void set_scale(double);
    void reset();

private:
    double scale;
    double initial_scale_guess;
    double fixed_scale;
    MSSMD5O<Two_scale> *lower;
    MSSMRHN<Two_scale> *upper;
    MSSMD5O_input_parameters inputPars;

    void make_initial_scale_guess();
    void update_scale();
    void invert_seesaw_formula
    (const Eigen::Matrix3d& WOp, const Eigen::Vector3d& YvDiag,
     Eigen::Matrix3d& Yv, Eigen::Matrix3d& Mv);
};

}

#endif
