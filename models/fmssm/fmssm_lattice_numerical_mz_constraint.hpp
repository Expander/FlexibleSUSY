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

#ifndef fmssm_lattice_numerical_mz_constraint_hpp
#define fmssm_lattice_numerical_mz_constraint_hpp

#include "lattice_compound_constraint.hpp"
#include "fmssm_lattice_numerical_constraints.hpp"


// auxiliary class for initializing own members before the base class
// see http://www.boost.org/doc/libs/1_53_0/libs/utility/base_from_member.html
struct Fmssm_mz_constraint_n_ {
    Fmssm_mz_constraint_n_();
    Fixed_t fix_scale_to_mz;
    Fmssm_constraint_on_gauge_couplings_n gcs;
    Fmssm_constraint_on_yukawas_n ycs;
};

class Fmssm_mz_constraint_n :
    public Fmssm_mz_constraint_n_,
    public CompoundConstraint<Lattice> {
public:
    Fmssm_mz_constraint_n(double tanBeta);
};

#endif // fmssm_lattice_numerical_mz_constraint_hpp
