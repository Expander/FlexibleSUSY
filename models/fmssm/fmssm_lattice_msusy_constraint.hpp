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

#ifndef FMSSM_LATTICE_MSUSY_CONSTRAINT_H
#define FMSSM_LATTICE_MSUSY_CONSTRAINT_H


#include "lattice_compound_constraint.hpp"
#include "fmssm_lattice_constraints.hpp"

namespace flexiblesusy {

// auxiliary class for initializing own members before the base class
// see http://www.boost.org/doc/libs/1_53_0/libs/utility/base_from_member.html
struct Fmssm_msusy_constraint_ {
    Fmssm_msusy_constraint_();
    Fmssm_constraint_on_ms msc;
    Fmssm_constraint_on_ewsb ewsb;
};

class Fmssm_msusy_constraint :
    public Fmssm_msusy_constraint_,
    public CompoundConstraint<Lattice> {
public:
    Fmssm_msusy_constraint(double tanBeta);
};

}

#endif
