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

#ifndef FMSSMN_LATTICE_MZ_CONSTRAINT_H
#define FMSSMN_LATTICE_MZ_CONSTRAINT_H

#include "lattice_compound_constraint.hpp"
#include "fmssmn_lattice_constraints.hpp"

class Fmssmn_mz_constraint : public CompoundConstraint<Lattice> {
public:
    Fmssmn_mz_constraint(double tanBeta);

    class Fixed_mz : public AnySingleSiteConstraint {
    public:
	Fixed_mz() :
	    AnySingleSiteConstraint(1,
		[=](AnySingleSiteConstraint *dummy) {
		    A(0,0) = u(0);
		    z(0) = std::log(mZ/f->scl0);
		})
	    {}
    };

    Fixed_mz fixed_mz;
    Fmssmn_constraint_on_gauge_couplings gcs;
    Fmssmn_constraint_on_yude ycs;
};

#endif // FMSSMN_LATTICE_MZ_CONSTRAINT_H
