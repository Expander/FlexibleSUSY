#include "fmssmn_lattice_mx_constraint.hpp"


Fmssmn_mx_constraint_::Fmssmn_mx_constraint_() :
    mxc(),
    ync(),
    mhc(),
    mgc(),
    mfc(),
    tfc()
{}

Fmssmn_mx_constraint::Fmssmn_mx_constraint() :
    CompoundConstraint<Lattice>::CompoundConstraint
    ({&mxc, &ync, &mhc, &mgc, &mfc, &tfc})
{}
