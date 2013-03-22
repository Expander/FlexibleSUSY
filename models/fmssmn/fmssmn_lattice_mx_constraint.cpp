#include "fmssmn_lattice_mx_constraint.hpp"

Fmssmn_mx_constraint::Fmssmn_mx_constraint() :
    mxc(),
    ync(),
    mhc(),
    mgc(),
    mfc(),
    tfc(),
    CompoundConstraint<Lattice>::CompoundConstraint
    ({&mxc, &ync, &mhc, &mgc, &mfc, &tfc})
{
}
