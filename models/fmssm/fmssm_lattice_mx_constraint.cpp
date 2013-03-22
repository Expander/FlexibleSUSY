#include "fmssm_lattice_mx_constraint.hpp"

Fmssm_mx_constraint::Fmssm_mx_constraint() :
    mxc(),
    mhc(),
    mgc(),
    mfc(),
    tfc(),
    CompoundConstraint<Lattice>::CompoundConstraint
    ({&mxc, &mhc, &mgc, &mfc, &tfc})
{
}
