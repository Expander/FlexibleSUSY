#include "fmssm_lattice_mx_constraint.hpp"


Fmssm_mx_constraint_::Fmssm_mx_constraint_() :
    mxc(),
    mhc(),
    mgc(),
    mfc(),
    tfc()
{}

Fmssm_mx_constraint::Fmssm_mx_constraint() :
    CompoundConstraint<Lattice>::CompoundConstraint
    ({&mxc, &mhc, &mgc, &mfc, &tfc})
{}
