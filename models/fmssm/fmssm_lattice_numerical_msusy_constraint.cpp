#include "consts.hpp"
#include "fmssm_lattice_numerical_msusy_constraint.hpp"


Fmssm_msusy_constraint_n_::Fmssm_msusy_constraint_n_() :
    msc(),
    ewsb()
{}

Fmssm_msusy_constraint_n::Fmssm_msusy_constraint_n(double tanBeta) :
    CompoundConstraint<Lattice>::CompoundConstraint({&msc, &ewsb})
{
    Real beta = atan(tanBeta);
    Real vu = vv * sin(beta);
    Real vd = vv * cos(beta);

    msc.vu = ewsb.vu = vu;
    msc.vd = ewsb.vd = vd;
}
