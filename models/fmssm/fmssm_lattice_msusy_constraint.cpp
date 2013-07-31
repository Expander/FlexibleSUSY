#include "consts.hpp"
#include "fmssm_lattice_msusy_constraint.hpp"


namespace flexiblesusy {

Fmssm_msusy_constraint_::Fmssm_msusy_constraint_() :
    msc(),
    ewsb()
{}

Fmssm_msusy_constraint::Fmssm_msusy_constraint(double tanBeta) :
    CompoundConstraint<Lattice>::CompoundConstraint({&msc, &ewsb})
{
    Real beta = atan(tanBeta);
    Real vu = vv * sin(beta);
    Real vd = vv * cos(beta);

    msc.vu = ewsb.vu = vu;
    msc.vd = ewsb.vd = vd;
}

}
