#include "consts.hpp"
#include "fmssmn_lattice_msusy_constraint.hpp"


namespace flexiblesusy {

Fmssmn_msusy_constraint_::Fmssmn_msusy_constraint_() :
    msc(),
    ewsb()
{}

Fmssmn_msusy_constraint::Fmssmn_msusy_constraint(double tanBeta) :
    CompoundConstraint<Lattice>::CompoundConstraint({&msc, &ewsb})
{
    Real beta = atan(tanBeta);
    Real vu = vv * sin(beta);
    Real vd = vv * cos(beta);

    msc.vu = ewsb.vu = vu;
    msc.vd = ewsb.vd = vd;
}

}
