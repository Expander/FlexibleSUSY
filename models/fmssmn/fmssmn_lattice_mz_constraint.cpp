#include "mathdefs.hpp"
#include "consts.hpp"
#include "SM.hpp"
#include "fmssm_oneloop.hpp"
#include "fmssmn_lattice_mz_constraint.hpp"


using namespace std;


Fmssmn_mz_constraint_::Fmssmn_mz_constraint_() :
    fix_scale_to_mz(mZ),
    gcs(),
    ycs()
{}

Fmssmn_mz_constraint::Fmssmn_mz_constraint(double tanBeta) :
    CompoundConstraint<Lattice>::CompoundConstraint
    ({&fix_scale_to_mz, &gcs, &ycs})
{
    gcs.g1 = g1L1(mZ);
    gcs.g2 = g1L2(mZ);
    gcs.g3 = g1L3(mZ);

    Real beta = atan(tanBeta);
    Real vu = vv * sin(beta);
    Real vd = vv * cos(beta);

    CM33 VCKM = standard_VCKM(60*deg);

    CM33 MUMW;
    MUMW << muMW, 0,    0,
	    0,    mcMW, 0,
	    0,    0,    mtMW;
    CM33 MDMW;
    MDMW << mdMW, 0,    0,
	    0,    msMW, 0,
	    0,    0,    mbMW;
    CM33 MEMW;
    MEMW << me,   0,    0,
	    0,    mmu,  0,
	    0,    0,    mtau;

    ycs.Yu = VCKM.transpose() * MUMW / vu;
    ycs.Yd = MDMW / vd;
    ycs.Ye = MEMW / vd;
}
