#include <vector>
#include "fmssm_lattice_numerical_constraints.hpp"


using namespace std;


#define fortran_fmssm_bc_n(name)			\
							\
extern "C" double name##_n_				\
(const Real& g1i, const Real& g2i, const Real& g3i,	\
 const Comp *Yui, const Comp *Ydi, const Comp *Yei,	\
 const Real& m2Hui, const Real& m2Hdi,			\
 const Comp *m2Qi, const Comp *m2Ui, const Comp *m2Di,	\
 const Comp *m2Li, const Comp *m2Ei,			\
 const Comp *Aui, const Comp *Adi, const Comp *Aei,	\
 const Comp& M1i, const Comp& M2i, const Comp& M3i,	\
 const Real& vu, const Real& vd,			\
 const Real& scale0, const Real *x, const int& i);

fortran_fmssm_bc_n(fmssm_ms)
fortran_fmssm_bc_n(fmssm_gauge_couplings)
fortran_fmssm_bc_n(fmssm_yukawas)
fortran_fmssm_bc_n(fmssm_ewsb)


Real Fmssm_constraint_on_ms_n::c(const Real *x) const
{
    return fmssm_ms_n_(0,0,0,
		       nullptr,nullptr,nullptr,
		       0,0,
		       nullptr,nullptr,nullptr,
		       nullptr,nullptr,
		       nullptr,nullptr,nullptr,
		       0,0,0,
		       vu,vd,
		       f->scl0, x, 0);
}

Fmssm_constraint_on_gauge_couplings_n_::
Fmssm_constraint_on_gauge_couplings_n_() :
    members(3)
{
    vector<vector<size_t>> nonzeros(dependence());
    for (size_t i = 0; i < 3; i++) {
	members[i] = new AnyNumericalConstraint(nonzeros[i],
	    [&,i](const AnyNumericalConstraint *self, const Real *x) {
		return fmssm_gauge_couplings_n_(g1,g2,g3,
						nullptr,nullptr,nullptr,
						0,0,
						nullptr,nullptr,nullptr,
						nullptr,nullptr,
						nullptr,nullptr,nullptr,
						0,0,0,
						0,0,
						self->f->scl0, x, i);
	    });
    }
}

Fmssm_constraint_on_gauge_couplings_n_::
~Fmssm_constraint_on_gauge_couplings_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalConstraint *>(m);
}

Fmssm_constraint_on_yukawas_n_::Fmssm_constraint_on_yukawas_n_() :
    members(54)
{
    vector<vector<size_t>> nonzeros(dependence());
    for (size_t i = 0; i < 54; i++) {
	members[i] = new AnyNumericalConstraint(nonzeros[i],
	    [&,i](const AnyNumericalConstraint *self, const Real *x) {
		return fmssm_yukawas_n_(0,0,0,
					Yu.data(),Yd.data(),Ye.data(),
					0,0,
					nullptr,nullptr,nullptr,
					nullptr,nullptr,
					nullptr,nullptr,nullptr,
					0,0,0,
					0,0,
					self->f->scl0, x, i);
	    });
    }
}

Fmssm_constraint_on_yukawas_n_::~Fmssm_constraint_on_yukawas_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalConstraint *>(m);
}

Fmssm_constraint_on_ewsb_n_::Fmssm_constraint_on_ewsb_n_() :
    members(4)
{
    vector<vector<size_t>> nonzeros(dependence());
    for (size_t i = 0; i < 4; i++) {
	members[i] = new AnyNumericalConstraint(nonzeros[i],
	    [&,i](const AnyNumericalConstraint *self, const Real *x) {
		return fmssm_ewsb_n_(0,0,0,
				     nullptr,nullptr,nullptr,
				     0,0,
				     nullptr,nullptr,nullptr,
				     nullptr,nullptr,
				     nullptr,nullptr,nullptr,
				     0,0,0,
				     vu,vd,
				     self->f->scl0, x, i);
	    });
    }
}

Fmssm_constraint_on_ewsb_n_::~Fmssm_constraint_on_ewsb_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalConstraint *>(m);
}
