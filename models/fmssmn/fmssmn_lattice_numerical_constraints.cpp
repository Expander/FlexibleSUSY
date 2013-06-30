#include <vector>
#include "lattice_numerical_constraint.hpp"
#include "fmssmn_lattice_numerical_constraints.hpp"


using namespace std;


#define fortran_fmssmn_bc_n(name)					\
									\
extern "C" double name##_n_						\
(const Real& g1i, const Real& g2i, const Real& g3i,			\
 const Comp *Yui, const Comp *Ydi, const Comp *Yni, const Comp *Yei,	\
 const Real& m2Hui, const Real& m2Hdi,					\
 const Comp *m2Qi, const Comp *m2Ui, const Comp *m2Di,			\
 const Comp *m2Li, const Comp *m2Ni, const Comp *m2Ei,			\
 const Comp *Aui, const Comp *Adi, const Comp *Ani, const Comp *Aei,	\
 const Comp& M1i, const Comp& M2i, const Comp& M3i,			\
 const Real& vu, const Real& vd,					\
 const Real& scale0, const Real *x, const int& i);

fortran_fmssmn_bc_n(fmssmn_yn)


Fmssmn_constraint_on_yn_n_::Fmssmn_constraint_on_yn_n_() :
    members(18)
{
    vector<vector<size_t>> nonzeros(dependence());
    for (size_t i = 0; i < 18; i++) {
	members[i] = new AnyNumericalConstraint(nonzeros[i],
	    [&,i](const AnyNumericalConstraint *self, const Real *x) {
		return fmssmn_yn_n_(0,0,0,
				    nullptr,nullptr,Yn.data(),nullptr,
				    0,0,
				    nullptr,nullptr,nullptr,
				    nullptr,nullptr,nullptr,
				    nullptr,nullptr,nullptr,nullptr,
				    0,0,0,
				    0,0,
				    self->f->scl0, x, i);
	    });
    }
}

Fmssmn_constraint_on_yn_n_::~Fmssmn_constraint_on_yn_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalConstraint *>(m);
}
