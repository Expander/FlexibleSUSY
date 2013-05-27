#include <vector>
#include "lattice_numerical_constraint.hpp"
#include "fmssm_fmssmn_lattice_numerical_matchings.hpp"


using namespace std;


#define fortran_fmssm_fmssmn_mc_n(name)					\
									\
extern "C" double name##_n_						\
(const Real& scale0, const Real *w, const Real *x, const int& i);

fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_gauge_couplings)
fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_yukawas)
fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_mu_b)
fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_gaugino_masses)
fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_higgs_masses)
fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_sfermion_masses)
fortran_fmssm_fmssmn_mc_n(fmssm_fmssmn_trilinears)


Fmssm_fmssmn_gauge_couplings_n_::Fmssm_fmssmn_gauge_couplings_n_() :
    members(3)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 3; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_gauge_couplings_n_(self->f->scl0,w,x,i);
	    });
    }
}

Fmssm_fmssmn_gauge_couplings_n_::~Fmssm_fmssmn_gauge_couplings_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}

Fmssm_fmssmn_yukawas_n_::Fmssm_fmssmn_yukawas_n_() :
    members(54)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 54; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_yukawas_n_(self->f->scl0,w,x,i);
	    });
    }
}

Fmssm_fmssmn_yukawas_n_::~Fmssm_fmssmn_yukawas_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}

Fmssm_fmssmn_mu_b_n_::Fmssm_fmssmn_mu_b_n_() :
    members(4)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 4; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_mu_b_n_(self->f->scl0,w,x,i);
	    });
    }
}
    
Fmssm_fmssmn_mu_b_n_::~Fmssm_fmssmn_mu_b_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}

Fmssm_fmssmn_gaugino_masses_n_::Fmssm_fmssmn_gaugino_masses_n_() :
    members(6)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 6; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_gaugino_masses_n_(self->f->scl0,w,x,i);
	    });
    }
}

Fmssm_fmssmn_gaugino_masses_n_::~Fmssm_fmssmn_gaugino_masses_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}

Fmssm_fmssmn_higgs_masses_n_::Fmssm_fmssmn_higgs_masses_n_() :
    members(2)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 2; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_higgs_masses_n_(self->f->scl0,w,x,i);
	    });
    }
}

Fmssm_fmssmn_higgs_masses_n_::~Fmssm_fmssmn_higgs_masses_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}

Fmssm_fmssmn_sfermion_masses_n_::Fmssm_fmssmn_sfermion_masses_n_() :
    members(45)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 45; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_sfermion_masses_n_(self->f->scl0,w,x,i);
	    });
    }
}

Fmssm_fmssmn_sfermion_masses_n_::~Fmssm_fmssmn_sfermion_masses_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}

Fmssm_fmssmn_trilinears_n_::Fmssm_fmssmn_trilinears_n_() :
    members(54)
{
    vector<vector<vector<size_t>>> nonzeros(dependence());
    for (size_t i = 0; i < 54; i++) {
	members[i] = new AnyNumericalMatching(nonzeros[i][0],nonzeros[i][1],
	    [&,i](const AnyNumericalMatching *self,
		  const Real *w, const Real *x) {
		return fmssm_fmssmn_trilinears_n_(self->f->scl0,w,x,i);
	    });
    }
}

Fmssm_fmssmn_trilinears_n_::~Fmssm_fmssmn_trilinears_n_()
{
    for (auto m: members)
	delete static_cast<AnyNumericalMatching *>(m);
}
