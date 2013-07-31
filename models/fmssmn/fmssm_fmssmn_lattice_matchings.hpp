#ifndef fmssm_fmssmn_lattice_matchings.hpp
#define fmssm_fmssmn_lattice_matchings.hpp


#include "lattice_foreign_constraint.hpp"


namespace flexiblesusy {

#define decl_fmssm_fmssmn_mc(name)					\
									\
extern "C" void name##_							\
(const Real& scale0, const Real *w, const Real *x, const int& i,	\
 Real *row, Real *rhs);

decl_fmssm_fmssmn_mc(fmssm_fmssmn_gauge_couplings)
decl_fmssm_fmssmn_mc(fmssm_fmssmn_yukawas)
decl_fmssm_fmssmn_mc(fmssm_fmssmn_mu_b)
decl_fmssm_fmssmn_mc(fmssm_fmssmn_gaugino_masses)
decl_fmssm_fmssmn_mc(fmssm_fmssmn_higgs_masses)
decl_fmssm_fmssmn_mc(fmssm_fmssmn_sfermion_masses)
decl_fmssm_fmssmn_mc(fmssm_fmssmn_trilinears)

class Fmssm_fmssmn_gauge_couplings : public ForeignMatching {
public:
    Fmssm_fmssmn_gauge_couplings() : ForeignMatching(3) {}
    void operator()() {
	for (size_t i = 0; i < 3; i++) {
	    mssm_mssmrhn_gauge_couplings_
	    (f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_yukawas : public ForeignMatching {
public:
    Fmssm_fmssmn_yukawas() : ForeignMatching(54) {}
    void operator()() {
	for (size_t i = 0; i < 54; i++) {
	    mssm_mssmrhn_yukawas_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_mu_b : public ForeignMatching {
public:
    Fmssm_fmssmn_mu_b() : ForeignMatching(4) {}
    void operator()() {
	for (size_t i = 0; i < 4; i++) {
	    mssm_mssmrhn_mu_b_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_gaugino_masses : public ForeignMatching {
public:
    Fmssm_fmssmn_gaugino_masses() : ForeignMatching(6) {}
    void operator()() {
	for (size_t i = 0; i < 6; i++) {
	    mssm_mssmrhn_gaugino_masses_
	    (f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_higgs_masses : public ForeignMatching {
public:
    Fmssm_fmssmn_higgs_masses() : ForeignMatching(2) {}
    void operator()() {
	for (size_t i = 0; i < 2; i++) {
	    mssm_mssmrhn_higgs_masses_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_sfermion_masses : public ForeignMatching {
public:
    Fmssm_fmssmn_sfermion_masses() : ForeignMatching(45) {}
    void operator()() {
	for (size_t i = 0; i < 45; i++) {
	    mssm_mssmrhn_sfermion_masses_
	    (f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

class Fmssm_fmssmn_trilinears : public ForeignMatching {
public:
    Fmssm_fmssmn_trilinears() : ForeignMatching(54) {}
    void operator()() {
	for (size_t i = 0; i < 54; i++) {
	    mssm_mssmrhn_trilinears_(f->scl0,nullptr,nullptr,i,&row[0],&rhs);
	    copy_row(i);
	}
    }
};

}

#endif // fmssm_fmssmn_lattice_matchings.hpp
