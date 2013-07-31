#ifndef fmssm_lattice_constraints_hpp
#define fmssm_lattice_constraints_hpp


#include "lattice_foreign_constraint.hpp"
#include "small_matrices.hpp"

namespace flexiblesusy {

#define fortran_fmssm_bc(name)				\
							\
extern "C" void name##_					\
(const Real& g1i, const Real& g2i, const Real& g3i,	\
 const Comp *Yui, const Comp *Ydi, const Comp *Yei,	\
 const Real& m2Hui, const Real& m2Hdi,			\
 const Comp *m2Qi, const Comp *m2Ui, const Comp *m2Di,	\
 const Comp *m2Li, const Comp *m2Ei,			\
 const Comp *Aui, const Comp *Adi, const Comp *Aei,	\
 const Comp& M1i, const Comp& M2i, const Comp& M3i,	\
 const Real& vu, const Real& vd,			\
 const Real& scale0, const Real *x, const int& i,	\
 Real *row, Real *rhs);

fortran_fmssm_bc(fmssm_mx)
fortran_fmssm_bc(fmssm_higgs_masses)
fortran_fmssm_bc(fmssm_gaugino_masses)
fortran_fmssm_bc(fmssm_sfermion_masses)
fortran_fmssm_bc(fmssm_trilinears)
fortran_fmssm_bc(fmssm_ms)
fortran_fmssm_bc(fmssm_gauge_couplings)
fortran_fmssm_bc(fmssm_yukawas)
fortran_fmssm_bc(fmssm_ewsb)

class Fmssm_constraint_on_mx : public ForeignConstraint {
public:
    Fmssm_constraint_on_mx() : ForeignConstraint(1) {}
    void operator()() {
	fmssm_mx_(0,0,0,
		  nullptr,nullptr,nullptr,
		  0,0,
		  nullptr,nullptr,nullptr,
		  nullptr,nullptr,
		  nullptr,nullptr,nullptr,
		  0,0,0,
		  0,0,
		  f->scl0, nullptr, 0,
		  &row[0], &rhs);
	copy_row(0);
    }
};

class Fmssm_constraint_on_ms : public ForeignConstraint {
public:
    Fmssm_constraint_on_ms() : ForeignConstraint(1) {}
    void operator()() {
	set_x();
	fmssm_ms_(0,0,0,
		  nullptr,nullptr,nullptr,
		  0,0,
		  nullptr,nullptr,nullptr,
		  nullptr,nullptr,
		  nullptr,nullptr,nullptr,
		  0,0,0,
		  vu,vd,
		  f->scl0, &x[0], 0,
		  &row[0], &rhs);
	copy_row(0);
    }
    Real vu, vd;
};

class Fmssm_constraint_on_gauge_couplings : public ForeignConstraint {
public:
    Fmssm_constraint_on_gauge_couplings() : ForeignConstraint(3) {}
    void operator()() {
	for (size_t i = 0; i < 3; i++) {
	    fmssm_gauge_couplings_(g1,g2,g3,
				   nullptr,nullptr,nullptr,
				   0,0,
				   nullptr,nullptr,nullptr,
				   nullptr,nullptr,
				   nullptr,nullptr,nullptr,
				   0,0,0,
				   0,0,
				   f->scl0, nullptr, i,
				   &row[0], &rhs);
	    copy_row(i);
	}
    }
    Real g1, g2, g3;
};

class Fmssm_constraint_on_yukawas : public ForeignConstraint {
public:
    Fmssm_constraint_on_yukawas() : ForeignConstraint(54) {}
    void operator()() {
	for (size_t i = 0; i < 54; i++) {
	    fmssm_yukawas_(0,0,0,
			   Yu.data(),Yd.data(),Ye.data(),
			   0,0,
			   nullptr,nullptr,nullptr,
			   nullptr,nullptr,
			   nullptr,nullptr,nullptr,
			   0,0,0,
			   0,0,
			   f->scl0, nullptr, i,
			   &row[0], &rhs);
	    copy_row(i);
	}
    }
    CM33 Yu, Yd, Ye;
};

class Fmssm_constraint_on_ewsb : public ForeignConstraint {
public:
    Fmssm_constraint_on_ewsb() : ForeignConstraint(4) {}
    void operator()() {
	set_x();
	for (size_t i = 0; i < 4; i++) {
	    fmssm_ewsb_(0,0,0,
		       nullptr,nullptr,nullptr,
		       0,0,
		       nullptr,nullptr,nullptr,
		       nullptr,nullptr,
		       nullptr,nullptr,nullptr,
		       0,0,0,
		       vu,vd,
		       f->scl0, &x[0], i,
		       &row[0], &rhs);
	    copy_row(i);
	}
    }
    Real vu, vd;
};

class Fmssm_constraint_on_higgs_masses : public ForeignConstraint {
public:
    Fmssm_constraint_on_higgs_masses() : ForeignConstraint(2) {}
    void operator()() {
	Real m2Hu = (1-f->a)*m2Hu_ini + f->a*m2Hu_fin;
	Real m2Hd = (1-f->a)*m2Hd_ini + f->a*m2Hd_fin;
	for (size_t i = 0; i < 2; i++) {
	    fmssm_higgs_masses_(0,0,0,
				nullptr,nullptr,nullptr,
				m2Hu,m2Hd,
				nullptr,nullptr,nullptr,
				nullptr,nullptr,
				nullptr,nullptr,nullptr,
				0,0,0,
				0,0,
				f->scl0, nullptr, i,
				&row[0], &rhs);
	    copy_row(i);
	}
    }
    Real m2Hu_ini, m2Hu_fin;
    Real m2Hd_ini, m2Hd_fin;
};

class Fmssm_constraint_on_gaugino_masses : public ForeignConstraint {
public:
    Fmssm_constraint_on_gaugino_masses() : ForeignConstraint(6) {}
    void operator()() {
	for (size_t i = 0; i < 6; i++) {
	    fmssm_gaugino_masses_(0,0,0,
				  nullptr,nullptr,nullptr,
				  0,0,
				  nullptr,nullptr,nullptr,
				  nullptr,nullptr,
				  nullptr,nullptr,nullptr,
				  M1,M2,M3,
				  0,0,
				  f->scl0, nullptr, i,
				  &row[0], &rhs);
	    copy_row(i);
	}
    }
    Comp M1, M2, M3;
};

class Fmssm_constraint_on_sfermion_masses : public ForeignConstraint {
public:
    Fmssm_constraint_on_sfermion_masses() : ForeignConstraint(45) {}
    void operator()() {
	for (size_t i = 0; i < 45; i++) {
	    fmssm_sfermion_masses_(0,0,0,
				   nullptr,nullptr,nullptr,
				   0,0,
				   m2Q.data(),m2U.data(),m2D.data(),
				   m2L.data(),m2E.data(),
				   nullptr,nullptr,nullptr,
				   0,0,0,
				   0,0,
				   f->scl0, nullptr, i,
				   &row[0], &rhs);
	    copy_row(i);
	}
    }
    CM33 m2Q, m2U, m2D, m2L, m2E;
};

class Fmssm_constraint_on_trilinears : public ForeignConstraint {
public:
    Fmssm_constraint_on_trilinears() : ForeignConstraint(54) {}
    void operator()() {
	for (size_t i = 0; i < 54; i++) {
	    fmssm_trilinears_(0,0,0,
			      nullptr,nullptr,nullptr,
			      0,0,
			      nullptr,nullptr,nullptr,
			      nullptr,nullptr,
			      Au.data(),Ad.data(),Ae.data(),
			      0,0,0,
			      0,0,
			      f->scl0, nullptr, i,
			      &row[0], &rhs);
	    copy_row(i);
	}
    }
    CM33 Au, Ad, Ae;
};

}

#endif // fmssm_lattice_constraints_hpp
