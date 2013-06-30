#ifndef fmssm_lattice_numerical_constraints_hpp
#define fmssm_lattice_numerical_constraints_hpp


#include "lattice_compound_constraint.hpp"
#include "lattice_numerical_constraint.hpp"
#include "small_matrices.hpp"


class Fmssm_constraint_on_ms_n : public NumericalConstraint {
public:
    Fmssm_constraint_on_ms_n() : NumericalConstraint(dependence()[0]) {}
    Real c(const Real *x) const;
    Real vu, vd;
private:
    std::vector<std::vector<size_t>> dependence();
};

struct Fmssm_constraint_on_gauge_couplings_n_ {
    Fmssm_constraint_on_gauge_couplings_n_();
    ~Fmssm_constraint_on_gauge_couplings_n_();
    std::vector<Constraint<Lattice>*> members;
    Real g1, g2, g3;
    std::vector<std::vector<size_t>> dependence();
};

class Fmssm_constraint_on_gauge_couplings_n :
    public Fmssm_constraint_on_gauge_couplings_n_,
    public CompoundConstraint<Lattice> {
public:
    Fmssm_constraint_on_gauge_couplings_n() : CompoundConstraint(members) {}
};

struct Fmssm_constraint_on_yukawas_n_ {
    Fmssm_constraint_on_yukawas_n_();
    ~Fmssm_constraint_on_yukawas_n_();
    std::vector<Constraint<Lattice>*> members;
    CM33 Yu, Yd, Ye;
    std::vector<std::vector<size_t>> dependence();
};

class Fmssm_constraint_on_yukawas_n :
    public Fmssm_constraint_on_yukawas_n_,
    public CompoundConstraint<Lattice> {
public:
    Fmssm_constraint_on_yukawas_n() : CompoundConstraint(members) {}
};

struct Fmssm_constraint_on_ewsb_n_ {
    Fmssm_constraint_on_ewsb_n_();
    ~Fmssm_constraint_on_ewsb_n_();
    std::vector<Constraint<Lattice>*> members;
    Real vu, vd;
    std::vector<std::vector<size_t>> dependence();
};

class Fmssm_constraint_on_ewsb_n :
    public Fmssm_constraint_on_ewsb_n_,
    public CompoundConstraint<Lattice> {
public:
    Fmssm_constraint_on_ewsb_n() : CompoundConstraint(members) {}
};


#endif // fmssm_lattice_numerical_constraints_hpp
