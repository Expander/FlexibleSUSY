#ifndef fmssmn_lattice_numerical_constraints_hpp
#define fmssmn_lattice_numerical_constraints_hpp


#include "lattice_compound_constraint.hpp"


struct Fmssmn_constraint_on_yn_n_ {
    Fmssmn_constraint_on_yn_n_();
    ~Fmssmn_constraint_on_yn_n_();
    std::vector<Constraint<Lattice>*> members;
    CM33 Yn;
    std::vector<std::vector<size_t>> dependence();
};

class Fmssmn_constraint_on_yn_n :
    public Fmssmn_constraint_on_yn_n_,
    public CompoundConstraint<Lattice> {
public:
    Fmssmn_constraint_on_yn_n() : CompoundConstraint(members) {}
};


#endif // fmssmn_lattice_numerical_constraints_hpp
