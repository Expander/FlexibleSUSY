#ifndef lattice_compound_constraint_hpp
#define lattice_compound_constraint_hpp


#include "compound_constraint.hpp"
#include "lattice_constraint.hpp"


class Lattice;


template<>
class CompoundConstraint<Lattice> : public Constraint<Lattice> {
public:
    CompoundConstraint(std::vector<Constraint<Lattice>*> cs) :
	constraints(cs)
	{}
    virtual void init(RGFlow<Lattice> *flow, size_t theory, size_t site)
    { for (auto c: constraints) c->init(flow, theory, site); }
    virtual void alloc_rows()
    { for (auto c: constraints) c->alloc_rows(); }
    virtual void free_rows()
    { for (auto c: constraints) c->free_rows(); }
    virtual void operator()() { for (auto c: constraints) (*c)(); }
protected:
    std::vector<Constraint<Lattice>*> constraints;
};

template<>
class CompoundMatching<Lattice> : public Matching<Lattice> {
public:
    CompoundMatching(std::vector<Matching<Lattice>*> ms) :
	matchings(ms)
	{}
    virtual void init(RGFlow<Lattice> *flow, size_t lower_theory)
    { for (auto c: matchings) c->init(flow, lower_theory); }
    virtual void alloc_rows()
    { for (auto c: matchings) c->alloc_rows(); }
    virtual void free_rows()
    { for (auto c: matchings) c->free_rows(); }
    virtual void operator()() { for (auto c: matchings) (*c)(); }
protected:
    std::vector<Matching<Lattice>*> matchings;
};


#endif // lattice_compound_constraint_hpp
