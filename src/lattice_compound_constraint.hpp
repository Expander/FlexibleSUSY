#ifndef lattice_compound_constraint_hpp
#define lattice_compound_constraint_hpp


#include "compound_constraint.hpp"
#include "lattice_constraint.hpp"

namespace flexiblesusy {

class Lattice;


template<>
class CompoundConstraint<Lattice> : public Constraint<Lattice> {
public:
    CompoundConstraint(std::vector<Constraint<Lattice>*> cs) :
	components(cs)
	{}
    virtual void init(RGFlow<Lattice> *flow, size_t theory, size_t site) {
	Constraint<Lattice>::init(flow, theory, site);
	for (auto c: components) c->init(flow, theory, site);
    }
    virtual void deactivate()
    { for (auto c: components) c->deactivate(); }
    virtual void alloc_rows()
    { for (auto c: components) c->alloc_rows(); }
    virtual void free_rows()
    { for (auto c: components) c->free_rows(); }
    virtual void operator()() { for (auto c: components) (*c)(); }
    virtual void relocate(const std::vector<size_t>& site_map)
    { for (auto c: components) c->relocate(site_map); }
protected:
    virtual void activate() {}
    std::vector<Constraint<Lattice>*> components;
};

template<>
class CompoundMatching<Lattice> : public Matching<Lattice> {
public:
    CompoundMatching(std::vector<Matching<Lattice>*> ms) :
	components(ms)
	{}
    virtual void init(RGFlow<Lattice> *flow, size_t lower_theory) {
	Matching<Lattice>::init(flow, lower_theory);
	for (auto c: components) c->init(flow, lower_theory);
    }
    virtual void deactivate()
    { for (auto c: components) c->deactivate(); }
    virtual void alloc_rows()
    { for (auto c: components) c->alloc_rows(); }
    virtual void free_rows()
    { for (auto c: components) c->free_rows(); }
    virtual void operator()() { for (auto c: components) (*c)(); }
protected:
    virtual void activate() {}
    std::vector<Matching<Lattice>*> components;
};

}

#endif // lattice_compound_constraint_hpp
