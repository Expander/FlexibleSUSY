#ifndef fmssm_fmssmn_lattice_numerical_matchings_hpp
#define fmssm_fmssmn_lattice_numerical_matchings_hpp


#include "lattice_compound_constraint.hpp"


namespace flexiblesusy {

struct Fmssm_fmssmn_gauge_couplings_n_ {
    Fmssm_fmssmn_gauge_couplings_n_();
    ~Fmssm_fmssmn_gauge_couplings_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_gauge_couplings_n :
    public Fmssm_fmssmn_gauge_couplings_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_gauge_couplings_n() : CompoundMatching(members) {}
};

struct Fmssm_fmssmn_yukawas_n_ {
    Fmssm_fmssmn_yukawas_n_();
    ~Fmssm_fmssmn_yukawas_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_yukawas_n :
    public Fmssm_fmssmn_yukawas_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_yukawas_n() : CompoundMatching(members) {}
};

struct Fmssm_fmssmn_mu_b_n_ {
    Fmssm_fmssmn_mu_b_n_();
    ~Fmssm_fmssmn_mu_b_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_mu_b_n :
    public Fmssm_fmssmn_mu_b_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_mu_b_n() : CompoundMatching(members) {}
};

struct Fmssm_fmssmn_gaugino_masses_n_ {
    Fmssm_fmssmn_gaugino_masses_n_();
    ~Fmssm_fmssmn_gaugino_masses_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_gaugino_masses_n :
    public Fmssm_fmssmn_gaugino_masses_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_gaugino_masses_n() : CompoundMatching(members) {}
};

struct Fmssm_fmssmn_higgs_masses_n_ {
    Fmssm_fmssmn_higgs_masses_n_();
    ~Fmssm_fmssmn_higgs_masses_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_higgs_masses_n :
    public Fmssm_fmssmn_higgs_masses_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_higgs_masses_n() : CompoundMatching(members) {}
};

struct Fmssm_fmssmn_sfermion_masses_n_ {
    Fmssm_fmssmn_sfermion_masses_n_();
    ~Fmssm_fmssmn_sfermion_masses_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_sfermion_masses_n :
    public Fmssm_fmssmn_sfermion_masses_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_sfermion_masses_n() : CompoundMatching(members) {}
};

struct Fmssm_fmssmn_trilinears_n_ {
    Fmssm_fmssmn_trilinears_n_();
    ~Fmssm_fmssmn_trilinears_n_();
    std::vector<Matching<Lattice>*> members;
    std::vector<std::vector<std::vector<size_t>>> dependence();
};

class Fmssm_fmssmn_trilinears_n :
    public Fmssm_fmssmn_trilinears_n_,
    public CompoundMatching<Lattice> {
public:
    Fmssm_fmssmn_trilinears_n() : CompoundMatching(members) {}
};

}

#endif // fmssm_fmssmn_lattice_numerical_matchings_hpp
