#ifndef fmssmn_lattice_hpp
#define fmssmn_lattice_hpp


#include "fmssmn.hpp"
#include "lattice_solver.hpp"

namespace flexiblesusy {

class Lattice;

template<>
class Fmssmn<Lattice>: public Lattice_model {
public:
    Fmssmn();
    virtual ~Fmssmn() {}
    virtual void calculate_spectrum();
    virtual std::string name() const { return "FMSSMN"; }
    virtual void print(std::ostream& s) const;

    Real dx(const Real a, const Real *x, size_t i) const;
    void ddx(const Real a, const Real *x, size_t i, Real *ddx) const;

    struct Translator : public RGFlow<Lattice>::Translator {
	Translator(RGFlow<Lattice> *f, size_t T, size_t m) :
	    RGFlow<Lattice>::Translator::Translator(f, T, m) {}
	#include "models/fmssmn/fmssmn_lattice_translator.inc"
    };

    Translator operator()(size_t m) const { return Translator(f, T, m); }
};

}

#endif // fmssmn_lattice_hpp
