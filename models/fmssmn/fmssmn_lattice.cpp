#include "fmssmn_lattice.hpp"


extern "C" double  dxfmssmn_(const double& a, const double *x, const int& i);
extern "C" double ddxfmssmn_(const double& a, const double *x, const int& i,
			     double *ddx);


namespace flexiblesusy {

Fmssmn<Lattice>::Fmssmn() : Lattice_model(214)
{
}

Real Fmssmn<Lattice>::dx(const Real a, const Real *x, size_t i) const
{
    int i_ = i;
    return dxfmssmn_(a, x, i_);
}

void Fmssmn<Lattice>::ddx(const Real a, const Real *x, size_t i, Real *ddx)const
{
    int i_ = i;
    ddxfmssmn_(a, x, i_, ddx);
}

void Fmssmn<Lattice>::calculate_spectrum()
{
}

void Fmssmn<Lattice>::print(std::ostream& /* s */) const
{
}

}
