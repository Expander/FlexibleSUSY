#include "fmssm_lattice.hpp"


extern "C" double  dxfmssm_(const double& a, const double *x, const int& i);
extern "C" double ddxfmssm_(const double& a, const double *x, const int& i,
			    double *ddx);


namespace flexiblesusy {

Fmssm<Lattice>::Fmssm() : Lattice_model(169)
{
}

Real Fmssm<Lattice>::dx(const Real a, const Real *x, size_t i) const
{
    int i_ = i;
    return dxfmssm_(a, x, i_);
}

void Fmssm<Lattice>::ddx(const Real a, const Real *x, size_t i, Real *ddx)const
{
    int i_ = i;
    ddxfmssm_(a, x, i_, ddx);
}

void Fmssm<Lattice>::calculate_spectrum()
{
}

void Fmssm<Lattice>::print(std::ostream& /* s */) const
{
}

}
