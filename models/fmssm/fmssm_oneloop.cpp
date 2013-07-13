#include "mathdefs.hpp"
#include "consts.hpp"
#include "fmssm_oneloop.hpp"


using namespace std;


const double b1 = -33/5.0;
const double b2 = -1;
const double b3 = 3;


double g1L(double tmt0, double g0, double b)
{
    return 1 / sqrt(1/sqr(g0) + tmt0*b/(8*sqr(pi)));
}

double g1L1(double mu)
{
    return g1L(log(mu/mZ), g1mZ, b1);
}

double g1L2(double mu)
{
    return g1L(log(mu/mZ), g2mZ, b2);
}

double g1L3(double mu)
{
    return g1L(log(mu/mtMW), g3mt, b3);
}

double MX1L()
{
    return mZ * exp((8*sqr(pi)*(1/sqr(g2mZ)-1/sqr(g1mZ)))/(b1-b2));
}
