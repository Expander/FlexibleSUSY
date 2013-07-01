
/** \file numerics.cpp
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html

   $Log: numerics.cpp,v $
   Revision 1.4  2006/04/25 11:32:22  allanach
   Fixed B22, B0, B1 - previous expression switched too soon to p=0 
   expression, leading to convergence failures where there shouldn't have been
   any (in particular, at high m0 and m12 for instance).

   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.8  2005/07/15 15:10:47  allanach
   Added analytic dilog routine

   Revision 1.7  2005/06/16 13:57:04  allanach
   Added a Cauchy-distributed random number generator

   Revision 1.6  2005/05/30 14:22:14  allanach
   Fixed stau mixing in ISAWIG interface to 7.64

   Revision 1.5  2005/05/13 16:07:27  allanach
   Edited precision due to double precision

   Revision 1.4  2005/04/13 14:53:54  allanach
   Corrected binning procedure so it does what you expect

   Revision 1.3  2005/04/12 10:44:20  allanach
   Added bin function to calculate the bin number of some data

   Revision 1.2  2005/04/11 14:06:36  allanach
   Added random number routines

   Revision 1.18  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.17  2003/08/19 14:26:22  allanach
   Changing lowOrg to be more sensible about gauge unification. Should now be
   called with POSITIVE mgut and a flag for gauge unification.

   Revision 1.16  2003/07/28 12:11:37  allanach
   More error trapping, and rearranging rpvsoftsusy to use correct Higgs VEV
   (which is sometimes called at MZ)

   Revision 1.15  2003/07/25 13:39:15  allanach
   Trapped errors properly rather than exiting

   Revision 1.14  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.13  2003/03/25 17:03:05  allanach
   Added extra case for b1(0,m1,0,q)

   Revision 1.12  2003/02/21 13:02:07  allanach
   Changed headings to new conventions

   Revision 1.10  2002/12/19 18:26:40  allanach
   Fixed numerical rounding error in d27 function

   Revision 1.7  2002/10/14 17:14:29  allanach
   Bug-fixed c0 function

   Revision 1.6  2002/09/09 10:42:54  allanach
   TOLERANCE replaces EPS as being more user-friendly

   Revision 1.5  2002/09/03 14:16:44  allanach
   Taken PRINTOUT, MIXING, TOLERANCE out of def.h to make it quicker to
   compile once they are changed.

   Revision 1.4  2002/07/30 12:57:31  allanach
   SOFTSUSY1.5

   Revision 1.3  2002/04/18 14:32:05  allanach
   Changed RGEs and anomalous dimensions to be compatible with new notation;
   started implementation of rewsb in R-parity violation

   Revision 1.6  2001/10/04 19:26:34  allanach
   New version deals with AMSB correctly

   Revision 1.4  2001/09/28 13:50:13  allanach
   More careful treatment of limits in Passarino-Veltman functions b0, b1.

   Revision 1.3  2001/08/08 09:52:33  allanach
   Added dilogarithm function - could be speeded up....

   Revision 1.2  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#include <algorithm>
#include "logger.hpp"
#include "rk.hpp"

using namespace std;
using namespace Eigen;

namespace runge_kutta {

/// Returns |a| with sign of b in front
double sign(double a, double b) 
{ return b >= 0 ? fabs(a) : -fabs(a); }

// returns >0 if there's a problem:
int integrateOdes(ArrayXd& ystart, double from, double to, double eps,
		  double h1, double hmin, Derivs derivs,
		  RungeKuttaQuinticStepper rkqs) {  
  int nvar =  ystart.size();
  int nstp, i;
  double x, hnext, hdid, h;
  ArrayXd yscal(nvar), y(ystart), dydx;
  
  x = from;
  h = sign(h1, to - from);
  
  const int MAXSTP = 400;
  const double TINY = 1.0e-16;

  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    dydx = derivs(x, y);
    yscal = y.abs() + (dydx * h).abs() + TINY;
    if ((x + h - to) * (x + h - from) > 0.0) h = to - x;
    int smallStep = rkqs(y, dydx, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (smallStep) return 1;

    if ((x - to) * (to - from) >= 0.0) {
      ystart = y;
      return 0;
    }
      
    if (fabs(hnext) <= hmin) {
      nstp = MAXSTP; // bail out
      {
	ERROR("Step size too small in rk.cpp:integrateOdes\n"
	      << "**********x = " << x << "*********");
	for (i = 0; i < nvar; i++) 
	    ERROR("y(" << i << ") = " << y(i) << " dydx(" << i <<
		  ") = " << dydx(i));
      }
    }
    
    h = hnext;
  }
  
  {
    ERROR("Bailed out of rk.cpp:too many steps in integrateOdes\n"
	    << "**********x = " << x << "*********");
    for (i = 0; i < nvar; i++) 
	ERROR("y(" << i << ") = " << y(i) << " dydx(" << i <<
	      ") = " << dydx(i));
  }
  
  return 1;
}

int odeStepper(ArrayXd& y, const ArrayXd& dydx, double *x, double htry,
	       double eps, ArrayXd& yscal, double *hdid, double *hnext,
	       Derivs derivs)
{
  const double SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;

  int n = y.size();
  double errmax, h, htemp, xnew;
  
  ArrayXd yerr(n), ytemp(n);
  h = htry;
  for (;;) {
    rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, derivs);
    errmax = (yerr / yscal).abs().maxCoeff();
    errmax  /= eps;
    if (errmax <= 1.0) break;
    htemp = SAFETY * h * pow(errmax, PSHRNK);
    h = (h >= 0.0 ? max(htemp ,0.1 * h) : min(htemp, 0.1 * h));
    xnew = (*x) + h;
    if (xnew == *x) 
      {
	{
	  ERROR("At x = " << *x
		<< ",stepsize underflow in odeStepper");
	}
	return 1;
      }
  }
  if (errmax > ERRCON) *hnext = SAFETY * h * pow(errmax,PGROW);
  else *hnext = 5.0 * h;
  *x += (*hdid = h);
  y = ytemp;
  return 0;
}

void rungeKuttaStep(const ArrayXd& y, const ArrayXd& dydx, double x,
		    double h, ArrayXd& yout, ArrayXd& yerr, Derivs derivs)
{
  const double a2 = 0.2,a3 = 0.3,a4 = 0.6,a5 = 1.0,a6 = 0.875,b21 =
    0.2,b31 = 3.0 / 40.0,b32 = 9.0 / 40.0,b41 = 0.3,b42 = -0.9,b43 = 1.2,
    b51 = -11.0 / 54.0, b52 = 2.5,b53 = -70.0 / 27.0,b54 = 35.0 / 27.0,
    b61 = 1631.0 / 55296.0,b62 = 175.0 / 512.0,b63 = 575.0 / 13824.0,
    b64 = 44275.0 / 110592.0,b65 = 253.0 / 4096.0,c1 = 37.0 / 378.0,
    c3 = 250.0 / 621.0,c4 = 125.0 / 594.0,c6 = 512.0 / 1771.0,
    dc5 = -277.00 / 14336.0;
  const double dc1 = c1-2825.0 / 27648.0,dc3 = c3-18575.0 / 48384.0,
    dc4 = c4-13525.0 / 55296.0,dc6 = c6-0.25;
  
  ArrayXd ytemp = b21 * h * dydx + y;
  ArrayXd ak2 = derivs(x + a2 * h, ytemp);

  // Allowing piece-wise calculating of ytemp for speed reasons
  ytemp = y + h * (b31 * dydx + b32 * ak2);
  ArrayXd ak3 = derivs(x + a3 * h, ytemp);

  ytemp = y + h * (b41 * dydx + b42 * ak2 + b43 * ak3);
  ArrayXd ak4 = derivs(x+a4*h,ytemp);

  ytemp = y + h * (b51 * dydx + b52 * ak2 + b53 * ak3 + b54 * ak4);
  ArrayXd ak5 = derivs(x + a5 * h, ytemp);

  ytemp = y + h * (b61 * dydx + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5);
  ArrayXd ak6 = derivs(x + a6 * h, ytemp);

  yout = y + h * (c1 * dydx + c3 * ak3 + c4 * ak4 + c6 * ak6);
  yerr = h * (dc1 * dydx + dc3 * ak3 + dc4 * ak4 + dc5 * ak5 + dc6 * ak6);
}

} // namespace runge_kutta
