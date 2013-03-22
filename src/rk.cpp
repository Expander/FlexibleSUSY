
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
#include <iostream>

#include "rk.hpp"

using std::cout;
using std::endl;
using std::flush;
using std::function;

/// Returns |a| with sign of b in front
inline Real sign(Real a, Real b) 
{ return ((b) >= 0.0 ? fabs(a) : -fabs(a)); }

// returns >0 if there's a problem:
int integrateOdes(BRVec& ystart, Real from, Real to, Real eps,
	      Real h1, Real hmin, 
	      Derivs derivs,
	      function<int(BRVec& y, const BRVec& dydx, Real
			   *x, Real htry, Real eps, BRVec& yscal,
			   Real *hdid, Real *hnext, 
			   Derivs derivs)> rkqs) {  
  int nvar =  ystart.size();
  int nstp, i;
  Real x, hnext, hdid, h;
  BRVec yscal(nvar), y(ystart), dydx;
  
  x = from;
  h = sign(h1, to - from);
  
  const int MAXSTP = 400;
  const Real TINY = 1.0e-16;

  for (nstp = 1; nstp <= MAXSTP; nstp++) {
    dydx = derivs(x, y);
    for (i = 0; i < nvar; i++)
      yscal(i) = fabs(y(i)) + fabs(dydx(i) * h) + TINY;
    if ((x + h - to) * (x + h - from) > 0.0) h = to - x;
    int smallStep = rkqs(y, dydx, &x, h, eps, yscal, &hdid, &hnext, derivs);
    if (smallStep) return 1;

    if ((x - to) * (to - from) >= 0.0) {
      for (i = 0; i< nvar; i++) ystart(i) = y(i);
      return 0;
    }
      
    if (fabs(hnext) <= hmin) {
      nstp = MAXSTP; // bail out
      {
	cout << "Step size too small in diffeq.cpp:integrateOdes\n";
	cout << "**********x = " << x << "*********\n";
	for (i = 0;i < nvar;i++) 
	  cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	    ") = " << dydx(i) << endl;
	cout.flush();
      }
    }
    
    h = hnext;
  }
  
  {
    cout << "Bailed out of diffeq.cpp:too many steps in integrateOdes\n";
    cout << "**********x = " << x << "*********\n";
    for (i = 0;i < nvar;i++) 
      cout << "y(" << i << ") = " << y(i) << " dydx(" << i <<
	") = " << dydx(i) << endl;
    cout.flush();
  }
  
  return 1;
}

int odeStepper(BRVec& y, const BRVec& dydx, Real *x, Real
		htry, Real eps, BRVec& yscal, Real *hdid, 
		Real *hnext,		
		Derivs derivs)
{
  const Real SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;

  int i, n = y.size();
  Real errmax, h, htemp, xnew;
  
  BRVec yerr(n), ytemp(n);
  h = htry;
  for (;;) {
    rungeKuttaStep(y, dydx, *x, h, ytemp, yerr, derivs);
    errmax = 0.0;
    for (i = 0; i < n;i++) errmax = std::max(errmax, fabs(yerr(i) / yscal(i)));
    errmax  /= eps;
    if (errmax <= 1.0) break;
    htemp = SAFETY * h * pow(errmax, PSHRNK);
    h = (h >= 0.0 ? std::max(htemp ,0.1 * h) : std::min(htemp, 0.1 * h));
    xnew = (*x) + h;
    if (xnew == *x) 
      {
	{
	cout << "At x = " << *x;
	cout << ",stepsize underflow in odeStepper" << flush << endl;
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

void rungeKuttaStep(const BRVec& y, const BRVec& dydx,
	     Real x, Real h, BRVec& yout, BRVec& yerr,
	     Derivs derivs) {
  int i;
  const Real a2 = 0.2,a3 = 0.3,a4 = 0.6,a5 = 1.0,a6 = 0.875,b21 =
    0.2,b31 = 3.0 / 40.0,b32 = 9.0 / 40.0,b41 = 0.3,b42 = -0.9,b43 = 1.2,
    b51 = -11.0 / 54.0, b52 = 2.5,b53 = -70.0 / 27.0,b54 = 35.0 / 27.0,
    b61 = 1631.0 / 55296.0,b62 = 175.0 / 512.0,b63 = 575.0 / 13824.0,
    b64 = 44275.0 / 110592.0,b65 = 253.0 / 4096.0,c1 = 37.0 / 378.0,
    c3 = 250.0 / 621.0,c4 = 125.0 / 594.0,c6 = 512.0 / 1771.0,
    dc5 = -277.00 / 14336.0;
  const Real dc1 = c1-2825.0 / 27648.0,dc3 = c3-18575.0 / 48384.0,
    dc4 = c4-13525.0 / 55296.0,dc6 = c6-0.25;
  
  int n = y.size();
  
  BRVec ytemp = y + b21 * h * dydx;
  BRVec ak2 = derivs(x + a2 * h, ytemp);

  // Allowing piece-wise calculating of ytemp for speed reasons
  for (i = 0; i < n; i++)
    ytemp(i) = y(i) + h * (b31 * dydx(i) + b32 * ak2(i));
  BRVec ak3 = derivs(x + a3 * h, ytemp);

  for (i = 0; i < n; i++)
    ytemp(i) = y(i) + h * (b41 * dydx(i) + b42 * ak2(i) + b43 * ak3(i));
  BRVec ak4 = derivs(x+a4*h,ytemp);

  for (i = 0; i < n; i++)
    ytemp(i) = y(i) + h * (b51 * dydx(i) + b52 * ak2(i) + b53
				   * ak3(i) + b54 * ak4(i));
  BRVec ak5 = derivs(x + a5 * h, ytemp);

  for (i = 0; i < n; i++)
    ytemp(i) = y(i) + h * (b61 * dydx(i) + b62 * ak2(i) + b63
				   * ak3(i) + b64 * ak4(i) + b65 * ak5(i));
  BRVec ak6 = derivs(x + a6 * h, ytemp);

  for (i = 0; i < n; i++)
    yout(i) = y(i) + h * (c1 * dydx(i) + c3 * ak3(i) + c4 *
				  ak4(i) + c6 * ak6(i));
  for (i = 0; i < n; i++)
    yerr(i) = h * (dc1 * dydx(i) + dc3 * ak3(i) + 
		   dc4 * ak4(i) + dc5 * ak5(i) + dc6 * ak6(i));
}
