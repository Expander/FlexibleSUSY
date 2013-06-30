/** \file numerics.h
   - Project:     SOFTSUSY 
   - Author:      Ben Allanach 
   - Manual:      hep-ph/0104145, Comp. Phys. Comm. 143 (2002) 305 
   - Webpage:     http://allanach.home.cern.ch/allanach/softsusy.html
   - Description: Integration of ODEs by Runge Kutta, minimum finding and
                derivative calculation

   $Log: numerics.h,v $
   Revision 1.3  2005/11/09 14:12:24  allanach
   Updated for SOFTSUSY 2.0.1 - cleaned everything up etc

   Revision 1.7  2005/07/15 15:10:47  allanach
   Added analytic dilog routine

   Revision 1.6  2005/06/16 13:57:04  allanach
   Added a Cauchy-distributed random number generator

   Revision 1.5  2005/05/13 16:07:27  allanach
   Edited precision due to double precision

   Revision 1.4  2005/04/13 14:53:55  allanach
   Corrected binning procedure so it does what you expect

   Revision 1.3  2005/04/12 10:44:20  allanach
   Added bin function to calculate the bin number of some data

   Revision 1.2  2005/04/11 14:06:36  allanach
   Added random number routines

   Revision 1.1.1.1  2004/11/19 16:18:31  allanach


   Revision 1.7  2004/01/15 13:54:54  allanach
   New heaer style implemented

   Revision 1.6  2003/08/14 09:25:49  allanach
   Used new standard convention for included files

   Revision 1.5  2003/05/20 15:19:40  allanach
   doxygen comment style implemented

   Revision 1.4  2002/10/30 16:22:48  allanach
   Fixed bug in f function

   Revision 1.3  2002/10/01 11:52:16  allanach
   Bug-fixed bound-finding routines. MGUT always determined.

   Revision 1.4  2001/08/08 09:52:33  allanach
   Added dilogarithm function - could be speeded up....

   Revision 1.3  2001/07/18 14:42:51  allanach
   Added proper header info
*/

#ifndef RK_HPP
#define RK_HPP

#include <functional>
#include <Eigen/Dense>

namespace runge_kutta {

typedef std::function<Eigen::ArrayXd(double, const Eigen::ArrayXd&)> Derivs;

/// A single step of Runge Kutta (5th order), input: 
/// y and dydx (derivative of y), x is independent variable. yout is value
/// after step. derivs is a user-supplied function
void rungeKuttaStep(const Eigen::ArrayXd& y, const Eigen::ArrayXd& dydx,
	double x, double h, Eigen::ArrayXd& yout, Eigen::ArrayXd& yerr,
	Derivs derivs);

/// organises the variable step-size for Runge-Kutta evolution
int odeStepper(Eigen::ArrayXd& y, const Eigen::ArrayXd& dydx, double *x,
	       double htry, double eps, Eigen::ArrayXd& yscal, double *hdid,
	       double *hnext, Derivs derivs);

typedef std::function<decltype(odeStepper)> RungeKuttaQuinticStepper;

/// Organises integration of 1st order system of ODEs
int integrateOdes(Eigen::ArrayXd& ystart, double x1, double x2, double eps,
		  double h1, double hmin, Derivs derivs,
		  RungeKuttaQuinticStepper rkqs);

} // namespace runge_kutta

#endif // RK_HPP
