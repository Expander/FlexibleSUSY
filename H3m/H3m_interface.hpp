#ifndef H3m_interface_HPP
#define H3m_interface_HPP

#include <complex>
#include <iostream>
#include <Eigen/Core>

namespace h3m {

typedef Eigen::Matrix<double,2,1> V2;
typedef Eigen::Matrix<double,2,2> RM22;
typedef Eigen::Matrix<double,3,3> RM33;

struct Parameters {
   // DR-bar parameters
   double scale{};         ///< renormalization scale
   double mu{};            ///< mu parameter
   double g3{};            ///< gauge coupling g3 SU(3)
   double vd{};            ///< VEV of down Higgs
   double vu{};            ///< VEV of up Higgs
   RM33 mq2{RM33::Zero()}; ///< soft-breaking squared left-handed squark mass parameters
   RM33 md2{RM33::Zero()}; ///< soft-breaking squared right-handed down-squark mass parameters
   RM33 mu2{RM33::Zero()}; ///< soft-breaking squared right-handed up-squark mass parameters
   double At{};
   double Ab{};

   // DR-bar masses
   double MG{};            ///< gluino
   double MW{};            ///< W
   double MZ{};            ///< Z
   double Mt{};            ///< top-quark
   double Mb{};            ///< down-quark
   double MA{};
   V2 MSt{V2::Zero()};     ///< stops
   V2 MSb{V2::Zero()};     ///< sbottoms

   // DR-bar mixing angles
   double s2t{};	   ///< sine of 2 times the stop mixing angle
   double s2b{};	   ///< sine of 2 times the sbot mixing angle
};

inline std::ostream& operator<<(std::ostream& ostr, const Parameters& pars)
{
   ostr <<
      "H3m parameters:\n"
      "  Q   = " << pars.scale << '\n' <<
      "  mu  = " << pars.mu << '\n' <<
      "  g3  = " << pars.g3 << '\n' <<
      "  vd  = " << pars.vd << '\n' <<
      "  vu  = " << pars.vu << '\n' <<
      "  mq2 = " << pars.mq2 << '\n' <<
      "  md2 = " << pars.md2 << '\n' <<
      "  mu2 = " << pars.mu2 << '\n' <<
      "  At  = " << pars.At << '\n' <<
      "  Ab  = " << pars.Ab << '\n' <<
      "  MG  = " << pars.MG << '\n' <<
      "  MW  = " << pars.MW << '\n' <<
      "  MZ  = " << pars.MZ << '\n' <<
      "  Mt  = " << pars.Mt << '\n' <<
      "  Mb  = " << pars.Mb << '\n' <<
      "  MA  = " << pars.MA << '\n' <<
      "  MSt = " << pars.MSt.transpose() << '\n' <<
      "  MSb = " << pars.MSb.transpose() << '\n' <<
      "  s2t = " << pars.s2t << '\n' <<
      "  s2b = " << pars.s2b << '\n';

   return ostr;
}

}	//	h3m

#endif	//	H3m_interface_HPP
