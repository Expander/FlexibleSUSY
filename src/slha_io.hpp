// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#ifndef SLHA_IO_H
#define SLHA_IO_H

#include <string>
#include <iosfwd>
#include <functional>
#include <Eigen/Core>
#include <boost/format.hpp>
#include "slhaea.h"

namespace softsusy {
   class QedQcd;
}

class DoubleMatrix;
class ComplexMatrix;

namespace flexiblesusy {

   namespace {
      /// SLHA line formatter for the MASS block entries
      const boost::format mass_formatter(" %9d   %16.8E   # %s\n");
      /// SLHA line formatter for the mixing matrix entries (NMIX, UMIX, VMIX, ...)
      const boost::format mixing_matrix_formatter(" %2d %2d   %16.8E   # %s\n");
      /// SLHA number formatter
      const boost::format number_formatter("%16.8E");
      /// SLHA line formatter for the one-element entries (HMIX, GAUGE, MSOFT, ...)
      const boost::format single_element_formatter(" %5d   %16.8E   # %s\n");
      /// SLHA line formatter for the SPINFO block entries
      const boost::format spinfo_formatter(" %5d   %s\n");
   }

#define PHYSICAL(p) model.get_physical().p
#define MODELPARAMETER(p) model.get_##p()
#define FORMAT_MASS(pdg,mass,name)                                      \
   boost::format(mass_formatter) % pdg % mass % name
#define FORMAT_MIXING_MATRIX(i,k,entry,name)                            \
   boost::format(mixing_matrix_formatter) % i % k % entry % name
#define FORMAT_ELEMENT(pdg,value,name)                                  \
   boost::format(single_element_formatter) % pdg % value % name
#define FORMAT_NUMBER(n)                                                \
   boost::format(number_formatter) % n
#define FORMAT_SPINFO(n,str)                                            \
   boost::format(spinfo_formatter) % n % str

class SLHA_io {
public:
   typedef std::function<void(int, double)> Tuple_processor;
   enum Position { front, back };

   SLHA_io();
   ~SLHA_io() {}

   // reading functions
   void fill(softsusy::QedQcd&) const;
   void read_from_file(const std::string&);
   void read_block(const std::string&, Tuple_processor) const;

   // writing functions
   void set_block(const std::ostringstream&, Position position = back);
   void set_block(const std::string&, double, const std::string&, double scale = 0.);
   void set_block(const std::string&, const Eigen::MatrixXd&, const std::string&, double scale = 0.);
   void set_block(const std::string&, const DoubleMatrix&, const std::string&, double scale = 0.);
   void set_block(const std::string&, const ComplexMatrix&, const std::string&, double scale = 0.);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& = std::cout);

private:
   SLHAea::Coll data;          ///< SHLA data
   static void process_sminputs_tuple(softsusy::QedQcd&, int, double);
};

} // namespace flexiblesusy

#endif
