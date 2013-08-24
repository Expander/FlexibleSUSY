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
#include "logger.hpp"
#include "error.hpp"

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
   struct Modsel {
      double parameter_output_scale; ///< key = 12
      Modsel() : parameter_output_scale(0.) {}
   };

   class ReadError : public Error {
   public:
      ReadError(const std::string& message_) : message(message_) {}
      virtual ~ReadError() {}
      virtual std::string what() const { return message; }
   private:
      std::string message;
   };

   SLHA_io();
   ~SLHA_io() {}

   // reading functions
   void fill(softsusy::QedQcd&) const;
   const Modsel& get_modsel() const { return modsel; }
   void read_from_file(const std::string&);
   void read_block(const std::string&, const Tuple_processor&) const;
   template <class Derived>
   void read_block(const std::string&, Eigen::MatrixBase<Derived>&) const;
   double read_entry(const std::string&, int) const;
   void read_modsel();

   // writing functions
   void set_data(const SLHAea::Coll& data_) { data = data_; }
   void set_block(const std::ostringstream&, Position position = back);
   void set_block(const std::string&, double, const std::string&, double scale = 0.);
   template <class Derived>
   void set_block(const std::string&, const Eigen::MatrixBase<Derived>&, const std::string&, double scale = 0.);
   void set_block(const std::string&, const DoubleMatrix&, const std::string&, double scale = 0.);
   void set_block(const std::string&, const ComplexMatrix&, const std::string&, double scale = 0.);
   void write_to_file(const std::string&);
   void write_to_stream(std::ostream& = std::cout);

private:
   SLHAea::Coll data;          ///< SHLA data
   Modsel modsel;              ///< data from block MODSEL
   static void process_sminputs_tuple(softsusy::QedQcd&, int, double);
   static void process_modsel_tuple(Modsel&, int, double);
};

template <class Derived>
void SLHA_io::read_block(const std::string& block_name, Eigen::MatrixBase<Derived>& matrix) const
{
   if (data.find(block_name) == data.cend()) {
      WARNING("block " << block_name << " not found");
      return;
   }

   const int cols = matrix.cols(), rows = matrix.rows();

   for (SLHAea::Block::const_iterator line = data.at(block_name).cbegin(),
        end = data.at(block_name).cend(); line != end; ++line) {
      if (!line->is_data_line())
         continue;

      if (line->size() >= 3) {
         const int i = SLHAea::to<int>((*line)[0]) - 1;
         const int k = SLHAea::to<int>((*line)[1]) - 1;
         if (0 <= i && i < rows && 0 <= k && k < cols) {
            const double value = SLHAea::to<double>((*line)[2]);
            matrix(i,k) = value;
         }
      } else {
         WARNING(block_name << " entry has less than 3 columns");
      }
   }
}

template <class Derived>
void SLHA_io::set_block(const std::string& name,
                        const Eigen::MatrixBase<Derived>& matrix,
                        const std::string& symbol, double scale)
{
   std::ostringstream ss;
   ss << "Block " << name;
   if (scale != 0.)
      ss << " Q= " << FORMAT_NUMBER(scale);
   ss << '\n';

   const int rows = matrix.rows();
   const int cols = matrix.cols();
   for (int i = 1; i <= rows; ++i)
      for (int k = 1; k <= cols; ++k) {
         ss << boost::format(mixing_matrix_formatter) % i % k % matrix(i-1,k-1)
            % (symbol + "(" + std::to_string(i) + "," + std::to_string(k) + ")");
      }

   set_block(ss);
}

} // namespace flexiblesusy

#endif
