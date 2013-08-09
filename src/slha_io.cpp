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

#include "slha_io.hpp"
#include "logger.hpp"
#include "lowe.h"
#include "config.h"

#include <fstream>

namespace flexiblesusy {

SLHA_io::SLHA_io()
   : input_filename()
   , data()
{
}

/**
 * Constructor
 * @brief opens SLHA input file and reads the content
 * @param input_filename_ SLHA input file name
 */
SLHA_io::SLHA_io(const std::string& input_filename_)
   : input_filename(input_filename_)
{
   std::ifstream ifs(input_filename_);
   data.read(ifs);
}

/**
 * @brief opens SLHA input file and reads the content
 * @param f SLHA input file name
 */
void SLHA_io::set_input_file(const std::string& f)
{
   input_filename = f;
   std::ifstream ifs(f);
   data.read(ifs);
}

void SLHA_io::fill(QedQcd& oneset)
{
   using namespace std::placeholders;
   SLHA_io::Tuple_processor sminputs_processor
      = std::bind(&SLHA_io::process_sminputs_tuple, oneset, _1, _2);

   read_block("SMINPUTS", sminputs_processor);
}

void SLHA_io::read_block(const std::string& block_name, Tuple_processor processor)
{
   if (data.find(block_name) == data.cend()) {
      WARNING("Block " << block_name << " not found in SLHA2 file "
              << input_filename);
      return;
   }

   for (SLHAea::Block::const_iterator line = data.at(block_name).cbegin(),
        end = data.at(block_name).cend(); line != end; ++line) {
      if (!line->is_data_line())
         continue;

      if (line->size() >= 2) {
         const int key = SLHAea::to<int>((*line)[0]);
         const double value = SLHAea::to<double>((*line)[1]);
         processor(key, value);
      } else {
         WARNING(block_name << " entry has not enough columns");
      }
   }
}

void SLHA_io::set_spinfo()
{
   data["SPINFO"][0][1] = PKGNAME;
   data["SPINFO"][0][2] = "# spectrum calculator";
   data["SPINFO"][1][1] = VERSION;
   data["SPINFO"][1][2] = "# version number of " PKGNAME;
}

void SLHA_io::write_to_file(const std::string& file_name)
{
   set_spinfo();

   std::ofstream ofs(file_name);
   ofs << data;
}

/**
 * fill oneset from given key - value pair
 *
 * @param oneset low-energy data set
 * @param key SLHA key in SMINPUTS
 * @param value value corresponding to key
 */
void SLHA_io::process_sminputs_tuple(QedQcd& oneset, int key, double value)
{
   switch (key) {
   case 1:
      oneset.setAlpha(ALPHA, 1.0 / value);
      break;
   case 2:
      // Gmu cannot be set yet
      break;
   case 3:
      oneset.setAlpha(ALPHAS, value);
      break;
   case 4:
      // MZ cannot be set yet
      // oneset.setMu(value);
      break;
   case 5:
      oneset.setMass(mBottom, value);
      oneset.setMbMb(value);
      break;
   case 6:
      oneset.setPoleMt(value);
      break;
   case 7:
      oneset.setMass(mTau, value);
      oneset.setPoleMtau(value);
      break;
   case 8:
      break;
   case 11:
      oneset.setMass(mElectron, value);
      break;
   case 12:
      break;
   case 13:
      oneset.setMass(mMuon, value);
      break;
   case 14:
      break;
   case 21:
      oneset.setMass(mDown, value);
      break;
   case 22:
      oneset.setMass(mUp, value);
      break;
   case 23:
      oneset.setMass(mStrange, value);
      break;
   case 24:
      oneset.setMass(mCharm, value);
      break;
   default:
      WARNING("Unrecognized key in SMINPUTS: " << key);
      break;
   }
}

} // namespace flexiblesusy
