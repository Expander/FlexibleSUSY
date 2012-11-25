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

#include "running_coupling.hpp"

#include <fstream>

Running_coupling::Running_coupling()
   : couplings(TData())
{
}

Running_coupling::~Running_coupling()
{
}

/**
 * Get the couplings at the largest scale
 *
 * @return a pair with the scale and a DoubleVector which contains the
 * couplings at this scale
 */
Running_coupling::TTouple Running_coupling::get_max_scale() const
{
   if (couplings.empty()) {
      ERROR("Data container is empty!");
      return TTouple(0.0, DoubleVector(1));
   }

   // find gauge couplings at the greatest scale
   TData::const_iterator maxScale
      = max_element(couplings.begin(), couplings.end(), TDataComp());

   return *maxScale;
}

/**
 * Delete all internal couplings.
 */
void Running_coupling::reset()
{
   couplings.clear();
}

/**
 * write help line which describes the written data
 *
 * @param coupling_name coupling name
 * @param fout output stream
 * @param number_of_couplings number of couplings at each scale
 * @param start_index index of first coupling (usually set to 1)
 */
void Running_coupling::write_comment_line(char coupling_name, std::ofstream& fout,
                                          std::size_t number_of_couplings,
                                          int start_index) const
{
   if (!fout.good())
      return;

   fout << "# scale [GeV]";

   for (std::size_t i = 0; i < number_of_couplings; ++i)
      fout << " | " << coupling_name << "(" << start_index + i << ")";

   fout << std::endl;
}

/**
 * Write all couplings to a text file.
 *
 * @param file_name name of file to write the data to
 */
void Running_coupling::write_to_file(const std::string& file_name) const
{
   const char coupling_name = 'g';

   if (couplings.empty())
      return;

   std::ofstream filestr(file_name.c_str(), ios::out);
   VERBOSE_MSG("opening file: " << file_name.c_str());
   if (filestr.fail()) {
      ERROR("can't open file " << file_name
            << " for writing running couplings");
      return;
   }

   write_comment_line(coupling_name, filestr,
                      couplings.front().second.size(),
                      couplings.front().second.displayStart());

   // write data
   for (TData::const_iterator it = couplings.begin();
        it != couplings.end(); ++it) {
      if (!filestr.good()) {
         ERROR("file " << file_name << " is corrupted");
         break;
      }

      filestr.width(16);
      filestr << left << it->first;

      // write all gauge couplings in order
      for (int i = it->second.displayStart();
           i <= it->second.displayEnd(); ++i) {
         filestr.width(16);
         filestr << left << it->second(i);
      }

      filestr << endl;
   }

   filestr.close();
   VERBOSE_MSG("file written: " << file_name.c_str());
}
