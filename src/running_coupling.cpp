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
#include "logger.hpp"

#include <fstream>

using namespace std;

Running_coupling::Running_coupling()
   : fGaugeCouplings(TData())
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
   const TData* data = &fGaugeCouplings;

   if (data->empty())
      ERROR("Data container is empty!");

   // find gauge couplings at the greatest scale
   TData::const_iterator maxScale
      = max_element(data->begin(), data->end(), TDataComp());

   return *maxScale;
}

/**
 * Delete all internal couplings.
 */
void Running_coupling::reset()
{
   fGaugeCouplings.clear();
}

/**
 * write help line which describes the written data
 *
 * @param couplingName coupling name
 * @param fout output stream
 * @param nElements number of couplings at each scale
 * @param startIndex index of first coupling (usually set to 1)
 */
void Running_coupling::write_comment_line(char couplingName, ofstream& fout,
                                        unsigned int nElements, int startIndex) const
{
   if (!fout.good())
      return;

   fout << "# scale [GeV]";

   for (unsigned int i = 0; i < nElements; ++i)
      fout << " | " << couplingName << "(" << startIndex + i << ")";

   fout << endl;
}

/**
 * Write all couplings to a text file.
 *
 * @param fileName name of file to write the data to
 */
void Running_coupling::write_to_file(const string& fileName) const
{
   const TData* data = &fGaugeCouplings;
   char couplingName = 'g';

   if (data->empty())
      return;

   ofstream filestr(fileName.c_str(), ios::out);
   if (filestr.fail()) {
      ERROR("can't open file " << fileName
            << " for writing running coupling " << couplingName);
      return;
   }

   write_comment_line(couplingName, filestr,
                      data->front().second.size(),
                      data->front().second.displayStart());

   // write data
   for (TData::const_iterator it = data->begin();
        it != data->end(); ++it) {
      if (!filestr.good()) {
         ERROR("file " << fileName << " is corrupted");
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
   VERBOSE_MSG("file written: " << fileName.c_str());
}
