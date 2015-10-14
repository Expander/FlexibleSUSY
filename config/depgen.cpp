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

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

/// returns directory from file name
std::string directory(const std::string& file_name)
{
   const std::size_t pos = file_name.find_last_of("/\\");
   if (pos == std::string::npos)
      return ".";
   return file_name.substr(0,pos);
}

/// returns file name w/o directory
std::string filename(const std::string& file_name)
{
   const std::size_t pos = file_name.find_last_of("/\\");
   return file_name.substr(pos+1);
}

/// checks if given file exists
bool file_exists(const std::string& name)
{
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

/// deletes duplicate elements from a vector
std::vector<std::string> delete_duplicates(const std::vector<std::string>& vec)
{
   std::vector<std::string> unique_vector(vec);

   std::sort(unique_vector.begin(), unique_vector.end());
   unique_vector.erase(std::unique(unique_vector.begin(), unique_vector.end()),
                       unique_vector.end());

   return unique_vector;
}

/// replace file name extension by `ext'
std::string replace_extension_by(const std::string& str, const std::string& ext)
{
   const std::string no_ext(str.substr(0, str.find_last_of(".")));
   return no_ext + '.' + ext;
}

/// tests whether `str' starts with `prefix'
bool starts_with(const std::string& str, const std::string& prefix)
{
   return !str.compare(0, prefix.size(), prefix);
}

/// removes whitespace from left side of string
std::string trim_left(const std::string& s)
{
   std::string str(s);

   str.erase(str.begin(),
             std::find_if(str.begin(), str.end(),
                          std::not1(std::ptr_fun<int, int>(std::isspace))));

   return str;
}

/// print usage message
void print_usage(const std::string& program_name)
{
   std::cout << "Usage: " << program_name << " [options] filename\n"
      "\n"
      "Options:\n"
      "  -I<path>      Search for header files in <path>\n"
      "  -MF <file>    Write dependencies to <file>\n"
      "  -MM           Ignore system headers enclosed by < and >\n"
      "  -MMD <file>   Equivalent to -MM -MF <file>\n"
      "  -MT <target>  Set name of the target\n"
      "  -o <file>     Equivalent to -MF <file>\n"
      "  --help,-h     Print this help message and exit\n"
      "\n"
      "Unsupported options:\n"
      "  -M            Add system headers to dependency list\n"
      "  -MD <file>    Equivalent to -M -MF <file>\n"
      "  -MG           Add missing headers to dependency list\n"
      "  -MP           Add phony target for each dependency other than main file\n"
      "  -MQ <target>  Same as -MT <traget> but quote characters special to make\n";
}

/// print dependency list
void print_dependencies(const std::string& target_name,
                        const std::vector<std::string>& dependencies,
                        std::ostream& ostr = std::cout)
{
   ostr << target_name << ':';

   for (std::vector<std::string>::const_iterator it = dependencies.begin(),
           end = dependencies.end(); it != end; ++it) {
      ostr << ' ' << *it;
   }

   ostr << '\n';
}

/// extract include statements from file (ignoring system headers)
std::vector<std::string> get_includes(const std::string& file_name)
{
   std::ifstream istr(file_name.c_str());
   std::vector<std::string> includes;
   std::string line;

   while (std::getline(istr, line)) {
      const std::string tline(trim_left(line));
      if (!tline.empty() && tline[0] == '#') {
         const std::string ttline(trim_left(tline.substr(1)));
         if (starts_with(ttline, "include")) {
            const std::string header(trim_left(ttline.substr(7)));
            if (!header.empty() && header[0] == '"')
               includes.push_back(header.substr(1, header.size()-2));
         }
      }
   }

   return includes;
}

/// returns files existing in `dir'
std::vector<std::string> get_existing(const std::string& dir,
                                      const std::vector<std::string>& files)
{
   std::vector<std::string> existing_files;

   for (std::vector<std::string>::const_iterator it = files.begin(),
           end = files.end(); it != end; ++it) {
      const std::string file_with_path(dir == "." ? *it : dir + '/' + *it);
      if (file_exists(file_with_path))
         existing_files.push_back(file_with_path);
   }

   return existing_files;
}

/// search recursively for include statments in `file_name'
/// taking into account only directories given in `paths'
std::vector<std::string> search_includes(const std::string& file_name,
                                         const std::vector<std::string>& paths,
                                         unsigned max_depth = 10)
{
   if (max_depth == 0)
      return std::vector<std::string>();

   // find included files from #include statements
   std::vector<std::string> includes(get_includes(file_name));

   // select only files that exist in paths
   std::vector<std::string> existing;
   for (std::vector<std::string>::const_iterator it = paths.begin(),
           end = paths.end(); it != end; ++it) {
      const std::vector<std::string> existing_in_path(get_existing(*it, includes));
      existing.insert(existing.end(), existing_in_path.begin(), existing_in_path.end());
   }

   // search recursively for included files in existing headers
   const std::vector<std::string> tmp_existing(existing);
   for (std::vector<std::string>::const_iterator it = tmp_existing.begin(),
           end = tmp_existing.end(); it != end; ++it) {
      const std::vector<std::string> sub_existing(search_includes(*it, paths, max_depth - 1));
      existing.insert(existing.end(), sub_existing.begin(), sub_existing.end());
   }

   return existing;
}

int main(int argc, const char* argv[])
{
   // include paths
   std::vector<std::string> paths;
   std::string file_name, target_name, output_file;

   for (int i = 1; i < argc; i++) {
      const std::string arg(argv[i]);
      if (starts_with(arg, "-I") && arg.length() > 2) {
         paths.push_back(arg.substr(2));
         continue;
      }
      if (arg == "-MM") {
         // ignore headers from system directories (default)
         continue;
      }
      if (arg == "-M" || arg == "-MG" || arg == "-MP") {
         std::cerr << "Warning: ignoring unsupported option " << arg << '\n';
         continue;
      }
      if (arg == "-MD" || arg == "-MQ") {
         std::cerr << "Warning: ignoring unsupported option " << arg << '\n';
         i++;
         continue;
      }
      if (arg == "-MT" && i + 1 < argc) {
         target_name = argv[++i];
         continue;
      }
      if ((arg == "-MF" || arg == "-MMD" || arg == "-o") && i + 1 < argc) {
         output_file = argv[++i];
         continue;
      }
      if (arg == "--help" || arg == "-h") {
         print_usage(argv[0]);
         return EXIT_SUCCESS;
      }
      // interpret last argument as file name
      if (i + 1 == argc) {
         file_name = arg;
         if (!file_exists(file_name)) {
            std::cerr << "Error: file does not exist: " << file_name << '\n';
            return EXIT_FAILURE;
         }
         continue;
      }

      std::cerr << "Error: unknown option: " << arg << '\n';
      print_usage(argv[0]);
      return EXIT_FAILURE;
   }

   // search for header inclusions
   std::vector<std::string> dependencies
      = delete_duplicates(search_includes(file_name, paths));

   // add file name to dependency list
   dependencies.insert(dependencies.begin(), file_name);

   if (target_name.empty())
      target_name = replace_extension_by(file_name, "o");

   if (output_file.empty()) {
      print_dependencies(target_name, dependencies);
   } else {
      std::ofstream ostr(output_file.c_str());
      print_dependencies(target_name, dependencies, ostr);
   }

   return EXIT_SUCCESS;
}
