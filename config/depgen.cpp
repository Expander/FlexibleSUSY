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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace flexiblesusy {
namespace depgen {

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

/// returns file names w/o directory
std::vector<std::string> filenames(const std::vector<std::string>& file_names)
{
   std::vector<std::string> result(file_names.size());
   std::transform(file_names.begin(), file_names.end(), result.begin(), filename);
   return result;
}

/// checks if given file exists
bool file_exists(const std::string& name)
{
   if (FILE *file = std::fopen(name.c_str(), "r")) {
      std::fclose(file);
      return true;
   }
   return false;
}

struct Is_not_duplicate {
   bool operator()(const std::string& element) {
      return s.insert(element).second;
   }
private:
   std::set<std::string> s;
};

struct Is_not_duplicate_ignore_path {
   bool operator()(const std::string& element) {
      return s.insert(filename(element)).second;
   }
private:
   std::set<std::string> s;
};

/// deletes duplicate elements from a vector (preseves order)
template <typename Predicate = decltype(Is_not_duplicate())>
std::vector<std::string> delete_duplicates(
   const std::vector<std::string>& vec,
   Predicate pred = Is_not_duplicate())
{
   std::vector<std::string> unique_vector;

   std::copy_if(vec.begin(), vec.end(), std::back_inserter(unique_vector),
                std::ref(pred));

   return unique_vector;
}

/// replace file name extension by `ext'
std::string replace_extension(const std::string& str, const std::string& ext)
{
   const std::string no_ext(str.substr(0, str.find_last_of('.')));
   return no_ext + '.' + ext;
}

/// tests whether `str' starts with `prefix'
bool starts_with(const std::string& str, const std::string& prefix)
{
   return str.compare(0, prefix.size(), prefix) == 0;
}

/// removes whitespace from left side of string
void trim_left(std::string& str)
{
   str.erase(str.begin(),
             std::find_if(str.begin(), str.end(),
                          [] (std::string::value_type c) { return std::isspace(c) == 0; }));
}

/// returns copy of s with whitespace removed from left side of string
std::string trim_left_copy(const std::string& s)
{
   std::string str(s);
   trim_left(str);
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
      "  -MG           Add missing headers to dependency list\n"
      "  -MI           Ignore errors of non-existing header(s)\n"
      "  -MM           Ignore system headers enclosed by < and >\n"
      "  -MMD <file>   Equivalent to -MM -MF <file>\n"
      "  -MP           Add phony target for each dependency other than main file\n"
      "  -MT <target>  Set name of the target\n"
      "  -o <file>     Equivalent to -MF <file>\n"
      "  --help,-h     Print this help message and exit\n"
      "\n"
      "Unsupported options:\n"
      "  -M            Add system headers to dependency list\n"
      "  -MD <file>    Equivalent to -M -MF <file>\n"
      "  -MQ <target>  Same as -MT <traget> but quote characters special to make\n";
}

/// print dependency list
void print_dependencies(std::ostream& ostr,
                        const std::string& target_name,
                        const std::vector<std::string>& dependencies)
{
   ostr << target_name << ':';

   for (const auto& d: dependencies)
      ostr << ' ' << d;

   ostr << '\n';
}

/// print empty phony targets for each dependency
void print_empty_phony_targets(std::ostream& ostr, const std::vector<std::string>& dependencies)
{
   for (const auto& d: dependencies)
      ostr << '\n' << d << ":\n";
}

/// returns file name from include "..." statement
std::string get_filename_from_include(std::string line)
{
   trim_left(line);

   if (line.empty() || line[0] != '#')
      return "";

   // skip `#' and following whitespace
   line = trim_left_copy(line.substr(std::strlen("#")));

   if (!starts_with(line, "include"))
      return "";

   // skip `include'
   line = trim_left_copy(line.substr(std::strlen("include")));

   // extract file name from "file-name"
   std::size_t pos1 = line.find_first_of('"');
   if (pos1 == std::string::npos)
      return "";

   pos1++;

   std::size_t pos2 = line.find_first_of('"', pos1);
   if (pos2 == std::string::npos)
      return "";

   pos2--;

   return line.substr(pos1, pos2);
}

/// extract include statements from file (ignoring system headers)
std::vector<std::string> get_included_files(const std::string& file_name)
{
   std::ifstream istr(file_name);
   std::vector<std::string> includes;
   std::string line;

   while (std::getline(istr, line)) {
      auto file = get_filename_from_include(line);
      if (!file.empty())
         includes.push_back(std::move(file));
   }

   return includes;
}

/// prepend `str' to all strings
std::vector<std::string> prepend(const std::string& str,
                                 const std::vector<std::string>& strings)
{
   std::vector<std::string> result;
   result.reserve(strings.size());

   for (auto& s: strings) {
      result.push_back(str + s);
   }

   return result;
}

/// insert `str' at the beginning of vector
std::vector<std::string> insert_at_front(
   const std::vector<std::string>& strings,
   const std::string& str)
{
   auto result = strings;
   result.insert(result.begin(), str);

   return result;
}

/// returns elements of `vec' for which pred(f) == true
template <class Predicate>
std::vector<std::string> filter(
   const std::vector<std::string>& vec,
   const Predicate& pred)
{
   std::vector<std::string> match;

   std::copy_if(vec.begin(), vec.end(),
                std::back_inserter(match), pred);

   return match;
}

/// returns files in directory `dir' for which pred(f) == true
template <class Predicate>
std::vector<std::string> filter_files(
   const std::string& dir,
   const std::vector<std::string>& files,
   const Predicate& pred)
{
   const std::string dirname(dir.empty() || dir == "." ? "" : dir + '/');
   const auto files_in_dir = prepend(dirname, files);

   return filter(files_in_dir, pred);
}

/// returns all elements of `v1', which are not `v2'
std::vector<std::string> complement(
   const std::vector<std::string>& v1,
   const std::vector<std::string>& v2)
{
   auto tv1 = v1;
   auto tv2 = v2;

   std::sort(tv1.begin(), tv1.end());
   std::sort(tv2.begin(), tv2.end());

   std::vector<std::string> diff;

   std::set_difference(tv1.begin(), tv1.end(),
                       tv2.begin(), tv2.end(),
                       std::back_inserter(diff));

   return diff;
}

/// concatenate strings with separator
template <typename T>
std::string concat(const std::vector<std::string>& strings, const T& separator)
{
   std::string result;

   for (const auto& s: strings)
      result += s + separator;

   return result;
}

/// search recursively for include statments in `file_name'
/// taking into account only directories given in `paths'
void search_includes(const std::string& file_name,
                     const std::vector<std::string>& paths,
                     std::vector<std::string>& result,
                     bool include_non_existing,
                     bool ignore_non_existing,
                     int max_depth)
{
   if (max_depth <= 0) {
      throw std::runtime_error(
         "Error: #include nested too deeply (maximum depth: "
         + std::to_string(max_depth) + "): " + file_name);
   }

   // find included files from #include statements, that are not
   // already in result
   const auto includes = complement(get_included_files(file_name), filenames(result));

   // select only files that exist in paths
   std::vector<std::string> existing;
   for (const auto& p: paths) {
      const auto existing_in_path = filter_files(p, includes, file_exists);
      existing.insert(existing.end(), existing_in_path.cbegin(), existing_in_path.cend());
      result.insert(result.end(), existing_in_path.cbegin(), existing_in_path.cend());
   }

   // search recursively for included files in existing headers
   for (const auto& f: existing)
      search_includes(f, paths, result, include_non_existing,
                      ignore_non_existing, max_depth - 1);

   // search for non-existing headers
   const auto non_existing = complement(filenames(includes), filenames(existing));

   if (!ignore_non_existing && !include_non_existing && !non_existing.empty()) {
      throw std::runtime_error(
         "Error: cannot find the following header file(s): "
         + concat(non_existing, ' '));
   }

   if (include_non_existing)
      result.insert(result.end(), non_existing.cbegin(), non_existing.cend());
}

/// search recursively for include statments in `file_name'
/// taking into account only directories given in `paths'
std::vector<std::string> search_includes(const std::string& file_name,
                                         const std::vector<std::string>& paths,
                                         bool include_non_existing,
                                         bool ignore_non_existing,
                                         int max_depth = 100)
{
   std::vector<std::string> result;
   search_includes(file_name, paths, result, include_non_existing,
                   ignore_non_existing, max_depth);
   return result;
}

} // namespace depgen
} // namespace flexiblesusy

int main(int argc, char* argv[])
{
   using namespace flexiblesusy::depgen;

   if (argc < 2) {
      std::cerr << "Error: no file given\n";
      print_usage(argv[0]);
      return EXIT_FAILURE;
   }

   // include paths
   std::vector<std::string> paths;
   std::string file_name, target_name, output_file;
   bool include_non_existing = false; // -MG
   bool ignore_non_existing = false; // -MI
   bool add_empty_phony_targets = false; // -MP

   for (int i = 1; i < argc; i++) {
      const std::string arg(argv[i]);
      if (starts_with(arg, "-D")) {
         continue;
      }
      if (starts_with(arg, "-I") && arg.length() > 2) {
         paths.push_back(arg.substr(std::strlen("-I")));
         continue;
      }
      if (arg == "-MG") {
         include_non_existing = true;
         continue;
      }
      if (arg == "-MI") {
         ignore_non_existing = true;
         continue;
      }
      if (arg == "-MM") {
         // ignore headers from system directories (default)
         continue;
      }
      if (arg == "-M") {
         std::cerr << "Warning: ignoring unsupported option " << arg << '\n';
         continue;
      }
      if (arg == "-MD" || arg == "-MQ") {
         std::cerr << "Warning: ignoring unsupported option " << arg << '\n';
         i++;
         continue;
      }
      if (arg == "-MP") {
         add_empty_phony_targets = true;
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

   // select output stream
   std::ofstream fstr(output_file);
   std::ostream* ostr = &std::cout;

   if (!output_file.empty()) {
      if (!fstr.good()) {
         std::cerr << "Error: cannot write to file " << output_file << '\n';
         return EXIT_FAILURE;
      }
      ostr = &fstr;
   }

   // include paths
   paths = insert_at_front(paths, directory(file_name));
   paths.emplace_back(".");
   paths = delete_duplicates(paths);

   try {
      // search for header inclusions in file
      const auto dependencies
         = delete_duplicates(
            search_includes(file_name, paths, include_non_existing,
                            ignore_non_existing),
            Is_not_duplicate_ignore_path());

      if (target_name.empty())
         target_name = replace_extension(filename(file_name), "o");

      // output
      print_dependencies(*ostr, target_name, insert_at_front(dependencies, file_name));

      if (add_empty_phony_targets)
         print_empty_phony_targets(*ostr, dependencies);
   } catch (const std::exception& e) {
      std::cerr << e.what() << '\n';
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
