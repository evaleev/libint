/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_codegenjobs_h_
#define _libint2_src_bin_libint_codegenjobs_h_

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iosfwd>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace libint2 {

class CodeContext;
class CompilationParameters;
class Libint2Iface;
class RecurrenceRelation;

/// The effects of generating the code for one shell set on the state of the
/// generator (the library interface and the parameters of its task), recorded
/// so that they can be applied in the same order however the shell sets were
/// generated (see CodegenJobs).
struct ShellSetJournal {
  /// lines for Libint2Iface::to_static_init()
  std::vector<std::string> static_init;
  /// lines for Libint2Iface::to_int_iface()
  std::vector<std::string> int_iface;
  /// arguments of TaskParameters::max_stack_size(am, size)
  std::vector<std::pair<unsigned int, unsigned int> > max_stack_size;
  /// arguments of TaskParameters::max_ntarget(ntarget)
  std::vector<unsigned int> max_ntarget;
};

/**
   CodegenJobs generates the code of the library: the builders (see
   build_libint.cc) add a job for each shell set of each task, and run()
   executes the jobs, then generates the code for the set-level recurrence
   relations (RRs) they use.

   With LIBINT_GENERATOR_JOBS=N > 1 (see generator_jobs()) the jobs are executed
   by N forked worker processes, which take them, largest first, from a counter
   in shared memory, and write the code into their own subdirectories of the
   source directory. Each worker sends back the journals of its jobs, the
   external symbols of the tasks, and the RRs that the tasks need, which this
   process combines; it then assigns each RR to one of the workers that has it,
   and that worker generates its code. Finally the files are moved into the
   source directory. This requires that the code generated for a shell set or an
   RR does not depend on what was generated before it in the same process (see
   src/bin/check_codegen_order.cmake).

   The journals are applied in the order in which the jobs were added, hence
   the library interface does not depend on the number of processes either.
   File names that differ only in case get suffixes (see
   CodeContext::label_to_filename()); these are assigned at the end in the
   lexicographic order of the names, rather than in the order in which they
   were seen, so that the file names do not depend on the number of processes.
*/
class CodegenJobs {
 public:
  using Job = std::function<void(ShellSetJournal&)>;

  /// the instance used by the builders
  static CodegenJobs& instance();

  /// adds the job that generates the code for a shell set of task @p task;
  /// @p cost estimates its cost (only the relative order of the costs matters)
  void add(const std::string& task, double cost, Job job);

  /// executes the jobs (in @p njobs processes) and generates the code for the
  /// set-level RRs; applies the journals to @p iface and to the tasks
  /// @p on_worker_start(k) is called in worker k when it starts (e.g. to
  /// adjust the progress reporting)
  void run(std::ostream& os,
           const std::shared_ptr<CompilationParameters>& cparams,
           Libint2Iface& iface, std::size_t njobs,
           const std::function<void(std::size_t)>& on_worker_start = {});

 private:
  struct Entry {
    std::string task;
    double cost;
    Job job;
  };
  std::vector<Entry> jobs_;

  void run_serially(std::ostream& os,
                    const std::shared_ptr<CompilationParameters>& cparams,
                    Libint2Iface& iface);
  void run_in_parallel(std::ostream& os,
                       const std::shared_ptr<CompilationParameters>& cparams,
                       Libint2Iface& iface, std::size_t njobs,
                       const std::function<void(std::size_t)>& on_worker_start);
};

/// @return the number of processes that generate the code, given by the
/// environment variable LIBINT_GENERATOR_JOBS (0 = the number of hardware
/// threads); 1 if it is not set
std::size_t generator_jobs();

/// @return the set-level RRs that the tasks need, in the order in which they
/// were created (or as LIBINT_RR_ORDER requests, see permute_for_testing())
std::vector<std::shared_ptr<RecurrenceRelation> > rrs_to_generate();

/// generates the code for @p rr into @p prefix + its file name, then removes
/// @p rr from the RRStack and purges the singleton stacks
void generate_rr(std::ostream& os, const std::shared_ptr<CodeContext>& context,
                 const std::string& prefix,
                 const std::shared_ptr<RecurrenceRelation>& rr);

/// generates the code for rrs_to_generate() in this process
void generate_rr_code(std::ostream& os,
                      const std::shared_ptr<CompilationParameters>& cparams);

/// The generated code must not depend on the order in which the shell sets
/// and the RRs are generated. To test this, the value of the environment
/// variable @p envvar (for testing only, see
/// src/bin/check_codegen_order.cmake), if set and not empty, permutes @p v:
/// - `reverse`: reverse order
/// - `label`: in the order of @p label(element)
/// - `shuffle[:<seed>]`: random permutation (seed 0 by default) that is the
///   same on every platform (std::shuffle is not: the standard does not
///   specify how it uses the random number engine)
template <typename T, typename Label>
void permute_for_testing(std::vector<T>& v, const char* envvar, Label&& label) {
  const char* order = std::getenv(envvar);
  if (order == nullptr || *order == '\0') return;
  const std::string o(order);
  const auto invalid = [&]() {
    return std::invalid_argument(std::string(envvar) + "=" + o +
                                 ": expected reverse, label, or "
                                 "shuffle[:<seed>] with an unsigned "
                                 "integer <seed>");
  };
  if (o == "reverse")
    std::reverse(v.begin(), v.end());
  else if (o == "label")
    std::stable_sort(v.begin(), v.end(), [&](const T& a, const T& b) {
      return label(a) < label(b);
    });
  else if (o == "shuffle" || o.compare(0, 8, "shuffle:") == 0) {
    unsigned long seed = 0;
    if (o.size() > 7) {  // shuffle:<seed>
      const std::string s = o.substr(8);
      if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw invalid();
      try {
        seed = std::stoul(s);
      } catch (const std::out_of_range&) {
        throw invalid();
      }
    }
    // Fisher-Yates on the raw output of mt19937, which the standard specifies
    std::mt19937 rng(static_cast<std::mt19937::result_type>(seed));
    for (auto i = v.size(); i > 1; --i) std::swap(v[i - 1], v[rng() % i]);
  } else
    throw invalid();
  std::cout << envvar << "=" << o
            << ": generating in a non-default order (for testing only)"
            << std::endl;
}

/// the registry of CodeContext::label_to_filename(): case-folded name ->
/// the distinct names with that folded name, in the order they were seen
std::map<std::string, std::vector<std::string> >& filename_registry();
/// @return @p name with all letters in lower case
std::string fold_case(const std::string& name);
/// @return the file name of the @p idx -th distinct name with the same
/// case-folded name
std::string disambiguated_filename(const std::string& name, std::size_t idx);

}  // namespace libint2

#endif  // header guard
