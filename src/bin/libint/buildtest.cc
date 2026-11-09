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

#include <context.h>
#include <dims.h>
#include <rr.h>
#include <task.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>
#if defined(__unix__) || defined(__APPLE__)
#include <sys/mman.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

using namespace std;
using namespace libint2;

namespace {

/// The generated code must not depend on the order in which the RRs are
/// generated. To test this, the environment variable LIBINT_RR_ORDER (for
/// testing only, see src/bin/check_rr_order.cmake) reorders @p rrs:
/// - `reverse`: reverse creation order
/// - `label`: lexicographic order of the labels
/// - `shuffle[:<seed>]`: random permutation (seed 0 by default) that is the
///   same on every platform (std::shuffle is not: the standard does not
///   specify how it uses the random number engine)
/// Unset or empty keeps the creation order.
void reorder_rrs(std::vector<std::shared_ptr<RecurrenceRelation> >& rrs,
                 std::ostream& os) {
  const char* order = std::getenv("LIBINT_RR_ORDER");
  if (order == nullptr || *order == '\0') return;
  const std::string o(order);
  const auto invalid = [&o]() {
    return std::invalid_argument("LIBINT_RR_ORDER=" + o +
                                 ": expected reverse, label, or "
                                 "shuffle[:<seed>] with an unsigned "
                                 "integer <seed>");
  };
  if (o == "reverse")
    std::reverse(rrs.begin(), rrs.end());
  else if (o == "label")
    std::sort(rrs.begin(), rrs.end(), [](const auto& a, const auto& b) {
      return a->label() < b->label();
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
    for (auto i = rrs.size(); i > 1; --i) std::swap(rrs[i - 1], rrs[rng() % i]);
  } else
    throw invalid();
  os << "LIBINT_RR_ORDER=" << o
     << ": generating the RRs in a non-default order (for testing only)"
     << endl;
}

}  // namespace

namespace libint2 {

#if defined(__unix__) || defined(__APPLE__)
#define LIBINT_GENERATOR_CAN_FORK 1
#else
#define LIBINT_GENERATOR_CAN_FORK 0
#endif

namespace {

#if LIBINT_GENERATOR_CAN_FORK
/// @return the peak resident memory of this process, in bytes
std::size_t peak_rss_bytes() {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0;
#if defined(__APPLE__)
  return static_cast<std::size_t>(usage.ru_maxrss);  // bytes
#else
  return static_cast<std::size_t>(usage.ru_maxrss) * 1024;  // kilobytes
#endif
}

/// @return the physical memory, in bytes (0 if unknown)
std::size_t physical_memory_bytes() {
  const long pages = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGESIZE);
  if (pages <= 0 || page_size <= 0) return 0;
  return static_cast<std::size_t>(pages) * static_cast<std::size_t>(page_size);
}
#endif

/// @return the number of processes that generate the set-level RRs, given by
/// the environment variable LIBINT_GENERATOR_JOBS; 1 if it is not set. 0
/// selects the number of hardware threads, limited so that the processes
/// fit in 3/4 of the physical memory: each process eventually needs about as
/// much memory as the generator (copy-on-write sharing does not last, since
/// generating code updates reference counts throughout the shared data)
std::size_t generator_jobs(std::ostream& os) {
  const char* env = std::getenv("LIBINT_GENERATOR_JOBS");
  if (env == nullptr || *env == '\0') return 1;
  const std::string value(env);
  long njobs = -1;
  try {
    std::size_t pos = 0;
    njobs = std::stol(value, &pos);
    if (pos != value.size()) njobs = -1;
  } catch (...) {
    njobs = -1;
  }
  if (njobs < 0)
    throw std::invalid_argument("LIBINT_GENERATOR_JOBS=" + value +
                                ": expected a non-negative integer");
  const bool automatic = (njobs == 0);
  if (automatic) njobs = std::max(1u, std::thread::hardware_concurrency());
#if LIBINT_GENERATOR_CAN_FORK
  const auto rss = peak_rss_bytes();
  const auto mem = physical_memory_bytes();
  if (rss > 0 && mem > 0) {
    const long fit = std::max<long>(1, static_cast<long>(mem / 4 * 3 / rss));
    if (njobs > fit) {
      if (automatic) {
        os << "generate_rr_code: limiting the number of processes to " << fit
           << " (each needs about " << (rss >> 20) << " MB of the "
           << (mem >> 20) << " MB of physical memory)" << std::endl;
        njobs = fit;
      } else
        os << "generate_rr_code: WARNING: LIBINT_GENERATOR_JOBS=" << njobs
           << " processes may need about " << ((rss >> 20) * njobs)
           << " MB, more than 3/4 of the " << (mem >> 20)
           << " MB of physical memory" << std::endl;
    }
  }
#endif
  return static_cast<std::size_t>(njobs);
}

/// generates the RRs in @p njobs forked processes that take the RRs, largest
/// first, from a shared counter. The output files of the RRs are written by
/// the workers; the only state that generating an RR adds to that of the
/// generator are the external symbols of the tasks that use it (see
/// RecurrenceRelation::generate_code), which the workers send back to this
/// process. Requires that the generated code does not depend on the order in
/// which the RRs are generated (see src/bin/check_rr_order.cmake).
template <typename Generate>
void generate_in_parallel(
    std::ostream& os,
    const std::vector<std::shared_ptr<RecurrenceRelation> >& rrs,
    std::size_t njobs, Generate&& generate) {
#if LIBINT_GENERATOR_CAN_FORK
  // largest targets first, so that the workers finish at about the same time
  std::vector<std::size_t> schedule(rrs.size());
  std::iota(schedule.begin(), schedule.end(), 0);
  std::stable_sort(
      schedule.begin(), schedule.end(), [&](std::size_t a, std::size_t b) {
        return rrs[a]->rr_target()->size() > rrs[b]->rr_target()->size();
      });

  // index of the next RR in schedule, shared by the workers
  using counter_type = std::atomic<std::size_t>;
  static_assert(counter_type::is_always_lock_free,
                "the RR counter shared by the generator processes must be "
                "lock-free");
  void* shm = mmap(nullptr, sizeof(counter_type), PROT_READ | PROT_WRITE,
                   MAP_SHARED | MAP_ANON, -1, 0);
  if (shm == MAP_FAILED)
    throw std::runtime_error(std::string("generate_rr_code: mmap failed: ") +
                             std::strerror(errno));
  counter_type* next = new (shm) counter_type(0);

  const std::string symfile_prefix =
      ".libint_rr_symbols." + std::to_string(getpid()) + ".";
  auto symfile = [&](std::size_t k) {
    return symfile_prefix + std::to_string(k);
  };

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  os << "generating the code for " << rrs.size() << " RRs in " << njobs
     << " processes" << std::endl;
  os.flush();
  std::cout.flush();
  std::cerr.flush();
  std::vector<pid_t> workers;
  bool fork_failed = false;
  int fork_errno = 0;
  for (std::size_t k = 0; k != njobs; ++k) {
    const pid_t pid = fork();
    if (pid < 0) {
      fork_failed = true;
      fork_errno = errno;
      break;
    }
    if (pid == 0) {  // worker
      int status = 0;
      try {
        for (std::size_t i; (i = next->fetch_add(1)) < rrs.size();)
          generate(rrs[schedule[i]]);
        std::ofstream sf(symfile(k));
        for (auto t = taskmgr.first(); t != taskmgr.plast(); ++t)
          for (auto& sym : t->symbols()->symbols())
            sf << t->label() << '\t' << sym << '\n';
        sf.close();
        if (!sf) throw std::runtime_error("failed to write " + symfile(k));
      } catch (std::exception& e) {
        std::cerr << "generate_rr_code: worker " << k << ": " << e.what()
                  << std::endl;
        status = 1;
      } catch (...) {
        std::cerr << "generate_rr_code: worker " << k << ": unknown exception"
                  << std::endl;
        status = 1;
      }
      os.flush();
      std::cout.flush();
      std::cerr.flush();
      _exit(status);  // skip the destructors of the singletons and the iface
    }
    workers.push_back(pid);
  }

  // wait for all workers, even if one failed, since they share the counter
  bool worker_failed = false;
  for (auto pid : workers) {
    int status = 0;
    pid_t result;
    do {
      result = waitpid(pid, &status, 0);
    } while (result < 0 && errno == EINTR);
    if (result < 0 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
      worker_failed = true;
  }
  munmap(shm, sizeof(counter_type));

  // merge the external symbols of the tasks
  std::map<std::string, std::shared_ptr<TaskExternSymbols> > tsymbols;
  for (auto t = taskmgr.first(); t != taskmgr.plast(); ++t)
    tsymbols[t->label()] = t->symbols();
  for (std::size_t k = 0; k != workers.size(); ++k) {
    std::ifstream sf(symfile(k));
    std::string line;
    while (!worker_failed && std::getline(sf, line)) {
      const auto tab = line.find('\t');
      tsymbols.at(line.substr(0, tab))
          ->add(TaskExternSymbols::SymbolList{line.substr(tab + 1)});
    }
    sf.close();
    std::remove(symfile(k).c_str());
  }

  if (fork_failed)
    throw std::runtime_error(
        std::string("generate_rr_code: fork failed (") +
        std::strerror(fork_errno) +
        "); set LIBINT_GENERATOR_JOBS to fewer processes, or to 1");
  if (worker_failed)
    throw std::runtime_error(
        "generate_rr_code: a process generating the set-level RRs failed");
#else
  os << "generate_rr_code: LIBINT_GENERATOR_JOBS=" << njobs
     << " is not supported on this platform, generating the RRs serially"
     << std::endl;
  for (auto& rr : rrs) generate(rr);
#endif
}

}  // namespace

void generate_rr_code(std::ostream& os,
                      const std::shared_ptr<CompilationParameters>& cparams,
                      std::deque<std::string>& decl_filenames,
                      std::deque<std::string>& def_filenames) {
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  ImplicitDimensions::set_default_dims(cparams);
  std::string prefix(cparams->source_directory());

  std::shared_ptr<RRStack> rrstack = RRStack::Instance();

  //
  // generate code for all recurrence relation actually used
  //
  // 1) merge RR lists from all tasks
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.plast();
  std::set<TaskExternSymbols::RRList::value_type> aggregate_rrlist;
  for (tciter t = taskmgr.first(); t != tend; ++t) {
    const std::shared_ptr<TaskExternSymbols> tsymbols = t->symbols();
    auto rrlist = tsymbols->rrlist();
    aggregate_rrlist.insert(rrlist.begin(), rrlist.end());
  }

  // resolve the RR ids in one pass over the stack (find_hashed is a linear
  // search, hence calling it per RR is O(N^2) in the number of RRs), in the
  // order in which the RRs were created
  std::vector<
      std::pair<TaskExternSymbols::RRid, std::shared_ptr<RecurrenceRelation> > >
      id_rrs;
  id_rrs.reserve(aggregate_rrlist.size());
  for (auto it = rrstack->begin(); it != rrstack->end(); ++it)
    if (aggregate_rrlist.count(it->second.first))
      id_rrs.emplace_back(it->second.first, it->second.second);
  if (id_rrs.size() != aggregate_rrlist.size())
    throw std::logic_error(
        "generate_rr_code: " +
        std::to_string(aggregate_rrlist.size() - id_rrs.size()) +
        " of the RRs used by the tasks are not on the RRStack");
  std::sort(id_rrs.begin(), id_rrs.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  std::vector<std::shared_ptr<RecurrenceRelation> > rrs;
  rrs.reserve(id_rrs.size());
  for (auto& id_rr : id_rrs) rrs.push_back(id_rr.second);
  reorder_rrs(rrs, os);

  // generates the code for one RR
  auto generate = [&](const std::shared_ptr<RecurrenceRelation>& rr) {
    os << "generating code for " << context->label_to_name(rr->label())
       << " target=" << rr->rr_target()->label() << endl;

    const std::string filename =
        prefix + context->label_to_filename(rr->label());
    std::basic_ofstream<char> declfile(filename + ".h");
    std::basic_ofstream<char> deffile(filename + ".cc");

    rr->generate_code(context, ImplicitDimensions::default_dims(), rr->label(),
                      declfile, deffile);

    declfile.close();
    deffile.close();
    if (!declfile || !deffile)
      throw std::runtime_error("generate_rr_code: failed to write " + filename +
                               ".{h,cc}");

    // Remove RR to save resources
    rrstack->remove(rr);
    // purge SingletonStacks, to save resources
    PurgeableStacks::Instance()->purge();
  };

  // the file names are assigned here, before any code is generated:
  // label_to_filename() disambiguates names that differ only in case in the
  // order in which it sees them, hence must see them in the same order
  // regardless of the order of generation and of the number of workers
  for (auto& rr : rrs) {
    const std::string filename =
        prefix + context->label_to_filename(rr->label());
    decl_filenames.push_back(filename + ".h");
    def_filenames.push_back(filename + ".cc");
  }

  const auto njobs = std::min<std::size_t>(generator_jobs(os), rrs.size());
  if (njobs <= 1) {
    for (auto& rr : rrs) generate(rr);
  } else {
    generate_in_parallel(os, rrs, njobs, generate);
  }
}

};  // namespace libint2
