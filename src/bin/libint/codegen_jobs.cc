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

#include <codegen_jobs.h>
#include <comp_deriv_gauss.h>
#include <context.h>
#include <default_params.h>
#include <dims.h>
#include <iface.h>
#include <rr.h>
#include <task.h>

#include <array>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <thread>

#if defined(__unix__) || defined(__APPLE__)
#define LIBINT_GENERATOR_CAN_FORK 1
#include <poll.h>
#include <signal.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#if defined(__linux__)
#include <sched.h>
#endif
#else
#define LIBINT_GENERATOR_CAN_FORK 0
#endif

namespace libint2 {

namespace {

/// @return the number of hardware threads that this process may use (on Linux
/// this respects the CPU affinity, e.g. taskset or the cpuset of a container)
std::size_t available_hardware_threads() {
#if defined(__linux__)
  cpu_set_t set;
  CPU_ZERO(&set);
  if (sched_getaffinity(0, sizeof(set), &set) == 0) {
    const int n = CPU_COUNT(&set);
    if (n > 0) return static_cast<std::size_t>(n);
  }
#endif
  return std::max(1u, std::thread::hardware_concurrency());
}

}  // namespace

std::size_t generator_jobs() {
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
  if (njobs == 0)
    return std::min(available_hardware_threads(), max_automatic_generator_jobs);
  return static_cast<std::size_t>(njobs);
}

CodegenJobs& CodegenJobs::instance() {
  static CodegenJobs jobs;
  return jobs;
}

void CodegenJobs::add(const std::string& task, double cost, Job job) {
  jobs_.push_back(Entry{task, cost, std::move(job)});
}

namespace detail {

NameMap final_filenames(const std::set<std::string>& names) {
  std::map<std::string, std::vector<std::string> > folded;
  for (auto& name : names) folded[fold_case(name)].push_back(name);  // sorted
  NameMap result;
  for (auto& f : folded)
    for (std::size_t i = 0; i != f.second.size(); ++i)
      result[f.second[i]] = disambiguated_filename(f.second[i], i);
  return result;
}

FilenameMapper::FilenameMapper(const std::string& d, const NameMap& local,
                               const NameMap& final)
    : dir(d) {
  for (auto& nf : local) {
    const auto& f = final.at(nf.first);
    local_to_final[nf.second] = f;
    if (f != nf.second) renames = true;
  }
}

std::string FilenameMapper::map_include(const std::string& line,
                                        const std::string& prefix) const {
  static const std::string open("#include <");
  if (line.compare(0, open.size(), open) != 0) return line;
  const auto close = line.find('>', open.size());
  if (close == std::string::npos) return line;
  std::string path = line.substr(open.size(), close - open.size());
  std::string head;  // directory part of the path
  const bool in_dir = !dir.empty() && path.compare(0, dir.size(), dir) == 0;
  if (in_dir) {
    path = path.substr(dir.size());
    head = prefix;
  } else {
    const auto slash = path.rfind('/');
    if (slash != std::string::npos) {
      head = path.substr(0, slash + 1);
      path = path.substr(slash + 1);
    }
  }
  const auto dot = path.rfind('.');
  const auto it = (dot == std::string::npos)
                      ? local_to_final.end()
                      : local_to_final.find(path.substr(0, dot));
  if (it == local_to_final.end()) {
    // not a file of this process: only the directory is mapped
    if (!in_dir) return line;
    return open + head + path + line.substr(close);
  }
  return open + head + it->second + path.substr(dot) + line.substr(close);
}

void finalize_files(const FilenameMapper& mapper, const NameMap& names,
                    const std::string& prefix) {
  const bool in_place = (mapper.dir == prefix);
  std::vector<std::pair<std::string, std::string> > moves;  // from, to
  for (auto& nf : names) {
    const auto& local = nf.second;
    const auto& final = mapper.local_to_final.at(local);
    for (const char* ext : {".h", ".cc"}) {
      const std::string from = mapper.dir + local + ext;
      std::ifstream test(from);
      if (!test.good()) continue;
      test.close();
      moves.emplace_back(from, prefix + final + ext);
    }
  }
  // rewrite the include directives if any file name changes
  if (mapper.renames) {
    for (auto& m : moves) {
      std::ifstream in(m.first);
      std::ostringstream out;
      std::string line;
      bool changed = false;
      while (std::getline(in, line)) {
        const auto mapped = mapper.map_include(line, "");
        changed = changed || (mapped != line);
        out << mapped << '\n';
      }
      in.close();
      if (changed) {
        std::ofstream o(m.first);
        o << out.str();
        if (!o) throw std::runtime_error("failed to rewrite " + m.first);
      }
    }
  }
  // in place, file names may be swapped, hence rename in two steps
  if (in_place) {
    std::vector<std::pair<std::string, std::string> > tmp;
    for (auto& m : moves)
      if (m.first != m.second) {
        const auto t = m.second + ".libint-rename";
        if (std::rename(m.first.c_str(), t.c_str()) != 0)
          throw std::runtime_error("failed to rename " + m.first);
        tmp.emplace_back(t, m.second);
      }
    moves = tmp;
  }
  for (auto& m : moves)
    if (m.first != m.second && std::rename(m.first.c_str(), m.second.c_str()))
      throw std::runtime_error("failed to rename " + m.first + " to " +
                               m.second);
}

}  // namespace detail

namespace {

using detail::FilenameMapper;
using detail::NameMap;

/// @return the file names of the names in the registry
NameMap registry_filenames() {
  NameMap result;
  for (auto& folded_names : filename_registry()) {
    const auto& names = folded_names.second;
    for (std::size_t i = 0; i != names.size(); ++i)
      result[names[i]] = disambiguated_filename(names[i], i);
  }
  return result;
}

void apply_journal(const ShellSetJournal& journal, const std::string& task,
                   Libint2Iface& iface, const FilenameMapper& mapper,
                   const std::string& prefix) {
  for (auto& s : journal.static_init) iface.to_static_init(s);
  for (auto& s : journal.int_iface) {
    // each line ends with a newline
    std::string line = s;
    std::string eol;
    while (!line.empty() && line.back() == '\n') {
      line.pop_back();
      eol += '\n';
    }
    iface.to_int_iface(mapper.map_include(line, prefix) + eol);
  }
  auto params = LibraryTaskManager::Instance().find(task)->params();
  for (auto& am_size : journal.max_stack_size)
    params->max_stack_size(am_size.first, am_size.second);
  for (auto n : journal.max_ntarget) params->max_ntarget(n);
}

//
// serialization of the reports of the workers: strings are written as
// <length> <characters>
//
void write_string(std::ostream& os, const std::string& s) {
  os << s.size() << ' ' << s;
}
std::string read_string(std::istream& is) {
  std::size_t n = 0;
  if (!(is >> n)) throw std::runtime_error("failed to read a worker report");
  is.get();  // the separator
  // read in chunks, so that a corrupt length cannot allocate more memory than
  // the report holds
  constexpr std::size_t chunk = 1 << 16;
  std::string s;
  while (s.size() < n) {
    const auto offset = s.size();
    const auto m = std::min(chunk, n - offset);
    s.resize(offset + m);
    if (!is.read(&s[offset], m))
      throw std::runtime_error("failed to read a worker report");
  }
  return s;
}
template <typename T>
T read_value(std::istream& is) {
  T v{};
  is >> v;
  if (!is) throw std::runtime_error("failed to read a worker report");
  return v;
}

void write_strings(std::ostream& os, const std::vector<std::string>& v) {
  os << v.size() << ' ';
  for (auto& s : v) write_string(os, s);
}
std::vector<std::string> read_strings(std::istream& is) {
  const auto n = read_value<std::size_t>(is);
  std::vector<std::string> v;  // not reserved: n may be corrupt
  for (std::size_t i = 0; i != n; ++i) v.push_back(read_string(is));
  return v;
}

void write_journal(std::ostream& os, const ShellSetJournal& j) {
  write_strings(os, j.static_init);
  write_strings(os, j.int_iface);
  os << j.max_stack_size.size() << ' ';
  for (auto& p : j.max_stack_size) os << p.first << ' ' << p.second << ' ';
  os << j.max_ntarget.size() << ' ';
  for (auto n : j.max_ntarget) os << n << ' ';
}
ShellSetJournal read_journal(std::istream& is) {
  ShellSetJournal j;
  j.static_init = read_strings(is);
  j.int_iface = read_strings(is);
  const auto ns = read_value<std::size_t>(is);
  for (std::size_t i = 0; i != ns; ++i) {
    const auto am = read_value<unsigned int>(is);
    const auto size = read_value<unsigned int>(is);
    j.max_stack_size.emplace_back(am, size);
  }
  const auto nt = read_value<std::size_t>(is);
  for (std::size_t i = 0; i != nt; ++i)
    j.max_ntarget.push_back(read_value<unsigned int>(is));
  return j;
}

void write_names(std::ostream& os, const NameMap& names) {
  os << names.size() << ' ';
  for (auto& nf : names) {
    write_string(os, nf.first);
    write_string(os, nf.second);
  }
}
NameMap read_names(std::istream& is) {
  NameMap names;
  const auto n = read_value<std::size_t>(is);
  for (std::size_t i = 0; i != n; ++i) {
    auto name = read_string(is);
    names[name] = read_string(is);
  }
  return names;
}

/// writes the external symbols of each task
void write_task_symbols(std::ostream& os) {
  auto& taskmgr = LibraryTaskManager::Instance();
  os << taskmgr.ntask() << ' ';
  for (auto t = taskmgr.first(); t != taskmgr.plast(); ++t) {
    write_string(os, t->label());
    const auto& symbols = t->symbols()->symbols();
    write_strings(os, std::vector<std::string>(symbols.begin(), symbols.end()));
  }
}
/// reads the external symbols of the tasks and adds them to the tasks
void merge_task_symbols(std::istream& is) {
  auto& taskmgr = LibraryTaskManager::Instance();
  const auto ntasks = read_value<std::size_t>(is);
  for (std::size_t i = 0; i != ntasks; ++i) {
    const auto label = read_string(is);
    const auto symbols = read_strings(is);
    taskmgr.find(label)->symbols()->add(
        TaskExternSymbols::SymbolList(symbols.begin(), symbols.end()));
  }
}

}  // namespace

void CodegenJobs::run(std::ostream& os,
                      const std::shared_ptr<CompilationParameters>& cparams,
                      Libint2Iface& iface, std::size_t njobs,
                      const std::function<void(std::size_t)>& on_worker_start) {
#if LIBINT_GENERATOR_CAN_FORK
  if (njobs > 1 && jobs_.size() > 1)
    run_in_parallel(os, cparams, iface, std::min(njobs, jobs_.size()),
                    on_worker_start);
  else
    run_serially(os, cparams, iface);
#else
  if (njobs > 1)
    os << "LIBINT_GENERATOR_JOBS=" << njobs
       << " is not supported on this platform, generating the code serially"
       << std::endl;
  run_serially(os, cparams, iface);
#endif
  jobs_.clear();
}

void CodegenJobs::run_serially(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    Libint2Iface& iface) {
  auto& taskmgr = LibraryTaskManager::Instance();
  const std::string prefix = cparams->source_directory();

  // the order in which the shell sets are generated must not matter;
  // LIBINT_SHELLSET_ORDER changes it to test this
  std::vector<std::size_t> order(jobs_.size());
  std::iota(order.begin(), order.end(), 0);
  permute_for_testing(order, "LIBINT_SHELLSET_ORDER",
                      [this](std::size_t i) { return i; });

  std::vector<ShellSetJournal> journals(jobs_.size());
  for (auto i : order) {
    taskmgr.current(jobs_[i].task);
    ImplicitDimensions::set_default_dims(cparams);
    jobs_[i].job(journals[i]);
  }

  generate_rr_code(os, cparams);

  const auto local = registry_filenames();
  std::set<std::string> names;
  for (auto& nf : local) names.insert(nf.first);
  const FilenameMapper mapper(prefix, local, detail::final_filenames(names));
  for (std::size_t i = 0; i != jobs_.size(); ++i)
    apply_journal(journals[i], jobs_[i].task, iface, mapper, prefix);
  detail::finalize_files(mapper, local, prefix);
}

#if LIBINT_GENERATOR_CAN_FORK
namespace {

/// writes one byte to fd; @return false on failure
bool send_byte(int fd, char c) {
  ssize_t n;
  do {
    n = write(fd, &c, 1);
  } while (n < 0 && errno == EINTR);
  return n == 1;
}
/// reads one byte from fd; @return 0 on EOF or failure
char receive_byte(int fd) {
  char c = 0;
  ssize_t n;
  do {
    n = read(fd, &c, 1);
  } while (n < 0 && errno == EINTR);
  return n == 1 ? c : 0;
}

/// the wait status of a process that could not be waited for
constexpr int unknown_wait_status = -1;

/// @return how a worker process ended, given its wait status
std::string describe_exit(int status) {
  std::ostringstream oss;
  if (status == unknown_wait_status)
    oss << "ended, but its exit status is unknown";
  else if (WIFEXITED(status)) {
    oss << "exited with status " << WEXITSTATUS(status);
    if (WEXITSTATUS(status) != 0) oss << " (see its error message above)";
  } else if (WIFSIGNALED(status)) {
    const int sig = WTERMSIG(status);
    oss << "was killed by signal " << sig;
    if (const char* name = strsignal(sig)) oss << " (" << name << ")";
    if (sig == SIGKILL)
      oss << ", possibly because the system ran out of memory; set "
             "LIBINT2_GENERATOR_JOBS (CMake) or LIBINT_GENERATOR_JOBS "
             "(environment) to fewer processes";
  } else
    oss << "ended with wait status " << status;
  return oss.str();
}

/// ignores SIGPIPE while it exists (and in the processes forked meanwhile), so
/// that writing to the pipe of a worker that died fails with EPIPE, which the
/// generator reports, rather than killing the generator
class IgnoreSIGPIPE {
 public:
  IgnoreSIGPIPE() {
    struct sigaction ignore;
    std::memset(&ignore, 0, sizeof(ignore));
    ignore.sa_handler = SIG_IGN;
    sigemptyset(&ignore.sa_mask);
    restore_ = (sigaction(SIGPIPE, &ignore, &old_) == 0);
  }
  ~IgnoreSIGPIPE() {
    if (restore_) sigaction(SIGPIPE, &old_, nullptr);
  }
  IgnoreSIGPIPE(const IgnoreSIGPIPE&) = delete;
  IgnoreSIGPIPE& operator=(const IgnoreSIGPIPE&) = delete;

 private:
  struct sigaction old_;
  bool restore_ = false;
};

/// For testing the handling of failed workers (see
/// src/bin/check_codegen_failure.cmake): LIBINT_GENERATOR_TEST_KILL_WORKER=
/// <k>:<point> makes worker k kill itself with SIGKILL (like the out-of-memory
/// killer) at <point>: shellsets (before generating its first shell set) or
/// reported (after reporting its shell sets, before generating the RRs)
bool test_kill_worker(std::size_t k, const char* point) {
  const char* env = std::getenv("LIBINT_GENERATOR_TEST_KILL_WORKER");
  if (env == nullptr) return false;
  return std::string(env) == std::to_string(k) + ":" + point;
}

}  // namespace
#endif

void CodegenJobs::run_in_parallel(
    std::ostream& os, const std::shared_ptr<CompilationParameters>& cparams,
    Libint2Iface& iface, std::size_t njobs,
    const std::function<void(std::size_t)>& on_worker_start) {
#if LIBINT_GENERATOR_CAN_FORK
  namespace fs = std::filesystem;
  auto& taskmgr = LibraryTaskManager::Instance();
  const std::string prefix = cparams->source_directory();
  const std::string tag = ".libint_worker." + std::to_string(getpid()) + ".";
  auto worker_dir = [&](std::size_t k) {
    return prefix + tag + std::to_string(k) + "/";
  };
  auto report_file = [&](std::size_t k, const char* what) {
    return prefix + tag + std::to_string(k) + "." + what;
  };

  // largest first, so that the workers finish at about the same time
  std::vector<std::size_t> schedule(jobs_.size());
  std::iota(schedule.begin(), schedule.end(), 0);
  std::stable_sort(schedule.begin(), schedule.end(),
                   [this](std::size_t a, std::size_t b) {
                     return jobs_[a].cost > jobs_[b].cost;
                   });

  // index of the next job in schedule, shared by the workers
  using counter_type = std::atomic<std::size_t>;
  static_assert(counter_type::is_always_lock_free,
                "the job counter shared by the generator processes must be "
                "lock-free");
  void* shm = mmap(nullptr, sizeof(counter_type), PROT_READ | PROT_WRITE,
                   MAP_SHARED | MAP_ANON, -1, 0);
  if (shm == MAP_FAILED)
    throw std::runtime_error(std::string("CodegenJobs: mmap failed: ") +
                             std::strerror(errno));
  counter_type* next = new (shm) counter_type(0);

  // pipes: worker -> generator ("phase 1 done"), generator -> worker ("go");
  // -1 = not open
  std::vector<std::array<int, 2> > to_parent(njobs, {-1, -1}),
      to_worker(njobs, {-1, -1});
  auto close_fd = [](int& fd) {
    if (fd >= 0) close(fd);
    fd = -1;
  };

  const IgnoreSIGPIPE ignore_sigpipe;

  std::vector<pid_t> workers;  // the process ids of the workers
  std::vector<int> wait_status;
  std::vector<bool> reaped;
  /// waits for worker k to end; @return its wait status
  auto reap = [&](std::size_t k) {
    if (!reaped[k]) {
      int status = 0;
      pid_t result;
      do {
        result = waitpid(workers[k], &status, 0);
      } while (result < 0 && errno == EINTR);
      wait_status[k] = (result == workers[k]) ? status : unknown_wait_status;
      reaped[k] = true;
    }
    return wait_status[k];
  };
  /// @return true if worker k succeeded
  auto succeeded = [&](std::size_t k) {
    const auto status = reap(k);
    return status != unknown_wait_status && WIFEXITED(status) &&
           WEXITSTATUS(status) == 0;
  };
  /// removes the files of the workers (on failure also their directories,
  /// which then hold partial output)
  auto cleanup = [&](bool remove_worker_dirs) {
    munmap(shm, sizeof(counter_type));
    std::error_code ec;
    for (std::size_t k = 0; k != njobs; ++k) {
      close_fd(to_parent[k][0]);
      close_fd(to_parent[k][1]);
      close_fd(to_worker[k][0]);
      close_fd(to_worker[k][1]);
      for (const char* what : {"shellsets", "assignment", "rrs"})
        fs::remove(report_file(k, what), ec);
      if (remove_worker_dirs) fs::remove_all(worker_dir(k), ec);
    }
  };

  std::string failure;  // why the generation failed
  try {
    for (std::size_t k = 0; k != njobs; ++k) {
      if (pipe(to_parent[k].data()) != 0 || pipe(to_worker[k].data()) != 0)
        throw std::runtime_error(std::string("pipe failed: ") +
                                 std::strerror(errno));
      // an earlier generator with the same process id (e.g. in a container)
      // may have left a directory behind if it was killed
      std::error_code ec;
      fs::remove_all(worker_dir(k), ec);
      if (!fs::create_directory(worker_dir(k), ec) || ec)
        throw std::runtime_error("failed to create " + worker_dir(k) + ": " +
                                 ec.message());
    }

    os << "generating the code for " << jobs_.size() << " shell sets in "
       << njobs << " processes" << std::endl;
    os.flush();
    std::cout.flush();
    std::cerr.flush();

    int fork_errno = 0;
    for (std::size_t k = 0; k != njobs; ++k) {
      const pid_t pid = fork();
      if (pid < 0) {
        fork_errno = errno;
        break;
      }
      if (pid > 0) {
        workers.push_back(pid);
        continue;
      }

      //
      // worker k
      //
      int status = 0;
      try {
        for (std::size_t l = 0; l != njobs; ++l) {
          close_fd(to_parent[l][0]);
          close_fd(to_worker[l][1]);
          if (l != k) {
            close_fd(to_parent[l][1]);
            close_fd(to_worker[l][0]);
          }
        }
        cparams->source_directory(worker_dir(k));
        if (on_worker_start) on_worker_start(k);

        // phase 1: shell sets
        if (test_kill_worker(k, "shellsets")) raise(SIGKILL);
        std::vector<std::pair<std::size_t, ShellSetJournal> > done;
        for (std::size_t i; (i = next->fetch_add(1)) < schedule.size();) {
          const auto j = schedule[i];
          taskmgr.current(jobs_[j].task);
          ImplicitDimensions::set_default_dims(cparams);
          done.emplace_back(j, ShellSetJournal{});
          jobs_[j].job(done.back().second);
        }
        {
          std::ofstream report(report_file(k, "shellsets"));
          report << done.size() << ' ';
          for (auto& d : done) {
            report << d.first << ' ';
            write_journal(report, d.second);
          }
          write_task_symbols(report);
          // the RRs that each task needs, with their cost
          auto rrstack = RRStack::Instance();
          std::map<RRStack::InstanceID, std::shared_ptr<RecurrenceRelation> >
              id_to_rr;
          for (auto it = rrstack->begin(); it != rrstack->end(); ++it)
            id_to_rr[it->second.first] = it->second.second;
          report << taskmgr.ntask() << ' ';
          for (auto t = taskmgr.first(); t != taskmgr.plast(); ++t) {
            write_string(report, t->label());
            const auto rrlist = t->symbols()->rrlist();
            report << rrlist.size() << ' ';
            for (auto& id : rrlist) {
              const auto& rr = id_to_rr.at(id);
              write_string(report, rr->label());
              report << rr->rr_target()->size() << ' ';
            }
          }
          report.close();
          if (!report) throw std::runtime_error("failed to write the report");
        }
        if (!send_byte(to_parent[k][1], 'd'))
          throw std::runtime_error("lost the generator");
        if (test_kill_worker(k, "reported")) raise(SIGKILL);
        if (receive_byte(to_worker[k][0]) != 'g')
          throw std::runtime_error("the generator stopped the worker");

        // phase 2: the RRs assigned to this worker
        {
          std::ifstream assignment(report_file(k, "assignment"));
          auto rrstack = RRStack::Instance();
          std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
          ImplicitDimensions::set_default_dims(cparams);
          const auto nrrs = read_value<std::size_t>(assignment);
          for (std::size_t r = 0; r != nrrs; ++r) {
            const auto label = read_string(assignment);
            const auto tasks = read_strings(assignment);
            const auto& id_rr = rrstack->find(label);
            if (!id_rr.second)
              throw std::runtime_error("RR " + label + " not found");
            // the tasks that use the RR get its external symbols
            for (auto& task : tasks)
              taskmgr.find(task)->symbols()->add(
                  TaskExternSymbols::RRList{id_rr.first});
            generate_rr(os, context, worker_dir(k), id_rr.second);
          }
        }
        {
          std::ofstream report(report_file(k, "rrs"));
          write_task_symbols(report);
          write_names(report, registry_filenames());
          // the instances of GenericGaussDeriv used by the code of this
          // worker; GenericGaussDeriv.cc is written by the parent at exit
          // (this process exits without running the destructors)
          const auto& generic_instances =
              CR_DerivGauss_GenericInstantiator::instance().instances();
          report << generic_instances.size() << ' ';
          for (auto& i : generic_instances)
            report << i.first << ' ' << (i.second ? 1 : 0) << ' ';
          report.close();
          if (!report) throw std::runtime_error("failed to write the report");
        }
      } catch (std::exception& e) {
        std::cerr << "code generator worker " << k << ": " << e.what()
                  << std::endl;
        status = 1;
      } catch (...) {
        std::cerr << "code generator worker " << k << ": unknown exception"
                  << std::endl;
        status = 1;
      }
      os.flush();
      std::cout.flush();
      std::cerr.flush();
      _exit(status);  // skip the destructors of the singletons and the iface
    }

    //
    // the generator
    //
    wait_status.assign(workers.size(), unknown_wait_status);
    reaped.assign(workers.size(), false);
    for (std::size_t k = 0; k != njobs; ++k) {
      close_fd(to_parent[k][1]);
      close_fd(to_worker[k][0]);
    }
    if (workers.size() != njobs)
      throw std::runtime_error(
          std::string("fork failed (") + std::strerror(fork_errno) +
          "); set LIBINT_GENERATOR_JOBS to fewer processes, or to 1");

    // phase 1: wait until all workers have reported their shell sets; a worker
    // that ends without reporting fails the generation at once
    {
      std::vector<bool> reported(njobs, false);
      std::size_t nreported = 0;
      while (nreported != njobs) {
        std::vector<pollfd> fds;
        std::vector<std::size_t> fd_worker;
        for (std::size_t k = 0; k != njobs; ++k)
          if (!reported[k]) {
            fds.push_back(pollfd{to_parent[k][0], POLLIN, 0});
            fd_worker.push_back(k);
          }
        if (poll(fds.data(), fds.size(), -1) < 0) {
          if (errno == EINTR) continue;
          throw std::runtime_error(std::string("poll failed: ") +
                                   std::strerror(errno));
        }
        for (std::size_t i = 0; i != fds.size(); ++i) {
          if (fds[i].revents == 0) continue;
          const auto k = fd_worker[i];
          if (receive_byte(to_parent[k][0]) != 'd')
            throw std::runtime_error("process " + std::to_string(k) +
                                     " generating the shell sets " +
                                     describe_exit(reap(k)));
          reported[k] = true;
          ++nreported;
        }
      }
    }

    // collect the journals, the symbols, and the RRs
    std::vector<ShellSetJournal> journals(jobs_.size());
    std::vector<std::size_t> journal_worker(jobs_.size(), njobs);
    struct RRInfo {
      std::size_t cost = 0;
      std::set<std::size_t> workers;
      std::vector<std::string> tasks;
    };
    std::map<std::string, RRInfo> rrs;
    for (std::size_t k = 0; k != njobs; ++k) {
      std::ifstream report(report_file(k, "shellsets"));
      const auto ndone = read_value<std::size_t>(report);
      for (std::size_t d = 0; d != ndone; ++d) {
        const auto j = read_value<std::size_t>(report);
        journals.at(j) = read_journal(report);
        journal_worker.at(j) = k;
      }
      merge_task_symbols(report);
      const auto ntasks = read_value<std::size_t>(report);
      for (std::size_t t = 0; t != ntasks; ++t) {
        const auto task = read_string(report);
        const auto nrrs = read_value<std::size_t>(report);
        for (std::size_t r = 0; r != nrrs; ++r) {
          const auto label = read_string(report);
          auto& info = rrs[label];
          info.cost = read_value<std::size_t>(report);
          info.workers.insert(k);
          if (std::find(info.tasks.begin(), info.tasks.end(), task) ==
              info.tasks.end())
            info.tasks.push_back(task);
        }
      }
    }
    for (auto k : journal_worker)
      if (k == njobs) throw std::runtime_error("a shell set was not generated");

    // phase 2: assign each RR to the least loaded worker that has it, largest
    // first
    std::vector<std::pair<std::string, const RRInfo*> > by_cost;
    for (auto& r : rrs) by_cost.emplace_back(r.first, &r.second);
    std::stable_sort(by_cost.begin(), by_cost.end(),
                     [](const auto& a, const auto& b) {
                       return a.second->cost > b.second->cost;
                     });
    std::vector<std::size_t> load(njobs, 0);
    std::vector<std::ostringstream> assignments(njobs);
    std::vector<std::size_t> nassigned(njobs, 0);
    for (auto& r : by_cost) {
      std::size_t best = njobs;
      for (auto k : r.second->workers)
        if (best == njobs || load[k] < load[best]) best = k;
      load[best] += r.second->cost + 1;
      ++nassigned[best];
      write_string(assignments[best], r.first);
      write_strings(assignments[best], r.second->tasks);
    }
    os << "generating the code for " << rrs.size() << " RRs in " << njobs
       << " processes" << std::endl;
    for (std::size_t k = 0; k != njobs; ++k) {
      std::ofstream assignment(report_file(k, "assignment"));
      assignment << nassigned[k] << ' ' << assignments[k].str();
      assignment.close();
      if (!assignment)
        throw std::runtime_error("failed to write " +
                                 report_file(k, "assignment"));
      // fails with EPIPE (SIGPIPE is ignored) if the worker died
      if (!send_byte(to_worker[k][1], 'g'))
        throw std::runtime_error("process " + std::to_string(k) +
                                 " generating the RRs " +
                                 describe_exit(reap(k)));
    }

    for (std::size_t k = 0; k != njobs; ++k)
      if (!succeeded(k))
        throw std::runtime_error("process " + std::to_string(k) +
                                 " generating the RRs " +
                                 describe_exit(reap(k)));

    // the external symbols and the file names of each worker
    std::vector<NameMap> local(njobs);
    for (std::size_t k = 0; k != njobs; ++k) {
      std::ifstream report(report_file(k, "rrs"));
      merge_task_symbols(report);
      local[k] = read_names(report);
      const auto ngeneric = read_value<std::size_t>(report);
      for (std::size_t i = 0; i != ngeneric; ++i) {
        const auto L = read_value<unsigned int>(report);
        const auto vectorize = read_value<int>(report);
        CR_DerivGauss_GenericInstantiator::instance().add(L, vectorize != 0);
      }
    }

    // file names that differ only in case get suffixes in the lexicographic
    // order of the names
    std::set<std::string> names;
    for (auto& l : local)
      for (auto& nf : l) names.insert(nf.first);
    const auto final = detail::final_filenames(names);
    std::vector<FilenameMapper> mappers;
    for (std::size_t k = 0; k != njobs; ++k)
      mappers.emplace_back(worker_dir(k), local[k], final);

    for (std::size_t j = 0; j != jobs_.size(); ++j)
      apply_journal(journals[j], jobs_[j].task, iface,
                    mappers[journal_worker[j]], prefix);
    for (std::size_t k = 0; k != njobs; ++k) {
      detail::finalize_files(mappers[k], local[k], prefix);
      // every file that a worker wrote must have been moved into the source
      // directory, else the library would silently lack it
      std::error_code ec;
      std::vector<std::string> left;
      for (auto& entry : fs::directory_iterator(worker_dir(k), ec))
        left.push_back(entry.path().filename().string());
      if (!left.empty()) {
        std::sort(left.begin(), left.end());
        std::string list;
        for (auto& f : left) list += " " + f;
        throw std::runtime_error("files written by process " +
                                 std::to_string(k) +
                                 " are not known to the generator:" + list);
      }
    }
  } catch (std::exception& e) {
    failure = e.what();
  }

  if (!failure.empty()) {
    // stop the workers that still run, and remove their partial output
    for (std::size_t k = 0; k != workers.size(); ++k)
      if (!reaped[k]) kill(workers[k], SIGTERM);
    for (std::size_t k = 0; k != workers.size(); ++k) reap(k);
    cleanup(true);
    throw std::runtime_error("CodegenJobs: " + failure);
  }
  cleanup(true);  // the worker directories are empty now
#else
  run_serially(os, cparams, iface);
#endif
}

}  // namespace libint2
