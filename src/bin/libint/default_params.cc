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

#include <default_params.h>
#include <libint2/config.h>
#include <task.h>

#include <cassert>

using namespace libint2;
using namespace std;

const std::string CompilationParameters::Defaults::source_directory("./");
const std::string CompilationParameters::Defaults::api_prefix("");
const std::string CompilationParameters::Defaults::realtype("double");
const std::string CompilationParameters::Defaults::task_name("default");

CompilationParameters::CompilationParameters()
    : default_task_name_(Defaults::task_name),
      max_vector_length_(Defaults::max_vector_length),
      vectorize_by_line_(Defaults::vectorize_by_line),
      align_size_(Defaults::align_size),
      unroll_threshold_(Defaults::unroll_threshold),
      source_directory_(Defaults::source_directory),
      api_prefix_(Defaults::api_prefix),
      single_evaltype_(Defaults::single_evaltype),
      use_C_linking_(Defaults::use_C_linking),
      count_flops_(Defaults::count_flops),
      profile_(Defaults::profile),
      accumulate_targets_(Defaults::accumulate_targets),
      realtype_(Defaults::realtype),
      contracted_targets_(Defaults::contracted_targets) {
  add_task(Defaults::task_name);
}

CompilationParameters::~CompilationParameters() {}

void CompilationParameters::print(std::ostream& os) const {
  using namespace std;
  os << "MAX_AM           = " << max_am(default_task_name()) << endl;
  os << "OPT_AM           = " << max_am_opt(default_task_name()) << endl;

  typedef LibraryTaskManager::TasksCIter citer;
  const LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  for (citer t = taskmgr.first(); t != taskmgr.plast(); ++t) {
    const std::string& tlabel = t->label();
    os << "Task " << tlabel << ":" << endl;
    os << "  MAX_AM         = " << max_am(tlabel) << endl;
    os << "  OPT_AM         = " << max_am_opt(tlabel) << endl;
    os << "  NUM_BF         = " << num_bf(tlabel) << endl;
  }

  os << "MAX_VECTOR_LENGTH    = " << max_vector_length() << endl;
  if (max_vector_length() > 1)
    os << "VECTORIZE_BY_LINE    = " << (vectorize_by_line() ? "true" : "false")
       << endl;
  if (align_size() > 0) os << "ALIGN_SIZE           = " << align_size() << endl;
  os << "UNROLL_THRESH        = " << unroll_threshold() << endl;
  os << "SOURCE_DIRECTORY     = " << source_directory() << endl;
  os << "API_PREFIX           = " << api_prefix() << endl;
  os << "USE_C_LINKING        = " << (use_C_linking() ? "true" : "false")
     << endl;
  os << "COUNT_FLOPS          = " << (count_flops() ? "true" : "false") << endl;
  os << "PROFILE              = " << (profile() ? "true" : "false") << endl;
  os << "ACCUMULATE_TARGETS   = " << (accumulate_targets() ? "true" : "false")
     << endl;
  os << "REALTYPE             = " << (realtype()) << endl;
  os << "CONTRACTED_TARGETS   = " << (contracted_targets() ? "true" : "false")
     << endl;
  os << endl;
}

void CompilationParameters::task_exists(const std::string& t) const {
  if (t != default_task_name()) {
    const LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    // Will throw if task manager doesn't know anything about this task
    taskmgr.find(t);
  }
}

unsigned int CompilationParameters::max_am(std::string t,
                                           unsigned int c) const {
  if (t.empty()) t = default_task_name();
  task_exists(t);

  typedef std::map<std::string, TaskParameters>::const_iterator citer;
  citer ti = task_params_.find(t);
  auto max_am = (ti != task_params_.end())
                    ? ti->second.max_am
                    : task_params_.find(default_task_name())->second.max_am;
  return (c < max_am.size()) ? max_am[c] : max_am[0];
}

unsigned int CompilationParameters::max_am_opt(std::string t) const {
  if (t.empty()) t = default_task_name();
  task_exists(t);

  typedef std::map<std::string, TaskParameters>::const_iterator citer;
  citer ti = task_params_.find(t);
  if (ti != task_params_.end())
    return ti->second.max_am_opt;
  else
    return task_params_.find(default_task_name())->second.max_am_opt;
}

unsigned int CompilationParameters::num_bf(std::string t) const {
  if (t.empty()) t = default_task_name();
  task_exists(t);

  typedef std::map<std::string, TaskParameters>::const_iterator citer;
  citer ti = task_params_.find(t);
  if (ti != task_params_.end())
    return ti->second.num_bf;
  else
    return task_params_.find(default_task_name())->second.num_bf;
}

void CompilationParameters::add_task(const std::string& t) {
  TaskParameters tp;
  // copy defaults from the default task
  if (t != default_task_name()) {
    assert(task_params_.find(default_task_name()) != task_params_.end());
    tp = TaskParameters(task_params_.find(default_task_name())->second);
  }
  task_params_.insert(std::make_pair(t, tp));
}

void CompilationParameters::max_am(const std::string& t, unsigned int ma,
                                   unsigned int c) {
  task_exists(t);

  typedef std::map<std::string, TaskParameters>::iterator iter;
  iter ti = task_params_.find(t);
  if (ti != task_params_.end()) {
    if (ti->second.max_am.size() <= c) ti->second.max_am.resize(c + 1);
    ti->second.max_am[c] = ma;
  } else {
    add_task(t);
    max_am(t, ma, c);
  }
  std::cout << "CompilationParameters::max_am: task=" << t << " max_am=" << ma
            << " center=" << c << std::endl;
}

void CompilationParameters::max_am_opt(const std::string& t, unsigned int v) {
  task_exists(t);

  typedef std::map<std::string, TaskParameters>::iterator iter;
  iter ti = task_params_.find(t);
  if (ti != task_params_.end())
    ti->second.max_am_opt = v;
  else {
    add_task(t);
    max_am_opt(t, v);
  }
  std::cout << "CompilationParameters::max_am_opt: task=" << t
            << " max_am_opt=" << v << std::endl;
}

void CompilationParameters::num_bf(const std::string& t, unsigned int nbf) {
  task_exists(t);

  typedef std::map<std::string, TaskParameters>::iterator iter;
  iter ti = task_params_.find(t);
  if (ti != task_params_.end())
    ti->second.num_bf = nbf;
  else {
    add_task(t);
    num_bf(t, nbf);
  }
}

//////////

TaskParameters::TaskParameters()
    : max_ntarget_(1),
      max_stack_size_(1, 1),
      max_vector_stack_size_(1, 0),
      max_hrr_hsrank_(1, 0),
      max_hrr_lsrank_(1, 0) {}

//////////

std::string libint2::label_to_funcname(const std::string& label) {
  // Do not prepend compute as it messes up the API prefix functionality.
#if 0
  std::string result("compute");
  result += label;
#else
  const std::string& result = label;
#endif
  return result;
}

bool libint2::condense_expr(unsigned int unroll_threshold, bool vectorize) {
  bool condense_expr = unroll_threshold > 0 && vectorize;
  return condense_expr;
}
