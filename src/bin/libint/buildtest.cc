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
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

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

void generate_rr_code(std::ostream& os,
                      const std::shared_ptr<CompilationParameters>& cparams,
                      std::deque<std::string>& decl_filenames,
                      std::deque<std::string>& def_filenames) {
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  ImplicitDimensions::set_default_dims(cparams);
  std::string prefix(cparams->source_directory());

  std::shared_ptr<RRStack> rrstack = RRStack::Instance();

#define GENERATE_ALL_RRS 0
#if GENERATE_ALL_RRS
  //
  // generate explicit code for all recurrence relation that were not inlined
  //
  RRStack::citer_type it = rrstack->begin();
  while (it != rrstack->end()) {
    std::shared_ptr<RecurrenceRelation> rr = (*it).second.second;
#else
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
  for (auto& rr : rrs) {
#endif
    os << "generating code for " << context->label_to_name(rr->label())
       << " target=" << rr->rr_target()->label() << endl;

    std::string decl_filename(prefix + context->label_to_filename(rr->label()));
    decl_filename += ".h";
    std::string def_filename(prefix + context->label_to_filename(rr->label()));
    def_filename += ".cc";
    std::basic_ofstream<char> declfile(decl_filename.c_str());
    std::basic_ofstream<char> deffile(def_filename.c_str());

    rr->generate_code(context, ImplicitDimensions::default_dims(), rr->label(),
                      declfile, deffile);

    declfile.close();
    deffile.close();
    decl_filenames.push_back(decl_filename);
    def_filenames.push_back(def_filename);

    // Remove RR to save resources
    rrstack->remove(rr);
    // purge SingletonStacks, to save resources
    PurgeableStacks::Instance()->purge();
  }
}

};  // namespace libint2
