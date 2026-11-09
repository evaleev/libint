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
#include <context.h>
#include <dims.h>
#include <purgeable.h>
#include <rr.h>
#include <task.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace libint2 {

std::vector<std::shared_ptr<RecurrenceRelation> > rrs_to_generate() {
  std::shared_ptr<RRStack> rrstack = RRStack::Instance();

  // merge RR lists from all tasks
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  std::set<TaskExternSymbols::RRList::value_type> aggregate_rrlist;
  for (auto t = taskmgr.first(); t != taskmgr.plast(); ++t) {
    auto rrlist = t->symbols()->rrlist();
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

  // the generated code must not depend on the order in which the RRs are
  // generated; LIBINT_RR_ORDER changes the order to test this (see
  // src/bin/check_codegen_order.cmake)
  permute_for_testing(rrs, "LIBINT_RR_ORDER",
                      [](const auto& rr) { return rr->label(); });
  return rrs;
}

void generate_rr(std::ostream& os, const std::shared_ptr<CodeContext>& context,
                 const std::string& prefix,
                 const std::shared_ptr<RecurrenceRelation>& rr) {
  os << "generating code for " << context->label_to_name(rr->label())
     << " target=" << rr->rr_target()->label() << std::endl;

  const std::string filename = prefix + context->label_to_filename(rr->label());
  std::basic_ofstream<char> declfile(filename + ".h");
  std::basic_ofstream<char> deffile(filename + ".cc");

  rr->generate_code(context, ImplicitDimensions::default_dims(), rr->label(),
                    declfile, deffile);

  declfile.close();
  deffile.close();
  if (!declfile || !deffile)
    throw std::runtime_error("generate_rr: failed to write " + filename +
                             ".{h,cc}");

  // Remove RR to save resources
  RRStack::Instance()->remove(rr);
  // purge SingletonStacks, to save resources
  PurgeableStacks::Instance()->purge();
}

void generate_rr_code(std::ostream& os,
                      const std::shared_ptr<CompilationParameters>& cparams) {
  std::shared_ptr<CodeContext> context(new CppCodeContext(cparams));
  ImplicitDimensions::set_default_dims(cparams);
  const std::string prefix(cparams->source_directory());
  for (auto& rr : rrs_to_generate()) generate_rr(os, context, prefix, rr);
}

};  // namespace libint2
