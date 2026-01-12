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

#include <deque>
#include <fstream>
#include <iostream>
#include <set>

using namespace std;
using namespace libint2;

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

  for (auto& rrid : aggregate_rrlist) {
    auto rr = rrstack->find_hashed(rrid).second;
    assert(rr);
#endif
    std::string rrlabel = cparams->api_prefix() + rr->label();
    os << "generating code for " << context->label_to_name(rrlabel)
       << " target=" << rr->rr_target()->label() << endl;

    std::string decl_filename(prefix + context->label_to_name(rrlabel));
    decl_filename += ".h";
    std::string def_filename(prefix + context->label_to_name(rrlabel));
    def_filename += ".cc";
    std::basic_ofstream<char> declfile(decl_filename.c_str());
    std::basic_ofstream<char> deffile(def_filename.c_str());

    rr->generate_code(context, ImplicitDimensions::default_dims(), rrlabel,
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
