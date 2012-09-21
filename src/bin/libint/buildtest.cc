
#include <iostream>
#include <fstream>
#include <deque>
#include <set>
#include <rr.h>
#include <context.h>
#include <dims.h>
#include <task.h>

using namespace libint2;

namespace libint2 {

  void
  generate_rr_code(std::ostream& os,
                   const SafePtr<CompilationParameters>& cparams,
                   std::deque<std::string>& decl_filenames,
                   std::deque<std::string>& def_filenames)
  {
    SafePtr<CodeContext> context(new CppCodeContext(cparams));
    ImplicitDimensions::set_default_dims(cparams);
    std::string prefix(cparams->source_directory());

    SafePtr<RRStack> rrstack = RRStack::Instance();

#define GENERATE_ALL_RRS 0
#if GENERATE_ALL_RRS
    //
    // generate explicit code for all recurrence relation that were not inlined
    //
    RRStack::citer_type it = rrstack->begin();
    while (it != rrstack->end()) {
      SafePtr<RecurrenceRelation> rr = (*it).second.second;
#else
    //
    // generate code for all recurrence relation actually used
    //
    // 1) merge RR lists from all tasks
    LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
    typedef LibraryTaskManager::TasksCIter tciter;
    const tciter tend = taskmgr.plast();
    std::set<TaskExternSymbols::RRList::value_type> aggregate_rrlist;
    for(tciter t=taskmgr.first(); t!=tend; ++t) {
      const SafePtr<TaskExternSymbols> tsymbols = t->symbols();
      typedef TaskExternSymbols::SymbolList SymbolList;
      auto rrlist = tsymbols->rrlist();
      aggregate_rrlist.insert(rrlist.begin(), rrlist.end());
    }

    for(auto& rrid: aggregate_rrlist) {
      auto rr = rrstack->find_hashed(rrid).second;
      assert(rr);
#endif
      std::string rrlabel = cparams->api_prefix() + rr->label();
      os << "generating code for " << context->label_to_name(rrlabel) << " target=" << rr->rr_target()->label() << endl;

      std::string decl_filename(prefix + context->label_to_name(rrlabel));  decl_filename += ".h";
      std::string def_filename(prefix + context->label_to_name(rrlabel));  def_filename += ".cc";
      std::basic_ofstream<char> declfile(decl_filename.c_str());
      std::basic_ofstream<char> deffile(def_filename.c_str());

      rr->generate_code(context,ImplicitDimensions::default_dims(),rrlabel,declfile,deffile);

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

};
