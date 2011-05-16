
#include <iostream>
#include <fstream>
#include <rr.h>
#include <context.h>
#include <dims.h>

using namespace libint2;

namespace libint2 {

  void
  generate_rr_code(std::ostream& os,
                   const SafePtr<CompilationParameters>& cparams,
                   std::deque<std::string>& decl_filenames,
                   std::deque<std::string>& def_filenames)
  {
    //
    // generate explicit code for all recurrence relation that were not inlined
    //
    SafePtr<CodeContext> context(new CppCodeContext(cparams));
    ImplicitDimensions::set_default_dims(cparams);
    std::string prefix(cparams->source_directory());

    SafePtr<RRStack> rrstack = RRStack::Instance();
    RRStack::citer_type it = rrstack->begin();
    while (it != rrstack->end()) {
      SafePtr<RecurrenceRelation> rr = (*it).second.second;
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
      // next RR
      it = rrstack->begin();
    }
  }

};
