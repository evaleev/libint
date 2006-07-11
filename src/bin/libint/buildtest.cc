
#include <iostream>
#include <fstream>
#include <rr.h>
#include <context.h>
#include <dims.h>

using namespace libint2;

namespace libint2 {

void
generate_rr_code(std::ostream& os, const SafePtr<CompilationParameters>& cparams)
{
  //
  // generate explicit code for all recurrence relation that were not inlined
  //
  SafePtr<CodeContext> context(new CppCodeContext(cparams));
  ImplicitDimensions::set_default_dims(cparams);
  std::string prefix(cparams->source_directory());

  SafePtr<RRStack> rrstack = RRStack::Instance();
  for(RRStack::citer_type it = rrstack->begin(); it!=rrstack->end(); it++) {
    SafePtr<RecurrenceRelation> rr = (*it).second.second;
    std::string rrlabel = cparams->api_prefix() + rr->label();
    os << " generating code for " << context->label_to_name(rrlabel) << " target=" << rr->rr_target()->label() << endl;
    
    std::string decl_filename(prefix + context->label_to_name(rrlabel));  decl_filename += ".h";
    std::string def_filename(prefix + context->label_to_name(rrlabel));  def_filename += ".cc";
    std::basic_ofstream<char> declfile(decl_filename.c_str());
    std::basic_ofstream<char> deffile(def_filename.c_str());
    
    rr->generate_code(context,ImplicitDimensions::default_dims(),rrlabel,declfile,deffile);
    
    declfile.close();
    deffile.close();
  }
}

};
