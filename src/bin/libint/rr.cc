
#include <fstream>

#include <rr.h>
#include <dg.h>
#include <dg.templ.h>
#include <strategy.h>
#include <code.h>
#include <graph_registry.h>
#include <extract.h>
#include <algebra.h>
#include <context.h>
#include <integral.h>
#include <task.h>
#include <prefactors.h>
#include <singl_stack.timpl.h>

using namespace libint2;
using namespace libint2::prefactor;

RecurrenceRelation::RecurrenceRelation() :
  nflops_(0), expr_()
{
}

RecurrenceRelation::~RecurrenceRelation()
{
}

//
// If there is no generic equivalent, generate explicit code for this recurrence relation:
// 1) append target and children to a DirectedGraph dg
// 2) set their code symbols
// 3) apply IntSet_to_Ints
// 4) Apply RRs such that no additional vertices appear
// 5) call dg->generate_code()
//
void
RecurrenceRelation::generate_code(const SafePtr<CodeContext>& context,
                                  const SafePtr<ImplicitDimensions>& dims,
                                  const std::string& funcname,
                                  std::ostream& decl, std::ostream& def)
{
  //
  // Check if there is a generic equivalent that can be used
  //
  if (this->has_generic(context->cparams())) {
    generate_generic_code(context,dims,funcname,decl,def);
    return;
  }

  const SafePtr<DGVertex> target_vptr = rr_target();
  std::cout << "RecurrenceRelation::generate_code: target = " << target_vptr->label() << std::endl;

  const SafePtr<CompilationParameters>& cparams = context->cparams();
  SafePtr<DirectedGraph> dg = generate_graph_();

  // Intermediates in RR code are either are automatic variables or have to go on vstack
  dg->registry()->stack_name("inteval->vstack");
  // No need to return the targets via inteval's targets
  dg->registry()->return_targets(false);

  // check if CSE to be performed
  typedef IntegralSet<IncableBFSet> ISet;
  SafePtr<ISet> target = dynamic_pointer_cast<ISet,DGVertex>(target_vptr);
  if (target) {
    //
    // do CSE only if max_am <= cparams->max_am_opt()
    //
    const unsigned int np = target->num_part();
    unsigned int max_am = 0;
    // bra
    for(unsigned int p=0; p<np; p++) {
      const unsigned int nf = target->num_func_bra(p);
      for(unsigned int f=0; f<nf; f++) {
	// Assuming shells here
	const unsigned int am = target->bra(p,f).norm();
	using std::max;
	max_am = max(max_am,am);
      }
    }
    // ket
    for(unsigned int p=0; p<np; p++) {
      const unsigned int nf = target->num_func_ket(p);
      for(unsigned int f=0; f<nf; f++) {
	// Assuming shells here
	const unsigned int am = target->ket(p,f).norm();
	using std::max;
	max_am = max(max_am,am);
      }
    }
    const bool need_to_optimize = (max_am <= cparams->max_am_opt());
    dg->registry()->do_cse(need_to_optimize);
  }
  dg->registry()->condense_expr(condense_expr(1000000000,cparams->max_vector_length()>1));
  dg->registry()->ignore_missing_prereqs(true);  // assume all prerequisites are available -- if some are not, something is VERY broken

#if DEBUG
  {
    std::basic_ofstream<char> dotfile("graph_rr.strat.dot");
    dg->print_to_dot(false,dotfile);
  }
#endif

  // Assign symbols for the target and source integral sets
  SafePtr<CodeSymbols> symbols(new CodeSymbols);
  assign_symbols_(symbols);
  // Traverse the graph
  dg->optimize_rr_out(context);
  dg->traverse();
#if DEBUG
    {
      std::basic_ofstream<char> dotfile("graph_rr.expr.dot");
      dg->print_to_dot(false,dotfile);
    }
#endif
  // Generate code
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
  SafePtr<ImplicitDimensions> localdims = adapt_dims_(dims);
  dg->generate_code(context,memman,localdims,symbols,funcname,decl,def);

  // extract all external symbols -- these will be members of the evaluator structure
  SafePtr<ExtractExternSymbols> extractor(new ExtractExternSymbols);
  dg->foreach(*extractor);
  const ExtractExternSymbols::Symbols& externsymbols = extractor->symbols();

#if 0
  // print out the symbols
  std::cout << "Recovered symbols from DirectedGraph for " << label() << std::endl;
  typedef ExtractExternSymbols::Symbols::const_iterator citer;
  citer end = externsymbols.end();
  for(citer t=externsymbols.begin(); t!=end; ++t)
    std::cout << *t << std::endl;
#endif

#if DEBUG
    {
      std::basic_ofstream<char> dotfile("graph_rr.symb.dot");
      dg->print_to_dot(false,dotfile);
    }
#endif

  // get this RR InstanceID
  RRStack::InstanceID myid = RRStack::Instance()->find(EnableSafePtrFromThis<this_type>::SafePtr_from_this()).first;

  // For each task which requires this RR:
  // 1) update max stack size
  // 2) append external symbols from this RR to its list
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.plast();
  for(tciter t=taskmgr.first(); t!=tend; ++t) {
    const SafePtr<TaskExternSymbols> tsymbols = t->symbols();
    if (tsymbols->find(myid)) {
      // update max stack size
      t->params()->max_vector_stack_size(memman->max_memory_used());
      // add external symbols
      tsymbols->add(externsymbols);
    }
  }

  dg->reset();
}

namespace libint2 {
  // generate_generic_code reuses this function from dg.cc:
  extern std::string declare_function(const SafePtr<CodeContext>& context, const SafePtr<ImplicitDimensions>& dims,
                                      const SafePtr<CodeSymbols>& args, const std::string& tlabel, const std::string& function_descr,
                                      std::ostream& decl);
}

void
RecurrenceRelation::generate_generic_code(const SafePtr<CodeContext>& context,
                                          const SafePtr<ImplicitDimensions>& dims,
                                          const std::string& funcname,
                                          std::ostream& decl, std::ostream& def)
{
  const SafePtr<DGVertex> target_vptr = rr_target();
  std::cout << "RecurrenceRelation::generate_generic_code: target = " << target_vptr->label() << std::endl;

  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  const std::string tlabel = taskmgr.current().label();
  SafePtr<ImplicitDimensions> localdims = adapt_dims_(dims);
  // Assign symbols for the target and source integral sets
  SafePtr<CodeSymbols> symbols(new CodeSymbols);
  assign_symbols_(symbols);

  // declare function
  const std::string func_decl = declare_function(context,localdims,symbols,tlabel,funcname,decl);

  //
  // Generate function's definition
  //

  // include standard headers
  def << context->std_header();
  //         + generic code declaration
  def << "#include <"
      << this->generic_header()
      << ">" << endl;
  def << endl;

  // start the body ...
  def << context->code_prefix();
  def << func_decl << context->open_block() << endl;
  def << context->std_function_header();

  // ... fill the body
  def << this->generic_instance(context,symbols) << endl;

  // ... end the body
  def << context->close_block() << endl;
  def << context->code_postfix();
}

SafePtr<DirectedGraph>
RecurrenceRelation::generate_graph_()
{
  SafePtr<DirectedGraph> dg(new DirectedGraph);
  dg->append_target(rr_target());
  for(unsigned int c=0; c<num_children(); c++)
    dg->append_vertex(rr_child(c));
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- the number of integral sets = " << dg->num_vertices() << endl;
#endif
  SafePtr<Strategy> strat(new Strategy);
  SafePtr<Tactic> ntactic(new NullTactic);
  // Always need to unroll integral sets first
  dg->registry()->unroll_threshold(1000000000);
  dg->apply(strat,ntactic);
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- the number of integral sets + integrals = " << dg->num_vertices() << endl;
#endif
  // Mark children sets and their descendants to not compute
  for(unsigned int c=0; c<num_children(); c++)
    dg->apply_at<&DGVertex::not_need_to_compute>(rr_child(c));
  // Apply recurrence relations using existing vertices on the graph (i.e.
  // such that no new vertices appear)
  SafePtr<Tactic> ztactic(new ZeroNewVerticesTactic(dg));
  dg->apply(strat,ztactic);
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- should be same as previous = " << dg->num_vertices() << endl;
#endif

  return dg;
}

void
RecurrenceRelation::assign_symbols_(SafePtr<CodeSymbols>& symbols)
{
  // Set symbols on the target and children sets
  rr_target()->set_symbol("target");
  symbols->append_symbol("target");
  for(unsigned int c=0; c<num_children(); c++) {
    ostringstream oss;
    oss << "src" << c;
    string symb = oss.str();
    rr_child(c)->set_symbol(symb);
    symbols->append_symbol(symb);
  }
}

SafePtr<ImplicitDimensions>
RecurrenceRelation::adapt_dims_(const SafePtr<ImplicitDimensions>& dims) const
{
  return dims;
}

const std::string&
RecurrenceRelation::label() const {
  if (label_.empty())
    label_ = generate_label();
  return label_;
}

std::string
RecurrenceRelation::description() const
{
  const std::string descr = label();
  return descr;
}

void
RecurrenceRelation::add_expr(const SafePtr<ExprType>& expr, int minus)
{
  if (expr_ == 0) {
    if (minus != -1) {
      expr_ = expr;
    }
    else {
      using libint2::prefactor::Scalar;
      SafePtr<ExprType> negative(new ExprType(ExprType::OperatorTypes::Times,expr,Scalar(-1.0)));
      expr_ = negative;
      ++nflops_;
    }
  }
  else {
    if (minus != -1) {
      SafePtr<ExprType> sum(new ExprType(ExprType::OperatorTypes::Plus,expr_,expr));
      expr_ = sum;
      ++nflops_;
    }
    else {
      SafePtr<ExprType> sum(new ExprType(ExprType::OperatorTypes::Minus,expr_,expr));
      expr_ = sum;
      ++nflops_;
    }
  }
}


bool
RecurrenceRelation::invariant_type() const {
  // By default, recurrence relations do not change the type of the functions, i.e. VRR applied to an integral over shells will produce integrals over shells
  return true;
}

std::string
RecurrenceRelation::spfunction_call(const SafePtr<CodeContext>& context, const SafePtr<ImplicitDimensions>& dims) const
{
  ostringstream os;
  os << context->label_to_name(label_to_funcname(context->cparams()->api_prefix() + label()))
    // First argument is the library object
     << "(inteval, "
    // Second is the target
     << context->value_to_pointer(rr_target()->symbol());
  // then come children
  const unsigned int nchildren = num_children();
  for(unsigned int c=0; c<nchildren; c++) {
    os << ", " << context->value_to_pointer(rr_child(c)->symbol());
  }
  os << ")" << context->end_of_stat() << endl;
  return os.str();
}

bool
RecurrenceRelation::has_generic(const SafePtr<CompilationParameters>& cparams) const {
  return false;
}

std::string
RecurrenceRelation::generic_header() const {
  throw std::logic_error("RecurrenceRelation::generic_header() -- should not be called! Check if DerivedRecurrenceRelation::generic_header() is implemented");
}

std::string
RecurrenceRelation::generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const {
  throw std::logic_error("RecurrenceRelation::generic_instance() -- should not be called! Check if DerivedRecurrenceRelation::generic_instance() is implemented");
}

namespace libint2 { namespace algebra {
  /// these operators are extremely useful to write compact expressions
  SafePtr<RecurrenceRelation::ExprType> operator+(const SafePtr<DGVertex>& A,
                                                  const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    return SafePtr<Oper>(new Oper(Oper::OperatorTypes::Plus,A,B));
  }
  SafePtr<RecurrenceRelation::ExprType> operator-(const SafePtr<DGVertex>& A,
                                                  const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    return SafePtr<Oper>(new Oper(Oper::OperatorTypes::Minus,A,B));
  }
  SafePtr<RecurrenceRelation::ExprType> operator*(const SafePtr<DGVertex>& A,
                                                  const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    return SafePtr<Oper>(new Oper(Oper::OperatorTypes::Times,A,B));
  }
  SafePtr<RecurrenceRelation::ExprType> operator/(const SafePtr<DGVertex>& A,
                                                  const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    return SafePtr<Oper>(new Oper(Oper::OperatorTypes::Divide,A,B));
  }
  const SafePtr<RecurrenceRelation::ExprType>& operator+=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                          const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    const SafePtr<Oper>& Sum = A + B;
    A = Sum;
    return A;
  }
  const SafePtr<RecurrenceRelation::ExprType>& operator-=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                          const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    const SafePtr<Oper>& Diff = A - B;
    A = Diff;
    return A;
  }
  const SafePtr<RecurrenceRelation::ExprType>& operator*=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                          const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    const SafePtr<Oper>& Product = A * B;
    A = Product;
    return A;
  }
  const SafePtr<RecurrenceRelation::ExprType>& operator/=(SafePtr<RecurrenceRelation::ExprType>& A,
                                                          const SafePtr<DGVertex>& B) {
    typedef RecurrenceRelation::ExprType Oper;
    const SafePtr<Oper>& Quotient = A / B;
    A = Quotient;
    return A;
  }
} } // namespace libint2::algebra

///////////////



///////////////
