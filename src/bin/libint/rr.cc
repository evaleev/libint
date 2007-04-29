
#include <rr.h>
#include <dg.h>
#include <dg.templ.h>
#include <strategy.h>
#include <code.h>
#include <graph_registry.h>
#include <extract.h>

using namespace libint2;

RecurrenceRelation::RecurrenceRelation()
{
}

RecurrenceRelation::~RecurrenceRelation()
{
}

//
// Generate code for this recurrence relation:
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
  const SafePtr<CompilationParameters>& cparams = context->cparams();
  SafePtr<DirectedGraph> dg = generate_graph_();
  
  // Intermediates in RR code are either are automatic variables or have to go on vstack
  dg->registry()->stack_name("inteval->vstack");
  // No need to return the targets via inteval's targets
  dg->registry()->return_targets(false);

  // check if CSE to be performed
  typedef IntegralSet<IncableBFSet> ISet;
  SafePtr<DGVertex> target_vptr = rr_target();
  std::cout << "RecurrenceRelation::generate_code: target = " << target_vptr->label() << std::endl;
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

  // Assign symbols for the target and source integral sets
  SafePtr<CodeSymbols> symbols(new CodeSymbols);
  assign_symbols_(symbols);
  // Traverse the graph
  dg->optimize_rr_out();
  dg->traverse();
#if DEBUG
  dg->debug_print_traversal(std::cout);
  cout << "The number of vertices = " << dg->num_vertices() << endl;
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

  // get this RR InstanceID
  RRStack::InstanceID myid = RRStack::Instance()->find(EnableSafePtrFromThis<this_type>::SafePtr_from_this()).first;

  // For each task which requires this RR:
  // 1) update max stack size
  // 2) append external symbols from this RR to its list
  LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
  typedef LibraryTaskManager::TasksCIter tciter;
  const tciter tend = taskmgr.last();
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

SafePtr<DirectedGraph>
RecurrenceRelation::generate_graph_()
{
  SafePtr<DirectedGraph> dg(new DirectedGraph);
  dg->append_target(rr_target());
  for(int c=0; c<num_children(); c++)
    dg->append_vertex(rr_child(c));
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- the number of integral sets = " << dg->num_vertices() << endl;
#endif
  // Always need to unroll integral sets
  SafePtr<Strategy> strat(new Strategy(1000000000));
  SafePtr<Tactic> ntactic(new NullTactic);
  dg->apply(strat,ntactic);
#if DEBUG
  cout << "RecurrenceRelation::generate_code -- the number of integral sets + integrals = " << dg->num_vertices() << endl;
#endif
  // Mark children sets and their descendants to not compute
  for(int c=0; c<num_children(); c++)
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
  for(int c=0; c<num_children(); c++) {
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
RecurrenceRelation::description() const
{
  return label();
}

///////////////

SafePtr<RRStack>
RRStack::rrstack_;

SafePtr<RRStack>&
RRStack::Instance() {
  if (!rrstack_) {
    SafePtr<RRStack> tmpstack(new RRStack);
    rrstack_ = tmpstack;
  }
  return rrstack_;
}

void
RRStack::add(const SafePtr<RRStack>& rrs)
{
  for(citer_type it=rrs->begin(); it != rrs->end(); it++) {
    find((*it).second.second);
  }
}

void
RRStack::remove(const data_type& rr)
{
  parent_type::remove(rr);
}

