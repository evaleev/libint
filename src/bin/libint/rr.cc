
#include <rr.h>
#include <dg.h>
#include <dg.templ.h>
#include <strategy.h>
#include <code.h>

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
                                  std::ostream& decl, std::ostream& def)
{
  SafePtr<DirectedGraph> dg = generate_graph_();
  // Assign symbols for the target and source integral sets
  SafePtr<CodeSymbols> symbols(new CodeSymbols);
  assign_symbols_(symbols);
  // Traverse the graph
  dg->optimize_rr_out();
  dg->traverse();
#if DEBUG
  dg->debug_print_traversal(std::cout);
#endif
  cout << "The number of vertices = " << dg->num_vertices() << endl;
  // Generate code
  SafePtr<MemoryManager> memman(new WorstFitMemoryManager());
  SafePtr<ImplicitDimensions> localdims = adapt_dims_(dims);
  dg->generate_code(context,memman,localdims,symbols,label(),decl,def);
  dg->reset();
}

SafePtr<DirectedGraph>
RecurrenceRelation::generate_graph_()
{
  SafePtr<DirectedGraph> dg(new DirectedGraph);
  dg->append_target(rr_target());
  for(int c=0; c<num_children(); c++)
    dg->append_vertex(rr_child(c));
  cout << "RecurrenceRelation::generate_code -- the number of integral sets = " << dg->num_vertices() << endl;
  // Always need to unroll integral sets
  SafePtr<Strategy> strat(new Strategy(1000000000));
  SafePtr<Tactic> ntactic(new NullTactic);
  dg->apply(strat,ntactic);
  cout << "RecurrenceRelation::generate_code -- the number of integral sets + integrals = " << dg->num_vertices() << endl;
  // Mark children sets and their descendants to not compute
  for(int c=0; c<num_children(); c++)
    dg->apply_at<&DGVertex::not_need_to_compute>(rr_child(c));
  // Apply recurrence relations using existing vertices on the graph (i.e.
  // such that no new vertices appear)
  SafePtr<Tactic> ztactic(new ZeroNewVerticesTactic(dg));
  dg->apply(strat,ztactic);
  cout << "RecurrenceRelation::generate_code -- should be same as previous = " << dg->num_vertices() << endl;
  
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

SafePtr<RRStack>
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

