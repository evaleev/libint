
#include <rr.h>
#include <dg.h>
#include <strategy.h>

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
                                  std::ostream& decl, std::ostream& def)
{
  SafePtr<DirectedGraph> dg(new DirectedGraph);
  dg->append_target(rr_target());
  for(int c=0; c<num_children(); c++)
    dg->append_vertex(rr_child(c));
  // Always need to unroll integral sets
  SafePtr<Strategy> strat(new Strategy(1000000000));
  SafePtr<Tactic> ztactic(new ZeroNewVerticesTactic(dg));
  dg->apply(strat,ztactic);
}

///////////////

void
RRStack::add(const SafePtr<RRStack>& rrs)
{
  for(citer_type it=rrs->begin(); it != rrs->end(); it++) {
    find((*it).second);
  }
}

