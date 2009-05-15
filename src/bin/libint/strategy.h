
#include <smart_ptr.h>
#include <tactic.h>
#include <global_macros.h>
#include <graph_registry.h>

#ifndef _libint2_src_bin_libint_strategy_h_
#define _libint2_src_bin_libint_strategy_h_

#define USE_HRR_FOR_TiG12 1

using namespace std;


namespace libint2 {

  class DGVertex;
  class DirectedGraph;

  /**
  Strategy specifies how to apply recurrence relations.
  */
  class Strategy {

  public:
    typedef SafePtr<RecurrenceRelation> RR;
    Strategy() {}
    ~Strategy() {}

    /// Returns the optimal recurrence relation for integral
    RR optimal_rr(const SafePtr<DirectedGraph>& graph,
                  const SafePtr<DGVertex>& integral,
                  const SafePtr<Tactic>& tactic);
  };

};

#endif
