
#include <smart_ptr.h>
#include <intset_to_ints.h>

#ifndef _libint2_src_bin_libint_strategy_h_
#define _libint2_src_bin_libint_strategy_h_

using namespace std;


namespace libint2 {

  class DGVertex;
  class DirectedGraph;
  
  /**
  Strategy specifies how to apply recurrence relations.
  */
  class Strategy {

  public:
    static const unsigned int max_size_to_unroll = 1111;
    typedef SafePtr<RecurrenceRelation> RR;
    Strategy() {}
    virtual ~Strategy() {}

    /// Returns the optimal recurrence relation for integral
    virtual RR optimal_rr(const SafePtr<DirectedGraph>& graph, const SafePtr<DGVertex>& integral);

  protected:

    /// Checks if need to unroll this integral set to individual integrals
    template <class I>
      RR unroll_intset(const SafePtr<I>& integral)
      {
        typedef IntegralSet_to_Integrals<I> ISet2I;
        SafePtr<ISet2I> rr(new ISet2I(integral));
        RR rr_cast = dynamic_pointer_cast<RecurrenceRelation,ISet2I>(rr);
        return rr_cast;
      }

    virtual RR optimal_rr_twoprep1111_sq(const SafePtr<DirectedGraph>& graph, const SafePtr< TwoPRep_11_11_sq >& integral);
    virtual RR optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph, const SafePtr< TwoPRep_11_11_int >& integral);

  };
    
};

#endif
