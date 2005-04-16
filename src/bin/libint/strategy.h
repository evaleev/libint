
#include <smart_ptr.h>
#include <tactic.h>
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

    static const unsigned int default_max_size_to_unroll = 1;

  public:
    typedef SafePtr<RecurrenceRelation> RR;
    Strategy(unsigned int max_size_to_unroll = default_max_size_to_unroll) :
      max_size_to_unroll_(max_size_to_unroll) {
        if (max_size_to_unroll == 0)
          throw std::runtime_error("Strategy::Strategy() -- max_size_to_unroll must be >= 1");
      }
    virtual ~Strategy() {}

    /// Returns the optimal recurrence relation for integral
    RR optimal_rr(const SafePtr<DirectedGraph>& graph,
                  const SafePtr<DGVertex>& integral,
                  const SafePtr<Tactic>& tactic);

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

    virtual RR optimal_rr_twoprep1111_sq(const SafePtr<DirectedGraph>& graph,
                                         const SafePtr< TwoPRep_11_11_sq >& integral,
                                         const SafePtr<Tactic>& tactic);
    virtual RR optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                          const SafePtr< TwoPRep_11_11_int >& integral,
                                          const SafePtr<Tactic>& tactic);

    /// Casts a smart pointer to rr of type RRType to a smart pointer to RecurrenceRelation. Utility function.
    template <class RRType>
      RR rr_cast(const SafePtr<RRType>& rr) const
      {
        RR rr_cast_ptr = static_pointer_cast<RecurrenceRelation,RRType>(rr);
        return rr_cast_ptr;
      }

    private:
    unsigned int max_size_to_unroll_;

  };
  
  
  //
  // IntSetRRStrategy describes how to compute a set of integrals from other,
  // precomputed sets of integrals. For example, if I need to figure out how
  // to compute a quartet of electron repulsion integrals from other quartets
  // of ERIs using a given RR (such as Obara-Saika).
  //
  class IntSetRRStrategy: public Strategy {
    public:
    typedef Strategy::RR RR;

    // actual computations are done at the integral level,
    // hence integral sets to be unrolled
    IntSetRRStrategy() : Strategy(1000000000) {}
    ~IntSetRRStrategy() {}
    
    private:
    // need to overload
    RR optimal_rr_twoprep1111_sq(const SafePtr<DirectedGraph>& graph,
                                 const SafePtr< TwoPRep_11_11_sq >& integral,
                                 const SafePtr<Tactic>& tactic);
    // need to overload
    RR optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                  const SafePtr< TwoPRep_11_11_sq >& integral,
                                  const SafePtr<Tactic>& tactic);
    
  };
  
};

#endif
