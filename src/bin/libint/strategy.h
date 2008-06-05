
#include <smart_ptr.h>
#include <tactic.h>
#include <intset_to_ints.h>
#include <vrr_11_twoprep_11.h>
#include <vrr_11_r12kg12_11.h>
#include <comp_11_tig12_11.h>
#include <comp_11_r1dotr1g12_11.h>
#include <comp_11_r2dotr2g12_11.h>
#include <comp_11_r1dotr2g12_11.h>
#include <hrr.h>
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

  protected:

#if 0
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
    virtual RR optimal_rr_R12kG121111_sq(const SafePtr<DirectedGraph>& graph,
                                         const SafePtr< R12kG12_11_11_sq >& integral,
                                         const SafePtr<Tactic>& tactic);
    virtual RR optimal_rr_R12kG121111_int(const SafePtr<DirectedGraph>& graph,
                                          const SafePtr< R12kG12_11_11_int >& integral,
                                          const SafePtr<Tactic>& tactic);
#if 0
    template <int K>
    RR optimal_rr_TiG121111_sq(const SafePtr<DirectedGraph>& graph,
                               const SafePtr< TiG12_11_11_base<CGShell> >& integral,
                               const SafePtr<Tactic>& tactic);
    template <int K>
    RR optimal_rr_TiG121111_int(const SafePtr<DirectedGraph>& graph,
                                const SafePtr< TiG12_11_11_base<CGF> >& integral,
                                const SafePtr<Tactic>& tactic);
    RR optimal_rr_R1dotR1G121111_sq(const SafePtr<DirectedGraph>& graph,
				    const SafePtr<R1dotR1G12_11_11_sq>& integral,
				    const SafePtr<Tactic>& tactic);
    RR optimal_rr_R1dotR1G121111_int(const SafePtr<DirectedGraph>& graph,
				     const SafePtr<R1dotR1G12_11_11_int>& integral,
				     const SafePtr<Tactic>& tactic);
    RR optimal_rr_R2dotR2G121111_sq(const SafePtr<DirectedGraph>& graph,
				    const SafePtr<R2dotR2G12_11_11_sq>& integral,
				    const SafePtr<Tactic>& tactic);
    RR optimal_rr_R2dotR2G121111_int(const SafePtr<DirectedGraph>& graph,
				     const SafePtr<R2dotR2G12_11_11_int>& integral,
				     const SafePtr<Tactic>& tactic);
#endif
    RR optimal_rr_R1dotR2G121111_sq(const SafePtr<DirectedGraph>& graph,
				    const SafePtr<R1dotR2G12_11_11_sq>& integral,
				    const SafePtr<Tactic>& tactic);
    RR optimal_rr_R1dotR2G121111_int(const SafePtr<DirectedGraph>& graph,
				     const SafePtr<R1dotR2G12_11_11_int>& integral,
				     const SafePtr<Tactic>& tactic);
    
    virtual RR optimal_rr_Dummy1111_sq(const SafePtr<DirectedGraph>& graph,
				       const SafePtr<DummySymmIntegral_11_11_sq>& integral,
				       const SafePtr<Tactic>& tactic);
    virtual RR optimal_rr_Dummy1111_int(const SafePtr<DirectedGraph>& graph,
					const SafePtr<DummySymmIntegral_11_11_int>& integral,
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
#endif
  };
  

#if 0  
  template <int K>
    SafePtr<RecurrenceRelation>
    Strategy::optimal_rr_TiG121111_sq(const SafePtr<DirectedGraph>& graph,
                                      const SafePtr< TiG12_11_11_base<CGShell> >& bptr,
                                      const SafePtr<Tactic>& tactic)
    {
      typedef TiG12_11_11<CGShell,K> inttype;
      const SafePtr<inttype> integral = dynamic_pointer_cast<inttype,TiG12_11_11_base<CGShell> >(bptr);
      
      //
      // This is a basic strategy for computing integral
      // 1) first see if should convert the SafePtr<RecurrenceRelation>
      //    set to infividual integrals
      // 2) otherwise compute it
      //
      const unsigned int size = integral->size();
      if (size == 1 || (size <= max_size_to_unroll_ && graph->registry()->can_unroll()))
        return unroll_intset<inttype>(integral);

#if USE_HRR_FOR_TiG12
      if (K == 1) {
	// shift from B to A
	typedef HRR<inttype,typename inttype::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
	SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
	if (rr_ptr->num_children())
	  return rr_cast(rr_ptr);
      }

      if (K == 0) {
	// shift from D to C
	typedef HRR<inttype,typename inttype::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
	SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
	if (rr_ptr->num_children())
	  return rr_cast(rr_ptr);
      }
#endif

      {
	typedef CR_11_TiG12_11<TiG12_11_11,CGShell,K> rr_type;
	SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
	if (rr_ptr->num_children())
	  return rr_cast(rr_ptr);
      }

      return SafePtr<RecurrenceRelation>();
    }

  
  // Generate all possible recurrence relations and then
  // use a Tactic object to decide which to use
  template <int K>
    SafePtr<RecurrenceRelation>
    Strategy::optimal_rr_TiG121111_int(const SafePtr<DirectedGraph>& graph,
                                       const SafePtr< TiG12_11_11_base<CGF> >& bptr,
                                       const SafePtr<Tactic>& tactic)
    {
      typedef TiG12_11_11<CGF,K> inttype;
      const SafePtr<inttype> integral = dynamic_pointer_cast<inttype,TiG12_11_11_base<CGF> >(bptr);
      vector<RR> rrstack;  // stack of all recurrence relations

#if USE_HRR_FOR_TiG12
      if (K == 1)
	// shift from B to A
	for(int xyz = 2; xyz >= 0; xyz--) {
	  typedef HRR<inttype,typename inttype::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
	  SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
	  if (rr_ptr->num_children())
	    rrstack.push_back(rr_cast(rr_ptr));
	}

      if (K == 0)
	// shift from D to C
	for(int xyz = 2; xyz >= 0; xyz--) {
	  typedef HRR<inttype,typename inttype::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
	  SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
	  if (rr_ptr->num_children())
	    rrstack.push_back(rr_cast(rr_ptr));
	}
#endif

      typedef CR_11_TiG12_11<TiG12_11_11,CGF,K> rr_type;
      SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
      if (rr_ptr->num_children())
        rrstack.push_back(rr_cast(rr_ptr));
      
      return tactic->optimal_rr(rrstack);
    }
#endif
  
};

#endif
