
#include <smart_ptr.h>
#include <tactic.h>
#include <intset_to_ints.h>
#include <vrr_11_r12kg12_11.h>
#include <comp_11_tig12_11.h>
#include <hrr.h>
#include <global_macros.h>

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
    template <int K>
    RR optimal_rr_R12kG121111_sq(const SafePtr<DirectedGraph>& graph,
                                 const SafePtr< R12kG12_11_11_base<CGShell> >& integral,
                                 const SafePtr<Tactic>& tactic);
    template <int K>
    RR optimal_rr_R12kG121111_int(const SafePtr<DirectedGraph>& graph,
                                  const SafePtr< R12kG12_11_11_base<CGF> >& integral,
                                  const SafePtr<Tactic>& tactic);
    template <int K>
    RR optimal_rr_TiG121111_sq(const SafePtr<DirectedGraph>& graph,
                               const SafePtr< TiG12_11_11_base<CGShell> >& integral,
                               const SafePtr<Tactic>& tactic);
    template <int K>
    RR optimal_rr_TiG121111_int(const SafePtr<DirectedGraph>& graph,
                                const SafePtr< TiG12_11_11_base<CGF> >& integral,
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

  };
  
  template <int K>
    SafePtr<RecurrenceRelation>
    Strategy::optimal_rr_R12kG121111_sq(const SafePtr<DirectedGraph>& graph,
                                        const SafePtr< R12kG12_11_11_base<CGShell> >& bptr,
                                        const SafePtr<Tactic>& tactic)
    {
      typedef R12kG12_11_11<CGShell,K> inttype;
      const SafePtr<inttype> integral = dynamic_pointer_cast<inttype,R12kG12_11_11_base<CGShell> >(bptr);
      
      //
      // This is a basic strategy for computing integral
      // 1) first see if should convert the SafePtr<RecurrenceRelation>
      //    set to infividual integrals
      // 2) if possible apply HRR
      // 3) else apply VRR
      //
      const unsigned int size = integral->size();
      if (size == 1 || (size <= max_size_to_unroll_ && graph->registry()->can_unroll()))
        return unroll_intset<inttype>(integral);
      
      {
        // AB HRR
        typedef HRR<R12kG12_11_11<CGShell,K>,CGShell,0,InBra,0,InKet,0> rr_type;
        SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
        if (rr_ptr->num_children())
          return rr_cast(rr_ptr);
      }
      
      {
        // CD HRR
        typedef HRR<R12kG12_11_11<CGShell,K>,CGShell,1,InBra,0,InKet,0> rr_type;
        SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
        if (rr_ptr->num_children())
          return rr_cast(rr_ptr);
      }
      
      {
        // A VRR
        typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,K,0,InBra> rr_type;
        SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
        if (rr_ptr->num_children())
          return rr_cast(rr_ptr);
      }
      
      {
        // C VRR
        typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,K,1,InBra> rr_type;
        SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
        if (rr_ptr->num_children())
          return rr_cast(rr_ptr);
      }
      
      return SafePtr<RecurrenceRelation>();
    }

  
  // Generate all possible recurrence relations and then
  // use a Tactic object to decide which to use
  template <int K>
    SafePtr<RecurrenceRelation>
    Strategy::optimal_rr_R12kG121111_int(const SafePtr<DirectedGraph>& graph,
                                         const SafePtr< R12kG12_11_11_base<CGF> >& bptr,
                                         const SafePtr<Tactic>& tactic)
    {
      typedef R12kG12_11_11<CGF,K> inttype;
      const SafePtr<inttype> integral = dynamic_pointer_cast<inttype,R12kG12_11_11_base<CGF> >(bptr);

      vector<RR> rrstack;  // stack of all recurrence relations
      
      // shift from B to A
      for(int xyz = 2; xyz >= 0; xyz--) {
        typedef HRR<R12kG12_11_11<CGF,K>,CGF,0,InBra,0,InKet,0> rr_type;
        SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
        if (rr_ptr->num_children())
          rrstack.push_back(rr_cast(rr_ptr));
      }
      
      // shift from D to C
      for(int xyz = 2; xyz >= 0; xyz--) {
        typedef HRR<R12kG12_11_11<CGF,K>,CGF,1,InBra,0,InKet,0> rr_type;
        SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
        if (rr_ptr->num_children())
          rrstack.push_back(rr_cast(rr_ptr));
      }

      // only apply VRR is AM on B and D is zero
#if USE_BRAKET_H
      if (integral->ket(0,0).zero() && integral->ket(1,0).zero()) {
#else
      if (integral->ket(0,0)->zero() && integral->ket(1,0)->zero()) {
#endif
        // decrease A
        for(int xyz = 2; xyz >= 0; xyz--) {
          typedef VRR_11_R12kG12_11<R12kG12_11_11,CGF,K,0,InBra> rr_type;
          SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
          if (rr_ptr->num_children())
            rrstack.push_back(rr_cast(rr_ptr));
        }
        
        // Decrease C
        for(int xyz = 2; xyz >= 0; xyz--) {
          typedef VRR_11_R12kG12_11<R12kG12_11_11,CGF,K,1,InBra> rr_type;
          SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
          if (rr_ptr->num_children())
            rrstack.push_back(rr_cast(rr_ptr));
        }
      }
      
      return tactic->optimal_rr(rrstack);
    }

  
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
      
      typedef CR_11_TiG12_11<TiG12_11_11,CGShell,K> rr_type;
      SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
      if (rr_ptr->num_children())
        return rr_cast(rr_ptr);
      
      throw std::logic_error("Strategy::optimal_rr_TiG121111_sq() -- did not find a way to compute");
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

      typedef CR_11_TiG12_11<TiG12_11_11,CGF,K> rr_type;
      SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
      if (rr_ptr->num_children())
        rrstack.push_back(rr_cast(rr_ptr));
      else
        throw std::logic_error("Strategy::optimal_rr_TiG121111_int() -- did not find a way to compute");
      
      return tactic->optimal_rr(rrstack);
    }
  
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
