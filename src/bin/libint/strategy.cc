
#define USE_ITR 0
#define USE_HRR 1

#include <vector>
#include <algorithm>
#include <strategy.h>
#include <dg.h>
#include <hrr.h>
#include <vrr_11_twoprep_11.h>
#include <itr_11_twoprep_11.h>
#include <vrr_11_r12kg12_11.h>
#include <comp_11_tig12_11.h>
#include <dummyintegral.h>

using namespace std;
using namespace libint2;

SafePtr<RecurrenceRelation>
Strategy::optimal_rr(const SafePtr<DirectedGraph>& graph,
                     const SafePtr<DGVertex>& integral,
                     const SafePtr<Tactic>& tactic)
{
  //
  // We must first determine the type of the integral
  //
  {
    SafePtr<TwoPRep_11_11_sq> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_sq,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_sq(graph,eri_ptr,tactic);
  }
  {
    SafePtr<TwoPRep_11_11_int> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_int,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_int(graph,eri_ptr,tactic);
  }
  {
    typedef R12kG12_11_11_base<CGShell> base_type;
    SafePtr<base_type> bptr = dynamic_pointer_cast<base_type,DGVertex>(integral);
    if (bptr != 0) {
      int k = R12kG12_11_11_Util::k<CGShell>(bptr);
      switch (k) {
        case 2:
          return optimal_rr_R12kG121111_sq<2>(graph,bptr,tactic);
        case 0:
          return optimal_rr_R12kG121111_sq<0>(graph,bptr,tactic);
        case -1:
          return optimal_rr_R12kG121111_sq<-1>(graph,bptr,tactic);
        default:
          throw logic_error("Strategy::optimal_rr() unable to determine K for R12kG12_11_11<CGShell,K> class");
      };
    }
  }
  {
    typedef R12kG12_11_11_base<CGF> base_type;
    SafePtr<base_type> bptr = dynamic_pointer_cast<base_type,DGVertex>(integral);
    if (bptr != 0) {
      int k = R12kG12_11_11_Util::k<CGF>(bptr);
      switch (k) {
        case 2:
          return optimal_rr_R12kG121111_int<2>(graph,bptr,tactic);
        case 0:
          return optimal_rr_R12kG121111_int<0>(graph,bptr,tactic);
        case -1:
          return optimal_rr_R12kG121111_int<-1>(graph,bptr,tactic);
        default:
          throw logic_error("Strategy::optimal_rr() unable to determine K for R12kG12_11_11<CGF,K> class");
      };
    }
  }
  {
    typedef TiG12_11_11_base<CGShell> base_type;
    SafePtr<base_type> bptr = dynamic_pointer_cast<base_type,DGVertex>(integral);
    if (bptr != 0) {
      int k = TiG12_11_11_Util::i<CGShell>(bptr);
      switch (k) {
        case 0:
          return optimal_rr_TiG121111_sq<0>(graph,bptr,tactic);
        case 1:
          return optimal_rr_TiG121111_sq<1>(graph,bptr,tactic);
        default:
          throw logic_error("Strategy::optimal_rr() unable to determine K for TiG12_11_11<CGShell,K> class");
      };
    }
  }
  {
    typedef TiG12_11_11_base<CGF> base_type;
    SafePtr<base_type> bptr = dynamic_pointer_cast<base_type,DGVertex>(integral);
    if (bptr != 0) {
      int k = TiG12_11_11_Util::i<CGF>(bptr);
      switch (k) {
        case 0:
          return optimal_rr_TiG121111_int<0>(graph,bptr,tactic);
        case 1:
          return optimal_rr_TiG121111_int<1>(graph,bptr,tactic);
        default:
          throw logic_error("Strategy::optimal_rr() unable to determine K for TiG12_11_11<CGF,K> class");
      };
    }
  }

  // Type insensitive RR can be applied to almost any integral
  {
    typedef DummySymmIntegral_11_11_sq int_type;
    SafePtr<int_type> int_ptr = dynamic_pointer_cast<int_type,DGVertex>(integral);
    if (int_ptr != 0)
      return optimal_rr_Dummy1111_sq(graph,int_ptr,tactic);
  }
  {
    typedef DummySymmIntegral_11_11_int int_type;
    SafePtr<int_type> int_ptr = dynamic_pointer_cast<int_type,DGVertex>(integral);
    if (int_ptr != 0)
      return optimal_rr_Dummy1111_int(graph,int_ptr,tactic);
  }

  // Don't know how to apply any RR
  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_sq(const SafePtr<DirectedGraph>& graph,
                                    const SafePtr<TwoPRep_11_11_sq>& integral,
                                    const SafePtr<Tactic>& tactic)
{
  //
  // This is a basic strategy for computing integral
  // 1) first see if should convert the set to infividual integrals
  // 2) if possible apply HRR
  // 3) else apply VRR
  //
  if (integral->size() <= max_size_to_unroll_)
    return unroll_intset<TwoPRep_11_11_sq>(integral);

#if USE_HRR
  {
    typedef HRR_ab_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    typedef HRR_cd_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
#endif
  
#if USE_ITR
  {
    typedef ITR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
#endif
  
  {
    typedef VRR_a_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
  
  {
    typedef VRR_b_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
  
  {
    typedef VRR_c_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
  
  {
    typedef VRR_d_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
  
  return SafePtr<RecurrenceRelation>();
}



// Generate all possible recurrence relations and then
// use a Tactic object to decide which to use

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr<TwoPRep_11_11_int>& integral,
                                     const SafePtr<Tactic>& tactic)
{
  vector<RR> rrstack;  // stack of all recurrence relations
  
#if USE_HRR
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_ab_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_cd_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
#endif

#if USE_ITR
  // shift from A to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef ITR_11_TwoPRep_11<TwoPRep_11_11,CGF,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
#endif

  // decrease A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  // decrease B
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,0,InKet> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  // Decrease C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,1,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // Decrease D
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,1,InKet> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  return tactic->optimal_rr(rrstack);
}


SafePtr<RecurrenceRelation>
Strategy::optimal_rr_Dummy1111_sq(const SafePtr<DirectedGraph>& graph,
				  const SafePtr<DummySymmIntegral_11_11_sq>& integral,
				  const SafePtr<Tactic>& tactic)
{
  if (integral->size() <= max_size_to_unroll_)
    return unroll_intset<DummySymmIntegral_11_11_sq>(integral);

  {
    typedef HRR_ab_11_Dummy_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    typedef HRR_cd_11_Dummy_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_Dummy1111_int(const SafePtr<DirectedGraph>& graph,
				   const SafePtr<DummySymmIntegral_11_11_int>& integral,
				   const SafePtr<Tactic>& tactic)
{
  vector<RR> rrstack;  // stack of all recurrence relations
  
#if USE_HRR
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_ab_11_Dummy_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_cd_11_Dummy_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
#endif

  return tactic->optimal_rr(rrstack);
}


