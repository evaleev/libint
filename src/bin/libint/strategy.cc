
#define USE_TACTIC 1

#include <vector>
#include <algorithm>
#include <strategy.h>
#include <dg.h>
#include <hrr.h>
#include <vrr_11_twoprep_11.h>

using namespace std;
using namespace libint2;

SafePtr<RecurrenceRelation>
Strategy::optimal_rr(const SafePtr<DirectedGraph>& graph,
                     const SafePtr<DGVertex>& integral,
                     const SafePtr<Tactic>& tactic)
{

  // We must first determine the type of the integral
  {
    SafePtr<TwoPRep_11_11_sq> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_sq,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_sq(graph,eri_ptr,tactic);
  }

  // We must first determine the type of the integral
  {
    SafePtr<TwoPRep_11_11_int> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_int,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_int(graph,eri_ptr,tactic);
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

  {
    typedef HRR_ab_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,0));
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    typedef HRR_cd_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,0));
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    typedef VRR_a_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,0));
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    typedef VRR_c_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,0));
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}


#if !USE_TACTIC

// This approach blindly seeks the first possible method to apply VRR

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr<TwoPRep_11_11_int>& integral,
                                     const SafePtr<Tactic>& tactic)
{
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_ab_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,xyz));
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_cd_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,xyz));
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  // decrease A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,0,InBra> vrr_type;
    SafePtr<vrr_type> vrr_ptr(new vrr_type(integral,xyz));
    if (vrr_ptr->num_children())
      return rr_cast(vrr_ptr);
  }
  
  // Else decrease C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,1,InBra> vrr_type;
    SafePtr<vrr_type> vrr_ptr(new vrr_type(integral,xyz));
    if (vrr_ptr->num_children())
      return rr_cast(vrr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}

#else

// This approach generates all possible recurrence relations and then
// uses a Tactic object to decide which to use

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr<TwoPRep_11_11_int>& integral,
                                     const SafePtr<Tactic>& tactic)
{
  vector<RR> rrstack;  // stack of all recurrence relations
  
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_ab_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,xyz));
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_cd_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,xyz));
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // decrease A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,xyz));
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  // Else decrease C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,1,InBra> rr_type;
    SafePtr<rr_type> rr_ptr(new rr_type(integral,xyz));
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  return tactic->optimal_rr(rrstack);
}

#endif


