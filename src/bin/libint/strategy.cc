
#define USE_STR1 1

#include <vector>
#include <algorithm>
#include <strategy.h>
#include <dg.h>

using namespace std;
using namespace libint2;

SafePtr<RecurrenceRelation>
Strategy::optimal_rr(const SafePtr<DirectedGraph>& graph, const SafePtr<DGVertex>& integral)
{

  // We must first determine the type of the integral
  {
    SafePtr<TwoPRep_11_11_sq> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_sq,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_sq(graph,eri_ptr);
  }

  // We must first determine the type of the integral
  {
    SafePtr<TwoPRep_11_11_int> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_int,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_int(graph,eri_ptr);
  }

  // Don't know how to apply any RR
  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_sq(const SafePtr<DirectedGraph>& graph,
                                    const SafePtr<TwoPRep_11_11_sq>& integral)
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


#if USE_STR1

// This approach blindly seeks the first possible method to apply VRR

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr<TwoPRep_11_11_int>& integral)
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

// This approach tries all possible methods and chooses the one which generates the fewest number of new vertices

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr<TwoPRep_11_11_int>& integral)
{
  // Apply VRR is every possible way and find the best one (which increases the number of integrals on graph the least)
  vector<unsigned int> nchildren; // number of children of each RR already on the graph
  vector<RR> rrvec;  // recurrence relations
  
  // decrease A
  for(int xyz = 2; xyz>= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,0,InBra> vrr_type;
    SafePtr<vrr_type> vrr_ptr(new vrr_type(integral,xyz));
    if (vrr_ptr->num_children() != 0) {
      rrvec.push_back(rr_cast(vrr_ptr));
      unsigned int nchildren_on_graph = graph->num_children_on(rr_ptr);
      nchildren.push_back(nchildren_on_graph);
    }
  }
  
  // Else decrease C
  for(int xyz = 2; xyz>= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGF,1,InBra> vrr_type;
    SafePtr<vrr_type> vrr_ptr(new vrr_type(integral,xyz));
    if (vrr_ptr->num_children() != 0) {
      rrvec.push_back(rr_cast(vrr_ptr));
      rrvec.push_back(rr_ptr);
      unsigned int nchildren_on_graph = graph->num_children_on(rr_ptr);
      nchildren.push_back(nchildren_on_graph);
    }
  }

  if (rrvec.size()) {
    // Search through nchildren and find the first largest number
    vector<unsigned int>::iterator max_elem = max_element(nchildren.begin(),nchildren.end());
    int pos = distance(nchildren.begin(),max_elem);
    return rrvec.at(pos);
  }
  else
    // Else return null pointer
    return SafePtr<RecurrenceRelation>();
}

#endif
