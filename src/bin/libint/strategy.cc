
#include <strategy.h>

using namespace std;
using namespace libint2;

SafePtr<RecurrenceRelation>
Strategy::optimal_rr(const SafePtr<DGVertex>& integral)
{

  // We must first determine the type of the integral
  {
    SafePtr<TwoPRep_11_11_sq > eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_sq,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_sq(eri_ptr);
  }

  // We must first determine the type of the integral
  {
    SafePtr<TwoPRep_11_11_int > eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_int,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_int(eri_ptr);
  }

  // Don't know how to apply any RR
  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_sq(const SafePtr<TwoPRep_11_11_sq>& integral)
{
  //
  // This is a basic strategy for computing integral
  // 1) first see if should convert the set to infividual integrals
  // 2) if possible apply HRR
  // 3) else apply VRR
  //
  if (integral->size() <= max_size_to_unroll)
    return unroll_intset<TwoPRep_11_11_sq>(integral);

  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<TwoPRep_11_11_int>& integral)
{
  return SafePtr<RecurrenceRelation>();
}

