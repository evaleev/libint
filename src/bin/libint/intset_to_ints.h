
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>
#include <iter.h>

#ifndef _libint2_src_bin_libint_intsettoints_h_
#define _libint2_src_bin_libint_intsettoints_h_

using namespace std;


namespace libint2 {

  /** IntegralSet_to_Integrals converts I, a set of integrals, to individual integrals. Although this is
  technically not a recurrence relation, it can be treated as one.
  */
  template <class I>
  class IntegralSet_to_Integrals : public RecurrenceRelation {
  public:
    typedef I TargetType;
    typedef typename I::iter_type ChildType;

    IntegralSet_to_Integrals(const SafePtr<I>&);
    ~IntegralSet_to_Integrals() {};

    const unsigned int num_children() const { return num_actual_children_; };
    /// target() returns points to the i-th child
    SafePtr<TargetType> target() { return target_; };
    /// child(i) returns points i-th child
    SafePtr<ChildType> child(unsigned int i);

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  private:

    SafePtr<TargetType> target_;
    vector< SafePtr<ChildType> > children_;

    unsigned int num_actual_children_;
  };
  
  
  template <class I>
    IntegralSet_to_Integrals<I>::IntegralSet_to_Integrals(const SafePtr<I>& Tint) :
    target_(Tint)
    {
      target_ = Tint;
      unsigned int m = Tint->m();

      typedef typename I::BraType IBraType;
      typedef typename I::KetType IKetType;
      IBraType* bra = new IBraType(Tint->bra());
      IKetType* ket = new IKetType(Tint->ket());

      // Construct a subiterator for I
      SubIteratorBase<I> siter(Tint);
      num_actual_children_ = siter.num_iter();
      
      // Set children pointers
      for(siter.init(); siter; ++siter)
        children_.push_back(siter.elem());
    };

  template <class I>
    SafePtr<typename I::iter_type>
    IntegralSet_to_Integrals<I>::child(unsigned int i)
    {
      assert(i>=0 && i<num_actual_children_);
      return children_.at(i);
    };
  

};

#endif
