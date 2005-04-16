
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

  /** IntegralSet_to_Integrals_base is dummy class used for dynamic casts only
   */
  class IntegralSet_to_Integrals_base {
  protected:
    virtual ~IntegralSet_to_Integrals_base() {}
  };

  /** IntegralSet_to_Integrals converts I, a set of integrals, to individual integrals. Although this is
  technically not a recurrence relation, it can be treated as one.
  */
  template <class I>
  class IntegralSet_to_Integrals : public RecurrenceRelation,
    public IntegralSet_to_Integrals_base {
  public:
    typedef I TargetType;
    typedef typename I::iter_type ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef AlgebraicOperator<DGVertex> ExprType;

    IntegralSet_to_Integrals(const SafePtr<I>&);
    ~IntegralSet_to_Integrals() {}

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return children_.size(); };
    /// Implementation of RecurrenceRelation::num_expr()
    const unsigned int num_expr() const { return 0; };
    /// target() returns pointer to target
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation's target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation's child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    /// Implementation of RecurrenceRelation::rr_expr()
    SafePtr<DGVertex> rr_expr(unsigned int i) const { return SafePtr<DGVertex>(); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return true;
    }
    /// Implementation of RecurrenceRelation::invariant_type()
    bool invariant_type() const {
      // Converts from one BFSet to another!
      return false;
    }
    /// Implementation of RecurrenceRelation::label()
    std::string label() const {
      throw std::runtime_error("IntegralSet_to_Integrals::label() -- code for this RR is never generated, so this function should never be used");
    }
    /// Implementation of RecurrenceRelation::nflops()
    unsigned int nflops() const { return 0; }

    const std::string cpp_function_name() {};
    const std::string cpp_source_name() {};
    const std::string cpp_header_name() {};
    std::ostream& cpp_source(std::ostream&) {};

  private:
    SafePtr<TargetType> target_;
    vector< SafePtr<ChildType> > children_;
  };
  
  
  template <class I>
    IntegralSet_to_Integrals<I>::IntegralSet_to_Integrals(const SafePtr<I>& Tint) :
    target_(Tint)
    {
      target_ = Tint;

      // Construct a subiterator for I
      SubIteratorBase<I> siter(Tint);
      
      // Set children pointers
      for(siter.init(); siter; ++siter)
        children_.push_back(siter.elem());
    };

  template <class I>
    SafePtr<typename I::iter_type>
    IntegralSet_to_Integrals<I>::child(unsigned int i) const
    {
      return children_.at(i);
    };
  

};

#endif
