

#ifndef _libint2_src_bin_libint_uncontract_h_
#define _libint2_src_bin_libint_uncontract_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>
#include <iter.h>
#include <algebra.h>
#include <context.h>
#include <dims.h>

using namespace std;


namespace libint2 {

  /** Uncontract_Integral_base is dummy class used for dynamic casts only
   */
  class Uncontract_Integral_base {
  protected:
    virtual ~Uncontract_Integral_base() {}
  };

  /** Uncontract_Integral converts (a set of) contracted integral(s) to its uncontracted counterpart. Although this is
  technically not a recurrence relation, it can be treated as one.
  */
  template <class I>
  class Uncontract_Integral : public RecurrenceRelation,
                              public Uncontract_Integral_base {
  public:
    typedef I TargetType;
    typedef TargetType ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    Uncontract_Integral(const SafePtr<I>&);
    virtual ~Uncontract_Integral() {}

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return children_.size(); };
    /// target() returns pointer to target
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation's target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation's child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    // can't inline it -- will always generate a call
    bool is_simple() const {
      return false;
    }

  private:
    SafePtr<TargetType> target_;
    vector< SafePtr<ChildType> > children_;

    // contracting integrals only depends on the number of integrals in a set
    // can simplify this to only refer to that number
    std::string generate_label() const {
      ostringstream os;
      os << "Generic Contract";
      return os.str();
    }

    //
    std::string spfunction_call(const SafePtr<CodeContext>& context,
                                const SafePtr<ImplicitDimensions>& dims) const;

  };


  template <class I>
  Uncontract_Integral<I>::Uncontract_Integral(const SafePtr<I>& Tint) :
    target_(Tint)
    {
      target_ = Tint;

      // create the uncontracted version
      typedef typename I::BraType bratype;
      typedef typename I::KetType kettype;
      typedef typename I::OperType opertype;
      bratype bra_unc = target_->bra();
      kettype ket_unc = target_->ket();
      // is it even contracted?
      bool target_is_contracted = false;
      {
        const unsigned int np = bra_unc.num_part();
        for(unsigned int p=0; p<np; ++p) {
          const unsigned int nf = bra_unc.num_members(p);
          for(unsigned int f=0; f<nf; ++f) {
            target_is_contracted |= bra_unc.member(p,f).contracted();
            bra_unc.member(p,f).uncontract();
          }
        }
      }
      {
        const unsigned int np = ket_unc.num_part();
        for(unsigned int p=0; p<np; ++p) {
          const unsigned int nf = ket_unc.num_members(p);
          for(unsigned int f=0; f<nf; ++f) {
            target_is_contracted |= ket_unc.member(p,f).contracted();
            ket_unc.member(p,f).uncontract();
          }
        }
      }
      opertype oper_unc = target_->oper();
      const bool oper_contracted = oper_unc.descr().contracted();
      target_is_contracted |= oper_contracted;
      oper_unc.descr().uncontract();

      if (target_is_contracted) {
#if DEBUG
        std::cout << "Uncontract_Integral: " << target_->description() << " is contracted" << std::endl;
#endif
        SafePtr<ChildType> c = ChildType::Instance(bra_unc, ket_unc,
                                                   target_->aux(),
                                                   oper_unc);
        children_.push_back(c);
      }
    };

  template <class I>
    SafePtr<typename Uncontract_Integral<I>::ChildType>
    Uncontract_Integral<I>::child(unsigned int i) const
    {
      return children_.at(i);
    };

  template <class I>
  std::string
  Uncontract_Integral<I>::spfunction_call(const SafePtr<CodeContext>& context,
                                          const SafePtr<ImplicitDimensions>& dims) const {
    ostringstream os;
    // contraction = reduction
    os << "_libint2_static_api_inc1_short_("
        << context->value_to_pointer(rr_target()->symbol()) << ","
        << context->value_to_pointer(rr_child(0)->symbol()) << ","
       << target_->size()
       << ")" << context->end_of_stat() << endl;
    unsigned int& nflops_ref = const_cast<unsigned int&>(nflops_);
    nflops_ref += target_->size();

    return os.str();
  }

};

#endif

