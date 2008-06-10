
#ifndef _libint2_src_bin_libint_genericrr_h_
#define _libint2_src_bin_libint_genericrr_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <dgvertex.h>
#include <rr.h>
#include <integral.h>
#include <algebra.h>
#include <flop.h>
#include <prefactors.h>
#include <context.h>
#include <default_params.h>
#include <util.h>

using namespace std;

namespace libint2 {

  /// RRImpl must inherit GenericRecurrenceRelation<RRImpl>
  template <typename RRImpl, typename F, typename Target, typename Child>
  class GenericRecurrenceRelation : public RecurrenceRelation {
    public:
      typedef F BasisFunctionType;
      typedef Target TargetType;
      typedef Child ChildType;
      typedef RecurrenceRelation ParentType;
      typedef ParentType::ExprType ExprType;
      
      static SafePtr<RRImpl> Instance(const SafePtr<TargetType>& Tint, unsigned int dir) {
        SafePtr<RRImpl> this_ptr(new RRImpl(Tint,dir));
        // Do post-construction duties
        if (this_ptr->num_children() != 0) {
          this_ptr->register_with_rrstack<RRImpl>();
        }
        return this_ptr;
      }

      /// Implementation of RecurrenceRelation::num_children()
      const unsigned int num_children() const { return children_.size(); };
      /// Implementation of RecurrenceRelation::rr_target()
      SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target_); }
      /// Implementation of RecurrenceRelation::rr_child()
      SafePtr<DGVertex> rr_child(unsigned int i) const { return dynamic_pointer_cast<DGVertex,ChildType>(children_.at(i)); }
      /// Implementation of RecurrenceRelation::is_simple()
      bool is_simple() const {
        return TrivialBFSet<BasisFunctionType>::result;
      }
      
      /// Implementation of RecurrenceRelation::generate_label()
      std::string generate_label() const
      {
        ostringstream os;
        os << RRImpl::descr() << " " << target_->label();
        return os.str();
      }

    protected:
      GenericRecurrenceRelation(const SafePtr<TargetType>& Tint, unsigned int dir) :
        target_(Tint), dir_(dir) {
        children_.reserve(RRImpl::max_nchildren);
      }

      /// add child
      const SafePtr<ChildType>& add_child(const SafePtr<ChildType>& child) {
        typedef std::vector< SafePtr<ChildType> > cvector;
        typedef typename cvector::const_iterator citer;
        const citer pos = std::find(children_.begin(),children_.end(),child);
        if (pos == children_.end())
          children_.push_back(child);
        return *(children_.rbegin());
      }

      SafePtr<TargetType> target_;
      
    private:
      unsigned int dir_;
      std::vector< SafePtr<ChildType> > children_;

  };
    
};

#endif
