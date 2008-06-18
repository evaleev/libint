
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

      /// make a child
      const SafePtr<ChildType>& make_child(const F& A, const F& B, const F& C, const F& D,
                                           const typename ChildType::AuxIndexType& aux = typename ChildType::AuxIndexType(),
                                           const typename ChildType::OperType& oper = typename ChildType::OperType()) {
        const SafePtr<ChildType>& i = ChildType::Instance(A,B,C,D,aux,oper);
        return add_child(i);
      }
      /// make a child from a wedge of physicists' brackets
      const SafePtr<ChildType>&
      make_child(const algebra::Wedge< BraketPair<F,PBra>, BraketPair<F,PKet> >& braket_wedge,
                 const typename ChildType::AuxIndexType& aux = typename ChildType::AuxIndexType(),
                 const typename ChildType::OperType& oper = typename ChildType::OperType()) {
        const SafePtr<ChildType>& i = ChildType::Instance(braket_wedge,aux,oper);
        return add_child(i);
      }
      /// make a child from a wedge of chemists' brackets
      const SafePtr<ChildType>&
      make_child(const algebra::Wedge< BraketPair<F,CBra>, BraketPair<F,CKet> >& braket_wedge,
                 const typename ChildType::AuxIndexType& aux = typename ChildType::AuxIndexType(),
                 const typename ChildType::OperType& oper = typename ChildType::OperType()) {
        const SafePtr<ChildType>& i = ChildType::Instance(braket_wedge,aux,oper);
        return add_child(i);
      }
      /// add child
      const SafePtr<ChildType>& add_child(const SafePtr<ChildType>& child) {
        typedef std::vector< SafePtr<ChildType> > cvector;
        typedef typename cvector::const_iterator citer;
        const citer pos = std::find(children_.begin(),children_.end(),child);
        if (pos == children_.end()) {
          children_.push_back(child);
          return *(children_.rbegin());
        }
        else
          return *pos;
      }
      
      /// take a wedge product of various (linear combinations of) brakets
      void wedge(const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PBra> >& bra_lc,
                 const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PKet> >& ket_lc,
                 const typename ChildType::AuxIndexType& aux = typename ChildType::AuxIndexType(),
                 const typename ChildType::OperType& oper = typename ChildType::OperType()) {
        using namespace libint2::algebra;
        typedef LinearCombination< SafePtr<DGVertex>,
                                   Wedge< BraketPair<F,PBra>,
                                          BraketPair<F,PKet>
                                        >
                                 > ProductLC;
        const ProductLC& product_lc = bra_lc ^ ket_lc;
        const size_t nprod = product_lc.size();
        for(unsigned int t=0; t<nprod; ++t) {
          const typename ProductLC::term_t& term = product_lc[t];
          const SafePtr<ChildType>& child = make_child(term.second,aux,oper);
          if (is_simple()) {
            if (expr_)
              expr_ += term.first * child;
            else
              expr_  = term.first * child;
          }
        }
      }
      void wedge(const BraketPair<F,PBra>& bra,
                 const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PKet> >& ket_lc,
                 const typename ChildType::AuxIndexType& aux = typename ChildType::AuxIndexType(),
                 const typename ChildType::OperType& oper = typename ChildType::OperType()) {
        using namespace libint2::prefactor;
        LinearCombination< SafePtr<DGVertex>, BraketPair<F,PBra> > bra_lc;
        bra_lc += make_pair(Scalar(1.0),bra);
        wedge(bra_lc,ket_lc,aux,oper);
      }
      void wedge(const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PBra> >& bra_lc,
                 const BraketPair<F,PKet>& ket,
                 const typename ChildType::AuxIndexType& aux = typename ChildType::AuxIndexType(),
                 const typename ChildType::OperType& oper = typename ChildType::OperType()) {
        using namespace libint2::prefactor;
        LinearCombination< SafePtr<DGVertex>, BraketPair<F,PKet> > ket_lc;
        ket_lc += make_pair(Scalar(1.0),ket);
        wedge(bra_lc,ket_lc,aux,oper);
      }
      
      SafePtr<TargetType> target_;
      
    private:
      unsigned int dir_;
      std::vector< SafePtr<ChildType> > children_;

  };
    
};

#endif
