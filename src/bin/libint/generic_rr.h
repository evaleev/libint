
#ifndef _libint2_src_bin_libint_genericrr_h_
#define _libint2_src_bin_libint_genericrr_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <boost/type_traits/is_same.hpp>

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

  /** RRImpl must inherit GenericRecurrenceRelation<RRImpl>
   */
  template <typename RRImpl, typename F, typename Target>
  class GenericRecurrenceRelation : public RecurrenceRelation {
    public:
      typedef F BasisFunctionType;
      typedef Target TargetType;
      typedef RecurrenceRelation ParentType;
      typedef ParentType::ExprType ExprType;
      
      /// Return an instance if applicable, or a null pointer otherwise
      static SafePtr<RRImpl> Instance(const SafePtr<TargetType>& Tint, unsigned int dir) {
        // screen out calls with nondefault extra parameters
        if (!RRImpl::directional() && dir != 0)
          return SafePtr<RRImpl>();
        // attempt to construct
        SafePtr<RRImpl> this_ptr(new RRImpl(Tint,dir));
        // if succeeded (nchildren > 0) do post-construction
        if (this_ptr->num_children() != 0) {
          this_ptr->register_with_rrstack<RRImpl>();
          return this_ptr;
        }
        // else return null pointer
        return SafePtr<RRImpl>();
      }

      /// Implementation of RecurrenceRelation::num_children()
      const unsigned int num_children() const { return children_.size(); };
      /// Implementation of RecurrenceRelation::rr_target()
      SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target_); }
      /// Implementation of RecurrenceRelation::rr_child()
      SafePtr<DGVertex> rr_child(unsigned int i) const {
        return children_.at(i);
      }
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
      /** is this recurrence relation parameterized by a direction.
          the default is false if BasisFunctionSet is CGShell,
          true otherwise. */
      static bool default_directional() {
        if (boost::is_same<BasisFunctionType,CGShell>::value)
          return false;
        return true;
      }

      /// add child
      const SafePtr<DGVertex>& add_child(const SafePtr<DGVertex>& child) {
        typedef std::vector< SafePtr<DGVertex> > cvector;
        typedef typename cvector::const_iterator citer;
        const citer pos = std::find(children_.begin(),children_.end(),child);
        if (pos == children_.end()) {
          children_.push_back(child);
          return *(children_.rbegin());
        }
        else
          return *pos;
      }

      /// use this helper to make children
      template<class RR, class C> friend class ChildFactory;
      /// make_child should really looks something like this, but gcc 4.3.0 craps out
      /// TODO test is this works
#if 0
      template <class RealChildType>
      const SafePtr<DGVertex>& make_child(const typename RealChildType::BasisFunctionType& A,
                                          const typename RealChildType::BasisFunctionType& B,
                                          const typename RealChildType::BasisFunctionType& C,
                                          const typename RealChildType::BasisFunctionType& D,
                                          const typename RealChildType::AuxIndexType& aux = typename RealChildType::AuxIndexType(),
                                          const typename RealChildType::OperType& oper = typename RealChildType::OperType()) {
        const SafePtr<DGVertex>& i = static_pointer_cast<DGVertex,RealChildType>(ChildType::Instance(A,B,C,D,aux,oper));
        return add_child(i);
      }
#endif
      
      SafePtr<TargetType> target_;
      
    private:
      unsigned int dir_;
      std::vector< SafePtr<DGVertex> > children_;

  };

  /// Helps GenericRecurrenceRelation to work around the compiler problem with make_child
  template <class GenRR, class ChildType>
  class ChildFactory {
    public:
    typedef typename ChildType::BasisFunctionType F;
    typedef typename ChildType::AuxIndexType AuxIndexType;
    typedef typename ChildType::OperType OperType;
      
    ChildFactory(GenRR* rr) : rr_(rr) {}
    
    /// make_child
    const SafePtr<DGVertex>& make_child(const F& A,
                                        const F& B,
                                        const F& C,
                                        const F& D,
                                        const AuxIndexType& aux = AuxIndexType(),
                                        const OperType& oper = OperType()) {
      auto i = static_pointer_cast<DGVertex,ChildType>(ChildType::Instance(A,B,C,D,aux,oper));
      return rr_->add_child(i);
    }
    /// make a child from a wedge of physicists' brackets
    const SafePtr<DGVertex>&
    make_child(const algebra::Wedge< BraketPair<F,PBra>, BraketPair<F,PKet> >& braket_wedge,
               const AuxIndexType& aux = AuxIndexType(),
               const OperType& oper = OperType()) {
      auto i = static_pointer_cast<DGVertex,ChildType>(ChildType::Instance(braket_wedge,aux,oper));
      return rr_->add_child(i);
    }
    /// make a child from a wedge of chemists' brackets
    const SafePtr<DGVertex>&
    make_child(const algebra::Wedge< BraketPair<F,CBra>, BraketPair<F,CKet> >& braket_wedge,
               const AuxIndexType& aux = AuxIndexType(),
               const OperType& oper = OperType()) {
      auto i = static_pointer_cast<DGVertex,ChildType>(ChildType::Instance(braket_wedge,aux,oper));
      return rr_->add_child(i);
    }
    /// take a wedge product of various (linear combinations of) brakets
    void wedge(const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PBra> >& bra_lc,
               const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PKet> >& ket_lc,
               const AuxIndexType& aux = AuxIndexType(),
               const OperType& oper = OperType()) {
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
        auto child = make_child(term.second,aux,oper);
        if (rr_->is_simple()) {
          if (rr_->expr_)
            rr_->expr_ += term.first * child;
          else
            rr_->expr_  = term.first * child;
        }
      }
    }
    void wedge(const BraketPair<F,PBra>& bra,
               const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PKet> >& ket_lc,
               const AuxIndexType& aux = AuxIndexType(),
               const OperType& oper = OperType()) {
      using namespace libint2::prefactor;
      LinearCombination< SafePtr<DGVertex>, BraketPair<F,PBra> > bra_lc;
      bra_lc += make_pair(Scalar(1.0),bra);
      wedge(bra_lc,ket_lc,aux,oper);
    }
    void wedge(const LinearCombination< SafePtr<DGVertex>, BraketPair<F,PBra> >& bra_lc,
               const BraketPair<F,PKet>& ket,
               const AuxIndexType& aux = AuxIndexType(),
               const OperType& oper = OperType()) {
      using namespace libint2::prefactor;
      LinearCombination< SafePtr<DGVertex>, BraketPair<F,PKet> > ket_lc;
      ket_lc += make_pair(Scalar(1.0),ket);
      wedge(bra_lc,ket_lc,aux,oper);
    }
    
    private:
      GenRR* rr_;
  };
  
};

#endif
