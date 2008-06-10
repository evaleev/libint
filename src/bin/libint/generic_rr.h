
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
#include <r12kg12_11_11.h>
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
  
  template <class BFSet, int part, FunctionPosition where>
    class VRR_new_11_R12kG12_11 : public GenericRecurrenceRelation< VRR_new_11_R12kG12_11<BFSet,part,where>,
                                                                    BFSet,
                                                                    GenIntegralSet_11_11<BFSet,R12kG12,mType>,
                                                                    GenIntegralSet_11_11<BFSet,R12kG12,mType> >
  {
  public:
    typedef BFSet BasisFunctionType;
    typedef GenIntegralSet_11_11<BFSet,R12kG12,mType> TargetType;
    typedef TargetType ChildType;
    typedef GenericRecurrenceRelation<VRR_new_11_R12kG12_11,BFSet,TargetType,ChildType> ParentType;
    friend class GenericRecurrenceRelation<VRR_new_11_R12kG12_11,BFSet,TargetType,ChildType>;
    static const unsigned int max_nchildren = 8;

    using ParentType::Instance;
    static std::string descr() { return "VRR"; }

  private:
    using RecurrenceRelation::expr_;
    using RecurrenceRelation::nflops_;
    using ParentType::target_;
    using ParentType::is_simple;

    /**
      dir specifies which quantum number is incremented.
      For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
      a Cartesian Gaussian.
     */
    VRR_new_11_R12kG12_11(const SafePtr<TargetType>&, unsigned int dir);
    
    const SafePtr<ChildType>& make_child(const BFSet& A, const BFSet& B, const BFSet& C, const BFSet& D, unsigned int m, int K = -1) {
      typedef typename TargetType::OperType OperType;
      const OperType K_oper = OperType(K);
      const SafePtr<ChildType>& i = ChildType::Instance(A,B,C,D,m,K_oper);
      return add_child(i);
    }

    /// Implementation of RecurrenceRelation::generate_label()
    std::string generate_label() const
    {
      typedef typename TargetType::AuxIndexType mType;
      static SafePtr<mType> aux0(new mType(0u));
      ostringstream os;
      // R12kG12 VRR recurrence relations codes are independent of m (it never appears anywhere in equations), hence
      // to avoid generating identical code make sure that the (unique) label does not contain m
      os << descr() << part << to_string(where)
         << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
      return os.str();
    }

#if LIBINT_ENABLE_GENERIC_CODE
    /// Implementation of RecurrenceRelation::has_generic()
    bool has_generic(const SafePtr<CompilationParameters>& cparams) const;
    /// Implementation of RecurrenceRelation::generic_header()
    std::string generic_header() const;
    /// Implementation of RecurrenceRelation::generic_instance()
    std::string generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const;
#endif
  };
  
  template <class F, int part, FunctionPosition where>
    VRR_new_11_R12kG12_11<F,part,where>::VRR_new_11_R12kG12_11(const SafePtr< TargetType >& Tint,
                                                               unsigned int dir) :
    ParentType(Tint,dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      const int K = Tint->oper()->descr().K();
      const unsigned int m = Tint->aux()->elem(0);
      const F _1 = unit<F>(dir);

      // if K is -1, the recurrence relation looks exactly as it would for ERI
      // thus generate the same code, and remember to use appropriate prefactors
      if (K == -1) {
        // Build on A
        if (part == 0 && where == InBra) {
          F a(Tint->bra(0,0) - _1);
          if (!exists(a)) return;
          F b(Tint->ket(0,0));
          F c(Tint->bra(1,0));
          F d(Tint->ket(1,0));

          const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
          const SafePtr<ChildType>& ABCD_mp1 = make_child(a,b,c,d,m+1);
          if (is_simple()) { expr_ = Vector("PA")[dir] * ABCD_m + Vector("WP")[dir] * ABCD_mp1;  nflops_+=3; }

          const F& am1 = a - _1;
          if (exists(am1)) {
            const SafePtr<ChildType>& Am1BCD_m = make_child(am1,b,c,d,m);
            const SafePtr<ChildType>& Am1BCD_mp1 = make_child(am1,b,c,d,m+1);
            if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2z") * (Am1BCD_m - Scalar("roz") * Am1BCD_mp1);  nflops_+=5; }
          }
          const F& bm1 = b - _1;
          if (exists(bm1)) {
            const SafePtr<ChildType>& ABm1CD_m = make_child(a,bm1,c,d,m);
            const SafePtr<ChildType>& ABm1CD_mp1 = make_child(a,bm1,c,d,m+1);
            if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2z") * (ABm1CD_m - Scalar("roz") * ABm1CD_mp1);  nflops_+=5; }
          }
          const F& cm1 = c - _1;
          if (exists(cm1)) {
            const SafePtr<ChildType>& ABCm1D_mp1 = make_child(a,b,cm1,d,m+1);
            if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2ze") * ABCm1D_mp1;  nflops_+=3; }
          }
          const F& dm1 = d - _1;
          if (exists(dm1)) {
            const SafePtr<ChildType>& ABCDm1_mp1 = make_child(a,b,c,dm1,m+1);
            if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2ze") * ABCDm1_mp1;  nflops_+=3; }
          }
          return;
        }
        // Build on C
        if (part == 1 && where == InBra) {
          F a(Tint->bra(0,0));
          F b(Tint->ket(0,0));
          F c(Tint->bra(1,0) - _1);
          if (!exists(c)) return;
          F d(Tint->ket(1,0));

          const SafePtr<ChildType>& ABCD_m = make_child(a,b,c,d,m);
          const SafePtr<ChildType>& ABCD_mp1 = make_child(a,b,c,d,m+1);
          if (is_simple()) { expr_ = Vector("QC")[dir] * ABCD_m + Vector("WQ")[dir] * ABCD_mp1;  nflops_+=3; }

          const F& cm1 = c - _1;
          if (exists(cm1)) {
            const SafePtr<ChildType>& ABCm1D_m = make_child(a,b,cm1,d,m);
            const SafePtr<ChildType>& ABCm1D_mp1 = make_child(a,b,cm1,d,m+1);
            if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("oo2e") * (ABCm1D_m - Scalar("roe") * ABCm1D_mp1);  nflops_+=5; }
          }
          const F& dm1 = d - _1;
          if (exists(dm1)) {
            const SafePtr<ChildType>& ABCDm1_m = make_child(a,b,c,dm1,m);
            const SafePtr<ChildType>& ABCDm1_mp1 = make_child(a,b,c,dm1,m+1);
            if (is_simple()) { expr_ += Vector(d)[dir] * Scalar("oo2e") * (ABCDm1_m - Scalar("roe") * ABCDm1_mp1);  nflops_+=5; }
          }
          const F& am1 = a - _1;
          if (exists(am1)) {
            const SafePtr<ChildType>& Am1BCD_mp1 = make_child(am1,b,c,d,m+1);
            if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("oo2ze") * Am1BCD_mp1;  nflops_+=3; }
          }
          const F& bm1 = b - _1;
          if (exists(bm1)) {
            const SafePtr<ChildType>& ABm1CD_mp1 = make_child(a,bm1,c,d,m+1);
            if (is_simple()) { expr_ += Vector(b)[dir] * Scalar("oo2ze") * ABm1CD_mp1;  nflops_+=3; }
          }
          return;
        }
        return;
      } // K == -1?
      else {
        // K != -1, the auxiliary quantum number is not used
        if (m != 0)
          throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- nonzero auxiliary quantum detected.");
        // B and D must be s functions -- general RR not implemented yet
        {
          F b(Tint->ket(0,0) - _1);
          F d(Tint->ket(1,0) - _1);
          if (exists(b) || exists(d))
            return;
        }
        
        // Build on A
        if (part == 0 && where == InBra) {
          F a(Tint->bra(0,0) - _1);
          if (!exists(a)) return;
          F b(Tint->ket(0,0));
          F c(Tint->bra(1,0));
          F d(Tint->ket(1,0));

          const SafePtr<ChildType>& ABCD_K = make_child(a,b,c,d,0u,K);
          if (is_simple()) { expr_ = Vector("R12kG12_pfac0_0")[dir] * ABCD_K;  nflops_+=1; }
          const F& am1 = a - _1;
          if (exists(am1)) {
            const SafePtr<ChildType>& Am1BCD_K = make_child(am1,b,c,d,0u,K);
            if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("R12kG12_pfac1_0") * Am1BCD_K;  nflops_+=3; }
          }
          const F& cm1 = c - _1;
          if (exists(cm1)) {
            const SafePtr<ChildType>& ABCm1D_K = make_child(a,b,cm1,d,0u,K);
            if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("R12kG12_pfac2") * ABCm1D_K;  nflops_+=3; }
          }
          if (K != 0) {
            const SafePtr<ChildType>& Ap1BCD_Km2 = make_child(a+_1,b,c,d,0u,K-2);
            const SafePtr<ChildType>& ABCp1D_Km2 = make_child(a,b,c+_1,d,0u,K-2);
            const SafePtr<ChildType>& ABCD_Km2 = make_child(a,b,c,d,0u,K-2);
            if (is_simple()) { expr_ += Scalar((double)K) * Scalar("R12kG12_pfac3_0")
                                        * (Ap1BCD_Km2 - ABCp1D_Km2 + Vector("R12kG12_pfac4_0")[dir] * ABCD_Km2);  nflops_+=6; }
          }
          return;
        }
        
        // Build on C
        if (part == 1 && where == InBra) {
          F a(Tint->bra(0,0));
          F b(Tint->ket(0,0));
          F c(Tint->bra(1,0) - _1);
          if (!exists(c)) return;
          F d(Tint->ket(1,0));

          const SafePtr<ChildType>& ABCD_K = make_child(a,b,c,d,0u,K);
          if (is_simple()) { expr_ = Vector("R12kG12_pfac0_1")[dir] * ABCD_K;  nflops_+=1; }
          const F& cm1 = c - _1;
          if (exists(cm1)) {
            const SafePtr<ChildType>& ABCm1D_K = make_child(a,b,cm1,d,0u,K);
            if (is_simple()) { expr_ += Vector(c)[dir] * Scalar("R12kG12_pfac1_1") * ABCm1D_K;  nflops_+=3; }
          }
          const F& am1 = a - _1;
          if (exists(am1)) {
            const SafePtr<ChildType>& Am1BCD_K = make_child(am1,b,c,d,0u,K);
            if (is_simple()) { expr_ += Vector(a)[dir] * Scalar("R12kG12_pfac2") * Am1BCD_K;  nflops_+=3; }
          }
          if (K != 0) {
            const SafePtr<ChildType>& ABCp1D_Km2 = make_child(a,b,c+_1,d,0u,K-2);
            const SafePtr<ChildType>& Ap1BCD_Km2 = make_child(a+_1,b,c,d,0u,K-2);
            const SafePtr<ChildType>& ABCD_Km2 = make_child(a,b,c,d,0u,K-2);
            if (is_simple()) { expr_ += Scalar((double)K) * Scalar("R12kG12_pfac3_1")
                                        * (ABCp1D_Km2 - Ap1BCD_Km2 + Vector("R12kG12_pfac4_1")[dir] * ABCD_Km2);  nflops_+=6; }
          }
          return;
        }
        return;
      } // K >= 0
    }

#if LIBINT_ENABLE_GENERIC_CODE
  template <class F, int part, FunctionPosition where>
    bool
    VRR_new_11_R12kG12_11<F,part,where>::has_generic(const SafePtr<CompilationParameters>& cparams) const
    {
      F sh_a(target_->bra(0,0));
      F sh_b(target_->ket(0,0));
      F sh_c(target_->bra(1,0));
      F sh_d(target_->ket(1,0));
      const unsigned int max_opt_am = cparams->max_am_opt();
      // generic code works for a0c0 classes where am(a) > 1 and am(c) > 1
      // to generate optimized code for xxxx integral need to generate specialized code for up to (x+x)0(x+x)0 integrals
      if (!TrivialBFSet<F>::result &&
          sh_b.zero() && sh_d.zero() &&
          sh_a.norm() > std::max(2*max_opt_am,1u) && sh_c.norm() > std::max(2*max_opt_am,1u))
        return true;
      else
        return false;
    }

  template <class F, int part, FunctionPosition where>
    std::string
    VRR_new_11_R12kG12_11<F,part,where>::generic_header() const
    {
      const int K = target_->oper()->descr().K();
      if (K == -1)
        return std::string("OSVRR_xs_xs.h");
      else
        return std::string("VRR_r12kg12_xs_xs.h");
    }

  template <class F, int part, FunctionPosition where>
    std::string
    VRR_new_11_R12kG12_11<F,part,where>::generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const
    {
      const int K = target_->oper()->descr().K();
      std::ostringstream oss;
      F sh_a(target_->bra(0,0));
      F sh_c(target_->bra(1,0));

      oss << "using namespace libint2;" << endl;
      if (K == -1) {
        oss << "libint2::OSVRR_xs_xs<" << part << "," << to_string(where) << "," << sh_a.norm() << "," << sh_c.norm() << ",";
        oss << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
      }
      else {
        oss << "libint2::VRR_r12kg12_xs_xs<" << part << "," << to_string(where) << "," << sh_a.norm() << "," << sh_c.norm() << ",";
        oss << K << "," << ((context->cparams()->max_vector_length() == 1) ? "false" : "true");
      }
      oss << ">::compute(inteval";
      
      const unsigned int nargs = args->n();
      for(unsigned int a=0; a<nargs; a++) {
        oss << "," << args->symbol(a);
      }
      
      // if K == 0 add dummy arguments so that the same generic function can be used for all K>=0 cases
      if (K == 0) {
        for(unsigned int a=0; a<3; ++a) {
          oss << ",(LIBINT2_REALTYPE*)0";
        }
      }
      
      oss << ");";
      
      return oss.str();
    }
#endif // #if !LIBINT_ENABLE_GENERIC_CODE
    
};

#endif
