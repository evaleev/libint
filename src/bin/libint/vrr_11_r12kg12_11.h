
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

#ifndef _libint2_src_bin_libint_vrr11r12kg1211_h_
#define _libint2_src_bin_libint_vrr11r12kg1211_h_

using namespace std;

#define USE_OLD_VRRR12kG12_CODE 0

namespace libint2 {

  /** VRR Recurrence Relation for 2-e integrals of the R12_k_G12 operators.
  part specifies the angular momentum of which particle is raised.
  bool bra specifies whether the angular momentum is raised in bra (true) or
  ket (false). I<BFSet,K> is the integral set specialization that describes the
  integrals of the R12_k_G12 operator.
  */
  template <template <typename,typename,typename> class I, class BFSet, int part, FunctionPosition where>
    class VRR_11_R12kG12_11 : public RecurrenceRelation
  {

  public:
    typedef RecurrenceRelation ParentType;
    typedef BFSet BasisFunctionType;
    typedef VRR_11_R12kG12_11<I,BFSet,part,where> ThisType;
    typedef I<BFSet,R12kG12,mType> TargetType;
    typedef TargetType ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).

        dir specifies which quantum number of a and b is shifted.
        For example, dir can be 0 (x), 1(y), or 2(z) if F is
        a Cartesian Gaussian.
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~VRR_11_R12kG12_11() { assert(part == 0 || part == 1); }

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return children_.size(); };
    /// Implementation of RecurrenceRelation::rr_target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target_); }
    /// Implementation of RecurrenceRelation::rr_child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return dynamic_pointer_cast<DGVertex,ChildType>(children_.at(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }
    
    const std::string cpp_function_name() {}
    const std::string cpp_source_name() {}
    const std::string cpp_header_name() {}
    std::ostream& cpp_source(std::ostream&) {}

  private:
    /**
      dir specifies which quantum number is incremented.
      For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
      a Cartesian Gaussian.
     */
    VRR_11_R12kG12_11(const SafePtr<TargetType>&, unsigned int dir);

    // Constructs the RR for K >= 0
    void children_and_expr_Kge0(const vector<BFSet>& bra, const vector<BFSet>& ket,
                                vector<BFSet>* bra_ptr, vector<BFSet>* ket_ptr);
    // Constructs the RR for K == -1
    void children_and_expr_Keqm1(const vector<BFSet>& bra, const vector<BFSet>& ket,
                                 vector<BFSet>* bra_ptr, vector<BFSet>* ket_ptr);

    unsigned int dir_;
    
    SafePtr<TargetType> target_;
    static const unsigned int max_nchildren_ = 8;
    std::vector< SafePtr<ChildType> > children_;
    const SafePtr<ChildType>& make_child(const BFSet& A, const BFSet& B, const BFSet& C, const BFSet& D, unsigned int m, int K = -1) {
      typedef typename TargetType::OperType OperType;
      const OperType K_oper = OperType(K);
      const SafePtr<ChildType>& i = ChildType::Instance(A,B,C,D,m,K_oper);
      children_.push_back(i);
      return *(children_.end()-1);
    }

    /// Implementation of RecurrenceRelation::generate_label()
    std::string generate_label() const
    {
      typedef typename TargetType::AuxIndexType mType;
      static SafePtr<mType> aux0(new mType(0u));
      ostringstream os;
      // R12kG12 VRR recurrence relations codes are independent of m (it never appears anywhere in equations), hence
      // to avoid generating identical code make sure that the (unique) label does not contain m
      os << "VRR Part" << part << " " << to_string(where) << genintegralset_label(target_->bra(),target_->ket(),aux0,target_->oper());
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
  
  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    SafePtr< VRR_11_R12kG12_11<I,F,part,where> >
    VRR_11_R12kG12_11<I,F,part,where>::Instance(const SafePtr< TargetType >& Tint,
                                                  unsigned int dir)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint,dir));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack<ThisType>();
      }
      return this_ptr;
    }
  
  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    VRR_11_R12kG12_11<I,F,part,where>::VRR_11_R12kG12_11(const SafePtr< TargetType >& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir)
    {
      using namespace libint2::algebra;
      using namespace libint2::prefactor;
      children_.reserve(max_nchildren_);
      const int K = target_->oper()->descr().K();
      const unsigned int m = target_->aux()->elem(0);
      const F _1 = unit<F>(dir);

#if !USE_OLD_VRRR12kG12_CODE
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
#else
      F sh_a(Tint->bra(0,0));
      F sh_b(Tint->ket(0,0));
      F sh_c(Tint->bra(1,0));
      F sh_d(Tint->ket(1,0));

      vector<F> bra;
      vector<F> ket;
      bra.push_back(sh_a);
      bra.push_back(sh_c);
      ket.push_back(sh_b);
      ket.push_back(sh_d);

      // Use indirection to choose bra or ket
      vector<F>* bra_ref = &bra;
      vector<F>* ket_ref = &ket;
      if (where == InKet) {
        bra_ref = &ket;
        ket_ref = &bra;
      }
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;

      // See if a-1 exists
      const F& sh_am1 = bra_ref->operator[](p_a) - _1;
      if (!exists(sh_am1)) {
        return;
      }
      
      if (K == -1)
        children_and_expr_Keqm1(bra,ket,bra_ref,ket_ref);
      else
        children_and_expr_Kge0(bra,ket,bra_ref,ket_ref);
#endif
    }

  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    void
    VRR_11_R12kG12_11<I,F,part,where>::children_and_expr_Keqm1(const vector<F>& bra, const vector<F>& ket,
                                                               vector<F>* bra_ref, vector<F>* ket_ref)
    {
#if USE_OLD_VRRR12kG12_CODE
      typedef typename TargetType::OperType OperType;
      const int K = target_->oper()->descr().K();
      const OperType K_oper = OperType(K);
      const unsigned int m = target_->aux()->elem(0);
      const F _1 = unit<F>(dir_);
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;
      // Get a-1
      bra_ref->operator[](p_a).dec(dir_);
      
      children_[0] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],mType(m),K_oper);
      children_[1] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],mType(m+1),K_oper);
      nchildren_ += 2;
      if (is_simple()) {
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.XY_X[part][where][dir_],rr_child(0)));
        SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.W_XY[part][dir_],rr_child(1)));
	SafePtr<ExprType> expr0p1_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr0_ptr,expr1_ptr));
	nflops_ += (3);
        add_expr(expr0p1_ptr);
      }

      // See if a-2 exists
      const F& sh_am2 = bra_ref->operator[](p_a) - _1;
      const bool a_minus_2_exists = exists(sh_am2);
      if (a_minus_2_exists) {
        bra_ref->operator[](p_a).dec(dir_);
        int next_child = nchildren_;
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],m,K_oper);
        children_[next_child+1] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],m+1,K_oper);
        bra_ref->operator[](p_a).inc(dir_);
        const unsigned int ni_a = bra_ref->operator[](p_a).qn(dir_);
        nchildren_ += 2;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.rho_o_alpha12[part], rr_child(next_child+1)));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Minus, rr_child(next_child), expr_intmd0));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_a], prefactors.one_o_2alpha12[part]));
          SafePtr<ExprType> expr2_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd2, expr_intmd1));
	  nflops_ += (5);
          add_expr(expr2_ptr);
        }
      }

      // See if b-1 exists
      const F& sh_bm1 = ket_ref->operator[](p_a) - _1;
      const bool b_minus_1_exists = exists(sh_bm1);
      if (b_minus_1_exists) {
        ket_ref->operator[](p_a).dec(dir_);
        int next_child = nchildren_;
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],m,K_oper);
        children_[next_child+1] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],m+1,K_oper);
        ket_ref->operator[](p_a).inc(dir_);
        const unsigned int ni_b = ket_ref->operator[](p_a).qn(dir_);
        nchildren_ += 2;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.rho_o_alpha12[part], rr_child(next_child+1)));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Minus, rr_child(next_child), expr_intmd0));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_b], prefactors.one_o_2alpha12[part]));
          SafePtr<ExprType> expr3_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd2, expr_intmd1));
	  nflops_ += (5);
          add_expr(expr3_ptr);
        }
      }

      // See if c-1 exists
      const F& sh_cm1 = bra_ref->operator[](p_c) - _1;
      const bool c_minus_1_exists = exists(sh_cm1);
      if (c_minus_1_exists) {
        bra_ref->operator[](p_c).dec(dir_);
        int next_child = nchildren_;
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],m+1,K_oper);
        bra_ref->operator[](p_c).inc(dir_);
        const unsigned int ni_c = bra_ref->operator[](p_c).qn(dir_);
        nchildren_ += 1;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_c], prefactors.one_o_2alphasum));
          SafePtr<ExprType> expr4_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, rr_child(next_child)));
	  nflops_ += (3);
          add_expr(expr4_ptr);
        }
      }

      // See if d-1 exists
      const F& sh_dm1 = ket_ref->operator[](p_c) - _1;
      const bool d_minus_1_exists = exists(sh_dm1);
      if (d_minus_1_exists) {
        ket_ref->operator[](p_c).dec(dir_);
        int next_child = nchildren_;
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],m+1,K_oper);
        ket_ref->operator[](p_c).inc(dir_);
        const unsigned int ni_d = ket_ref->operator[](p_c).qn(dir_);
        nchildren_ += 1;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_d], prefactors.one_o_2alphasum));
          SafePtr<ExprType> expr5_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, rr_child(next_child)));
	  nflops_ += (3);
          add_expr(expr5_ptr);
        }
      }
#endif
    }
    
  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    void
    VRR_11_R12kG12_11<I,F,part,where>::children_and_expr_Kge0(const vector<F>& bra, const vector<F>& ket,
                                                                vector<F>* bra_ref, vector<F>* ket_ref)
    {
#if USE_OLD_VRRR12kG12_CODE
      typedef typename TargetType::OperType OperType;
      const int K = target_->oper()->descr().K();
      const OperType K_oper = OperType(K);
      const OperType Km2_oper = OperType(K-2);
      const unsigned int m = target_->aux()->elem(0);
      const F _1 = unit<F>(dir_);
      if (m != 0)
        throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- nonzero auxiliary quantum detected.");
      
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;
      // Get a-1
      bra_ref->operator[](p_a).dec(dir_);
      
      children_[0] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0u,K_oper);
      nchildren_ += 1;
      if (is_simple()) {
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.R12kG12VRR_pfac0[part][dir_],rr_child(0)));
	nflops_ += (1);
        add_expr(expr0_ptr);
      }

      // See if a-2 exists
      const F& sh_am2 = bra_ref->operator[](p_a) - _1;
      const bool a_minus_2_exists = exists(sh_am2);
      if (a_minus_2_exists) {
        bra_ref->operator[](p_a).dec(dir_);
        children_[1] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0u,K_oper);
        bra_ref->operator[](p_a).inc(dir_);
        const unsigned int ni_a = bra_ref->operator[](p_a).qn(dir_);
        nchildren_ += 1;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.R12kG12VRR_pfac1[part], rr_child(1)));
          SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_a], expr_intmd0));
	  nflops_ += (3);
          add_expr(expr1_ptr);
        }
      }

      // See if b-1 exists
      const F& sh_bm1 = ket_ref->operator[](p_a) - _1;
      const bool b_minus_1_exists = exists(sh_bm1);
      if (b_minus_1_exists) {
        throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- AM on centers b and d must be zero, general RR is not yet implemented");
      }

      // See if c-1 exists
      const F& sh_cm1 = bra_ref->operator[](p_c) - _1;
      const bool c_minus_1_exists = exists(sh_cm1);
      if (c_minus_1_exists) {
        bra_ref->operator[](p_c).dec(dir_);
        int next_child = nchildren_;
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0u,K_oper);
        bra_ref->operator[](p_c).inc(dir_);
        const unsigned int ni_c = bra_ref->operator[](p_c).qn(dir_);
        nchildren_ += 1;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_c], prefactors.R12kG12VRR_pfac2));
          SafePtr<ExprType> expr2_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, rr_child(next_child)));
	  nflops_ += (3);
          add_expr(expr2_ptr);
        }
      }

      // See if d-1 exists
      const F& sh_dm1 = ket_ref->operator[](p_c) - _1;
      const bool d_minus_1_exists = exists(sh_dm1);
      if (d_minus_1_exists) {
        throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- AM on centers b and d must be zero, general RR is not yet implemented");
      }
      
      if (K != 0) {
        int next_child = nchildren_;
        bra_ref->operator[](p_a).inc(dir_);
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0u,Km2_oper);
        bra_ref->operator[](p_a).dec(dir_);
        bra_ref->operator[](p_c).inc(dir_);
        children_[next_child+1] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0u,Km2_oper);
        bra_ref->operator[](p_c).dec(dir_);
        children_[next_child+2] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0u,Km2_oper);
        nchildren_ += 3;
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Minus, rr_child(next_child), rr_child(next_child+1)));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Times, prefactors.R12kG12VRR_pfac4[part][dir_], rr_child(next_child+2)));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Plus, expr_intmd0, expr_intmd1));
          SafePtr<ExprType> expr_intmd3(new ExprType(ExprType::OperatorTypes::Times, prefactors.R12kG12VRR_pfac3[part],expr_intmd2));
          SafePtr<ExprType> expr_intmd4(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[K], expr_intmd3));
	  nflops_ += (6);
          add_expr(expr_intmd4);
        }
      }
#endif
    }

#if LIBINT_ENABLE_GENERIC_CODE
  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    bool
    VRR_11_R12kG12_11<I,F,part,where>::has_generic(const SafePtr<CompilationParameters>& cparams) const
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

  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    std::string
    VRR_11_R12kG12_11<I,F,part,where>::generic_header() const
    {
      const int K = target_->oper()->descr().K();
      if (K == -1)
        return std::string("OSVRR_xs_xs.h");
      else
        return std::string("VRR_r12kg12_xs_xs.h");
    }

  template <template <typename,typename,typename> class I, class F, int part, FunctionPosition where>
    std::string
    VRR_11_R12kG12_11<I,F,part,where>::generic_instance(const SafePtr<CodeContext>& context, const SafePtr<CodeSymbols>& args) const
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
