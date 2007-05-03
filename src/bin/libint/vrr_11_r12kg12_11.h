
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

#ifndef _libint2_src_bin_libint_vrr11r12kg1211_h_
#define _libint2_src_bin_libint_vrr11r12kg1211_h_

using namespace std;


namespace libint2 {

  /** VRR Recurrence Relation for 2-e integrals of the R12_k_G12 operators.
  part specifies the angular momentum of which particle is raised.
  bool bra specifies whether the angular momentum is raised in bra (true) or
  ket (false). I<BFSet,K> is the integral set specialization that describes the
  integrals of the R12_k_G12 operator.
  */
  template <template <class,int> class I, class BFSet, int K, int part, FunctionPosition where>
    class VRR_11_R12kG12_11 : public RecurrenceRelation
  {

  public:
    typedef RecurrenceRelation ParentType;
    typedef BFSet BasisFunctionType;
    typedef VRR_11_R12kG12_11<I,BFSet,K,part,where> ThisType;
    typedef I<BFSet,K> TargetType;
    typedef R12kG12_11_11_base<BFSet> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).

        dir specifies which quantum number of a and b is shifted.
        For example, dir can be 0 (x), 1(y), or 2(z) if F is
        a Cartesian Gaussian.
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~VRR_11_R12kG12_11();

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return nchildren_; };
    /// target() returns pointer to the i-th child
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer to the i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation::rr_target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation::rr_child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return dynamic_pointer_cast<DGVertex,ChildType>(child(i)); }
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
    SafePtr<ChildType> children_[max_nchildren_];
    unsigned int nchildren_;

    std::string generate_label() const;
  };
  
  /** class VRR_11_R12kG12_11_util implements the logic of the VRR. For any
      K >= 0 the same RR is used. K=-1 is a special case and must be handled separately (via
      specialization of this template, see vrr_11_r12kg12_11.cc).
  */
  template <int K>
    struct
    VRR_11_R12kG12_11_util {
      static const bool K_eq_m1 = (K==-1);
    };

  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    SafePtr< VRR_11_R12kG12_11<I,F,K,part,where> >
    VRR_11_R12kG12_11<I,F,K,part,where>::Instance(const SafePtr<I<F,K> >& Tint,
                                                  unsigned int dir)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint,dir));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack<ThisType>();
      }
      return this_ptr;
    }
  
  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    VRR_11_R12kG12_11<I,F,K,part,where>::VRR_11_R12kG12_11(const SafePtr<I<F,K> >& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir), nchildren_(0)
    {
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
      if (!bra_ref->operator[](p_a).dec(dir)) {
        return;
      }
      // turn a-1 back to a
      bra_ref->operator[](p_a).inc(dir);
      
      if (VRR_11_R12kG12_11_util<K>::K_eq_m1)
        children_and_expr_Keqm1(bra,ket,bra_ref,ket_ref);
      else
        children_and_expr_Kge0(bra,ket,bra_ref,ket_ref);
    };
  
  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    void
    VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Keqm1(const vector<F>& bra, const vector<F>& ket,
                                                                 vector<F>* bra_ref, vector<F>* ket_ref)
    {
      unsigned int m = target_->m();
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;
      // Get a-1
      bra_ref->operator[](p_a).dec(dir_);
      
      children_[0] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m);
      children_[1] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      nchildren_ += 2;
      if (is_simple()) {
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.XY_X[part][where][dir_],rr_child(0)));
        SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.W_XY[part][dir_],rr_child(1)));
	nflops_ += (3);
        add_expr(expr1_ptr);
      }

      // See if a-2 exists
      bool a_minus_2_exists = true;
      if (!bra_ref->operator[](p_a).dec(dir_)) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        int next_child = nchildren_;
        children_[next_child] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m);
        children_[next_child+1] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
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
      bool b_minus_1_exists = true;
      if (!ket_ref->operator[](p_a).dec(dir_)) {
        b_minus_1_exists = false;
      }
      if (b_minus_1_exists) {
        int next_child = nchildren_;
        children_[next_child] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m);
        children_[next_child+1] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
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
      bool c_minus_1_exists = true;
      if (!bra_ref->operator[](p_c).dec(dir_)) {
        c_minus_1_exists = false;
      }
      if (c_minus_1_exists) {
        int next_child = nchildren_;
        children_[next_child] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
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
      bool d_minus_1_exists = true;
      if (!ket_ref->operator[](p_c).dec(dir_)) {
        d_minus_1_exists = false;
      }
      if (d_minus_1_exists) {
        int next_child = nchildren_;
        children_[next_child] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
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
    }
  
  template <int K>
  struct EqualZero {
    static const bool result = (K==0);
  };
  
  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    void
    VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0(const vector<F>& bra, const vector<F>& ket,
                                                                vector<F>* bra_ref, vector<F>* ket_ref)
    {
      unsigned int m = target_->m();
      if (m != 0)
        throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- nonzero auxiliary quantum detected.");
      
      // On which particle to act
      int p_a = part;
      int p_c = (p_a == 0) ? 1 : 0;
      // Get a-1
      bra_ref->operator[](p_a).dec(dir_);
      
      children_[0] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],0);
      nchildren_ += 1;
      if (is_simple()) {
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.R12kG12VRR_pfac0[part][dir_],rr_child(0)));
	nflops_ += (1);
        add_expr(expr0_ptr);
      }

      // See if a-2 exists
      bool a_minus_2_exists = true;
      if (!bra_ref->operator[](p_a).dec(dir_)) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        children_[1] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],0);
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
      bool b_minus_1_exists = true;
      if (!ket_ref->operator[](p_a).dec(dir_)) {
        b_minus_1_exists = false;
      }
      if (b_minus_1_exists) {
        throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- AM on centers b and d must be zero, general RR is not yet implemented");
      }

      // See if c-1 exists
      bool c_minus_1_exists = true;
      if (!bra_ref->operator[](p_c).dec(dir_)) {
        c_minus_1_exists = false;
      }
      if (c_minus_1_exists) {
        int next_child = nchildren_;
        children_[next_child] = I<F,K>::Instance(bra[0],ket[0],bra[1],ket[1],0);
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
      bool d_minus_1_exists = true;
      if (!ket_ref->operator[](p_c).dec(dir_)) {
	d_minus_1_exists = false;
      }
      if (d_minus_1_exists) {
	throw std::logic_error("VRR_11_R12kG12_11<I,F,K,part,where>::children_and_expr_Kge0() -- AM on centers b and d must be zero, general RR is not yet implemented");
      }
      
      if (EqualZero<K>::result == false) {
        int next_child = nchildren_;
        bra_ref->operator[](p_a).inc(dir_);
        children_[next_child] = I<F,K-2>::Instance(bra[0],ket[0],bra[1],ket[1],0);
        bra_ref->operator[](p_a).dec(dir_);
        bra_ref->operator[](p_c).inc(dir_);
        children_[next_child+1] = I<F,K-2>::Instance(bra[0],ket[0],bra[1],ket[1],0);
        bra_ref->operator[](p_c).dec(dir_);
        children_[next_child+2] = I<F,K-2>::Instance(bra[0],ket[0],bra[1],ket[1],0);
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
    }

  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    VRR_11_R12kG12_11<I,F,K,part,where>::~VRR_11_R12kG12_11()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    SafePtr< typename VRR_11_R12kG12_11<I,F,K,part,where>::ChildType >
    VRR_11_R12kG12_11<I,F,K,part,where>::child(unsigned int i) const
    {
      assert(i>=0 && i<nchildren_);
      unsigned int nc=0;
      for(int c=0; c<max_nchildren_; c++) {
        if (children_[c] != 0) {
          if (nc == i)
            return children_[c];
          nc++;
        }
      }
    };

  template <template <class,int> class I, class F, int K, int part, FunctionPosition where>
    std::string
    VRR_11_R12kG12_11<I,F,K,part,where>::generate_label() const
    {
      ostringstream os;
      
      os << "VRR Part" << part << " " <<
      (where == InBra ? "bra" : "ket") << " ( ";
      F sh_a(target_->bra(0,0)); os << sh_a.label() << " ";
      F sh_b(target_->ket(0,0)); os << sh_b.label() << " | R12^" << (K < 0 ? "m" : "") << abs(K) << " * G12 | ";
      F sh_c(target_->bra(1,0)); os << sh_c.label() << " ";
      F sh_d(target_->ket(1,0)); os << sh_d.label() << " )";
      
      return os.str();
    }
    
  /*
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,0,InBra> VRR_a_11_TwoPRep_11_sh;
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,1,InBra> VRR_c_11_TwoPRep_11_sh;
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,0,InKet> VRR_b_11_TwoPRep_11_sh;
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,1,InKet> VRR_d_11_TwoPRep_11_sh;
  */

};

#endif
