
#ifndef _libint2_src_bin_libint_vrr11twoprep11_h_
#define _libint2_src_bin_libint_vrr11twoprep11_h_

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>
#include <algebra.h>
#include <flop.h>
#include <prefactors.h>

using namespace std;


namespace libint2 {

  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
  the angular momentum is raised. bool bra specifies whether the angular momentum
  is raised in bra (true) or ket (false). Class ERI specifies which particular implementation
  of ERI to use.
  */
  template <template <class> class ERI, class BFSet, int part, FunctionPosition where>
  class VRR_11_TwoPRep_11 : public RecurrenceRelation {

  public:
    typedef ERI<BFSet> TargetType;
    typedef ERI<BFSet> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef AlgebraicOperator<DGVertex> ExprType;

    /**
      dir specifies which quantum number is incremented.
      For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
      a Cartesian Gaussian.
     */
    VRR_11_TwoPRep_11(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~VRR_11_TwoPRep_11();

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return nchildren_; };
    /// Implementation of RecurrenceRelation::num_expr()
    const unsigned int num_expr() const { return nexpr_; };
    /// target() returns pointer to the i-th child
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer to the i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// expr(i) returns pointer to the expression for the i-th child
    SafePtr<ExprType> expr(unsigned int i) const;
    /// Implementation of RecurrenceRelation::rr_target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation::rr_child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    /// Implementation of RecurrenceRelation::rr_expr()
    SafePtr<DGVertex> rr_expr(unsigned int i) const { return static_pointer_cast<DGVertex,ExprType>(expr(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }

    /// Implementation of RecurrenceRelation::nflops()
    unsigned int nflops() const { return nflops_; }
    
    const std::string cpp_function_name() {}
    const std::string cpp_source_name() {}
    const std::string cpp_header_name() {}
    std::ostream& cpp_source(std::ostream&) {}

  private:
    static const unsigned int max_nchildren_ = 8;
    static const unsigned int max_nexpr_ = 6;
    unsigned int dir_;

    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[max_nchildren_];
    SafePtr<ExprType> expr_[max_nexpr_];

    unsigned int nchildren_;
    unsigned int nexpr_;
    unsigned int nflops_;

  };
  
  template <template <class> class ERI, class F, int part, FunctionPosition where>
    VRR_11_TwoPRep_11<ERI,F,part,where>::VRR_11_TwoPRep_11(const SafePtr<ERI<F> >& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir)
    {
      target_ = Tint;

      F sh_a(Tint->bra(0,0));
      F sh_b(Tint->ket(0,0));
      F sh_c(Tint->bra(1,0));
      F sh_d(Tint->ket(1,0));
      unsigned int m = Tint->m();

      vector<F> bra;
      vector<F> ket;

      bra.push_back(sh_a);
      bra.push_back(sh_c);
      ket.push_back(sh_b);
      ket.push_back(sh_d);

      // Zero out children pointers
      nchildren_ = 0;

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
      try {
        bra_ref->operator[](p_a).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[0] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
      children_[1] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      nchildren_ += 2;
      nflops_ += ConvertNumFlops<F>(3);
      if (is_simple()) {
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.XY_X[part][InKet][dir],children_[0]));
        SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.W_XY[part][dir],children_[1]));
        expr_[0] = expr0_ptr;
        expr_[1] = expr1_ptr;
      }
      nexpr_ += 2;

      // See if a-2 exists
      bool a_minus_2_exists = true;
      try {
        bra_ref->operator[](p_a).dec(dir);
      }
      catch (InvalidDecrement) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        children_[2] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
        children_[3] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
        const unsigned int ni_a = bra_ref->operator[](p_a).qn(dir);
        bra_ref->operator[](p_a).inc(dir);
        nchildren_ += 2;
        nflops_ += ConvertNumFlops<F>(5);
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.rho_o_alpha12[part], children_[3]));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Minus, children_[2], expr_intmd0));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Times, prefactors.one_o_2alpha12[part], expr_intmd1));
          SafePtr<ExprType> expr2_ptr(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_a], expr_intmd2));
          expr_[2] = expr2_ptr;
        }
        nexpr_ += 1;
      }

      // See if b-1 exists
      bool b_minus_1_exists = true;
      try {
        ket_ref->operator[](p_a).dec(dir);
      }
      catch (InvalidDecrement) {
        b_minus_1_exists = false;
      }
      if (b_minus_1_exists) {
        children_[4] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
        children_[5] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
        const unsigned int ni_b = ket_ref->operator[](p_a).qn(dir);
        ket_ref->operator[](p_a).inc(dir);
        nchildren_ += 2;
        nflops_ += ConvertNumFlops<F>(5);
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.rho_o_alpha12[part], children_[5]));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Minus, children_[4], expr_intmd0));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Times, prefactors.one_o_2alpha12[part], expr_intmd1));
          SafePtr<ExprType> expr3_ptr(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_b], expr_intmd2));
          expr_[3] = expr3_ptr;
        }
        nexpr_ += 1;
      }

      // See if c-1 exists
      try {
        bra_ref->operator[](p_c).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[6] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      const unsigned int ni_c = bra_ref->operator[](p_c).qn(dir);
      bra_ref->operator[](p_c).inc(dir);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(3);
      if (is_simple()) {
        SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.one_o_2alphasum, children_[6]));
        SafePtr<ExprType> expr4_ptr(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_c], expr_intmd0));
        expr_[4] = expr4_ptr;
      }
      nexpr_ += 1;

      // See if d-1 exists
      try {
        ket_ref->operator[](p_c).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[7] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      const unsigned int ni_d = ket_ref->operator[](p_c).qn(dir);
      ket_ref->operator[](p_c).inc(dir);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(3);
      if (is_simple()) {
        SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.one_o_2alphasum, children_[7]));
        SafePtr<ExprType> expr5_ptr(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_d], expr_intmd0));
        expr_[5] = expr5_ptr;
      }
      nexpr_ += 1;

    };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    VRR_11_TwoPRep_11<ERI,F,part,where>::~VRR_11_TwoPRep_11()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    SafePtr< ERI<F> >
    VRR_11_TwoPRep_11<ERI,F,part,where>::child(unsigned int i) const
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

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    SafePtr< AlgebraicOperator<DGVertex> >
    VRR_11_TwoPRep_11<ERI,F,part,where>::expr(unsigned int i) const
    {
      assert(i>=0 && i<nexpr_);
      unsigned int ne=0;
      for(int e=0; e<max_nexpr_; e++) {
        if (expr_[e] != 0) {
          if (ne == i)
            return expr_[e];
          ne++;
        }
      }
    };

  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InBra> VRR_a_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InBra> VRR_c_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InKet> VRR_b_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InKet> VRR_d_11_TwoPRep_11_sh;


};

#endif
