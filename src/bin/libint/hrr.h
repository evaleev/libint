
#ifndef _libint2_src_bin_libint_hrr_h_
#define _libint2_src_bin_libint_hrr_h_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <rr.h>
#include <integral.h>

using namespace std;


namespace libint2 {

  /** A generic Horizontal Recurrence Relation:

  |a b) = |a+1 b-1) + AB |a b-1)

  Int is the integral class. part specifies for which particle
  the angular momentum is shifted. Function a is assumed to gain quanta,
  function b loses quanta. loc_a and loc_b specify where
  functions a and b are located (bra or ket). pos_a and pos_b
  specify which function to be used (usually pos_a and pos_b are set
  to 0 to refer to the first function for this particle in this location).

*/
  template <template <class> class I, class BFSet, int part,
  FunctionPosition loc_a, unsigned int pos_a,
  FunctionPosition loc_b, unsigned int pos_b>
  class HRR : public RecurrenceRelation {

  public:
    typedef I<BFSet> TargetType;
    typedef I<BFSet> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef AlgebraicOperator<DGVertex> ExprType;

    /**
      dir specifies which quantum number of a and b is shifted.
      For example, dir can be 0 (x), 1(y), or 2(z) if F is
      a Cartesian Gaussian.
      */
    HRR(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~HRR();

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return nchildren_; };
    /// Implementation of RecurrenceRelation::num_expr()
    const unsigned int num_expr() const { return nexpr_; };
    /// returns pointer to the target
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer to i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// expr(i) returns pointer to the expression for the i-th child
    SafePtr<ExprType> expr(unsigned int i) const;
    /// Implementation of RecurrenceRelation's target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation's child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    /// Implementation of RecurrenceRelation::rr_expr()
    SafePtr<DGVertex> rr_expr(unsigned int i) const { return static_pointer_cast<DGVertex,ExprType>(expr(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }
    /// Implementation of RecurrenceRelation::label()
    std::string label() const { return label_; }
    /// Implementation of RecurrenceRelation::nflops()
    unsigned int nflops() const { return nflops_; }

    const std::string cpp_function_name() {}
    const std::string cpp_source_name() {}
    const std::string cpp_header_name() {}
    std::ostream& cpp_source(std::ostream&) {}

  private:
    static const unsigned int max_nchildren_ = 2;
    static const unsigned int max_nexpr_ = 2;
    unsigned int dir_;

    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[max_nchildren_];
    SafePtr<ExprType> expr_[max_nexpr_];

    unsigned int nchildren_;
    unsigned int nexpr_;
    unsigned int nflops_;
    void oper_checks() const;

    std::string label_;
    std::string generate_label(const SafePtr<TargetType>& target) const;
  };

  
  
  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::HRR(const SafePtr<TargetType>& Tint, unsigned int dir) :
    target_(Tint), dir_(dir), nchildren_(0), nexpr_(0), nflops_(0), label_(generate_label(Tint))
    {
      target_ = Tint;
      typename I<F>::AuxQuantaType aux = Tint->aux();

      typedef typename I<F>::BraType IBraType;
      typedef typename I<F>::KetType IKetType;
      IBraType* bra = new IBraType(Tint->bra());
      IKetType* ket = new IKetType(Tint->ket());

      //
      // InBra and InKet cases have to treated explicitly since BraType and KetType don't have to match
      //
      if (loc_b == InBra) {
        // See if b-1 exists
        F sh_b(bra->member(part,pos_b));
        try {
          sh_b.dec(dir_);
        }
        catch (InvalidDecrement) {
          delete bra;
          delete ket;
          return;
        }
        bra->set_member(sh_b,part,pos_b);
        children_[1] = I<F>::Instance(*bra,*ket,aux);

        if (loc_a == InBra) {  // a in bra
          F sh_a(bra->member(part,pos_a));
          sh_a.inc(dir_);
          bra->set_member(sh_a,part,pos_a);
        }
        else {  // a in ket
          F sh_a(ket->member(part,pos_a));
          sh_a.inc(dir_);
          ket->set_member(sh_a,part,pos_a);
        }
        children_[0] = I<F>::Instance(*bra,*ket,aux);
        nchildren_ += 2;

        if (is_simple()) {
          SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[1],children_[0]));
          SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.X_Y[part][dir],children_[1]));
          if (loc_a == InBra && loc_b == InKet) {
            SafePtr<ExprType> sum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr0_ptr,expr1_ptr));
            expr_[0] = sum_ptr;
            nexpr_ += 1;
          }
          else
            throw std::runtime_error("HRR::HRR() -- geometric prefactor is not general enough. Please, contact main developer.");
        }
      }
      else {
        // See if b-1 exists
        F sh_b(ket->member(part,pos_b));
        try {
          sh_b.dec(dir_);
        }
        catch (InvalidDecrement) {
          delete bra;
          delete ket;
          return;
        }
        ket->set_member(sh_b,part,pos_b);
        children_[1] = I<F>::Instance(*bra,*ket,aux);

        if (loc_a == InBra) {  // a in bra
          F sh_a(bra->member(part,pos_a));
          sh_a.inc(dir_);
          bra->set_member(sh_a,part,pos_a);
        }
        else {  // a in ket
          F sh_a(ket->member(part,pos_a));
          sh_a.inc(dir_);
          ket->set_member(sh_a,part,pos_a);
        }
        children_[0] = I<F>::Instance(*bra,*ket,aux);
        nchildren_ += 2;

        if (is_simple()) {
          SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[1],children_[0]));
          SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.X_Y[part][dir],children_[1]));
          if (loc_a == InBra && loc_b == InKet) {
            SafePtr<ExprType> sum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr0_ptr,expr1_ptr));
            expr_[0] = sum_ptr;
            nexpr_ += 1;
          }
          else
            throw std::runtime_error("HRR::HRR() -- geometric prefactor is not general enough. Please, contact main developer.");
        }
      }

      delete bra;
      delete ket;
    }

  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::~HRR()
    {
      oper_checks();
    }

  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    void
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::oper_checks() const
    {
      //
      // Here we check basic HRR applicability requirements on the integral class
      //

      // part is within the range
      typedef typename I<F>::OperatorType Oper;
      if (part < 0 || part >= Oper::Properties::np) {
        assert(false);
      }

      // can move across operator only if it's multiplicative
      if (loc_a != loc_b && !Oper::Properties::multiplicative) {
        assert(false);
      }

      // Cannot apply when a and b are the same
      if (loc_a == loc_b && pos_a == pos_b) {
        assert(false);
      }
    }
          
  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    SafePtr<typename HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::ChildType>
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::child(unsigned int i) const
    {
      assert(i>=0 && i<nchildren_);

      unsigned int nc=0;
      for(int c=0; c<max_nchildren_; c++) {
        if (children_[c]) {
          if (nc == i)
            return children_[c];
          nc++;
        }
      }
    };

  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    SafePtr< AlgebraicOperator<DGVertex> >
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::expr(unsigned int i) const
    {
      assert(i>=0 && i<nexpr_);

      unsigned int ne=0;
      for(int e=0; e<max_nexpr_; e++) {
        if (expr_[e]) {
          if (ne == i)
            return expr_[e];
          ne++;
        }
      }
    };

  template <template <class> class I, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    std::string
    HRR<I,F,part,loc_a,pos_a,loc_b,pos_b>::generate_label(const SafePtr<TargetType>& target) const
    {
      ostringstream os;
      
      os << "HRR Part " << part << " "
      << (loc_a == InBra ? "bra" : "ket") << " " << pos_a << "  "
      << (loc_b == InBra ? "bra" : "ket") << " " << pos_b << " ";
      
      if (loc_a == InBra) {
        F sh_a(target->bra(part,pos_a));
        os << sh_a.label() << " ";
        
        if (loc_b == InBra) {
          F sh_b(target->bra(part,pos_b));
          os << sh_b.label();
        }
        else {
          F sh_b(target->ket(part,pos_b));
          os << sh_b.label();
        }
      }
      else {
        F sh_a(target->ket(part,pos_a));
        os << sh_a.label() << " ";
        
        if (loc_b == InBra) {
          F sh_b(target->bra(part,pos_b));
          os << sh_b.label();
        }
        else {
          F sh_b(target->ket(part,pos_b));
          os << sh_b.label();
        }
      }
      
      return os.str();
    }
    
  typedef HRR<TwoPRep_11_11,CGShell,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,0,InKet,0,InBra,0> HRR_ba_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11,CGShell,1,InKet,0,InBra,0> HRR_dc_11_TwoPRep_11_sh;

  typedef HRR<TwoPRep_11_11,CGF,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_int;
  typedef HRR<TwoPRep_11_11,CGF,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_int;

};

#endif
