
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

#ifndef _libint2_src_bin_libint_itr11twoprep11_h_
#define _libint2_src_bin_libint_itr11twoprep11_h_

using namespace std;


namespace libint2 {

  /** ITR (Interelectron Transfer Relation) for 2-e ERI. part specifies for which particle
  the angular momentum is raised. where specifies whether the angular momentum
  is shifted in bra or ket. Class ERI specifies which particular implementation
  of ERI to use.
  */
  template <template <class> class ERI, class BFSet, int part, FunctionPosition where>
    class ITR_11_TwoPRep_11 : public RecurrenceRelation
    {

  public:
    typedef RecurrenceRelation ParentType;
    typedef ITR_11_TwoPRep_11 ThisType;
    typedef ERI<BFSet> TargetType;
    typedef ERI<BFSet> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).

        dir specifies which quantum number of a and b is shifted.
        For example, dir can be 0 (x), 1(y), or 2(z) if F is
        a Cartesian Gaussian.
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~ITR_11_TwoPRep_11();

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return nchildren_; };
    /// target() returns pointer to the i-th child
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer to the i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation::rr_target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation::rr_child()
    SafePtr<DGVertex> rr_child(unsigned int i) const { return static_pointer_cast<DGVertex,ChildType>(child(i)); }
    /// Implementation of RecurrenceRelation::rr_expr()
    SafePtr<ExprType> rr_expr() const { return expr_; }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }
    /// Implementation of RecurrenceRelation::invariant_type()
    bool invariant_type() const {
      return true;
    }
    /// Implementation of RecurrenceRelation::label()
    const std::string& label() const { return label_; }

    /// Implementation of RecurrenceRelation::nflops()
    unsigned int nflops() const { return nflops_; }
    
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
    ITR_11_TwoPRep_11(const SafePtr<TargetType>&, unsigned int dir);

    /// registers this RR with the stack, if needed
    bool register_with_rrstack() const;

    static const unsigned int max_nchildren_ = 4;
    unsigned int dir_;

    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[max_nchildren_];
    SafePtr<ExprType> expr_;

    unsigned int nchildren_;
    unsigned int nflops_;

    std::string label_;
    std::string generate_label(const SafePtr<TargetType>& target) const;
  };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    SafePtr< ITR_11_TwoPRep_11<ERI,F,part,where> >
    ITR_11_TwoPRep_11<ERI,F,part,where>::Instance(const SafePtr<TargetType>& Tint,
                                                  unsigned int dir)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint,dir));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack();
      }
      return this_ptr;
    }
  
  template <template <class> class ERI, class F, int part, FunctionPosition where>
    bool
    ITR_11_TwoPRep_11<ERI,F,part,where>::register_with_rrstack() const
    {
      // only register RRs with for shell sets
      if (TrivialBFSet<F>::result)
        return false;
      SafePtr<RRStack> rrstack = RRStack::Instance();
      SafePtr<ThisType> this_ptr =
	const_pointer_cast<ThisType,const ThisType>(
	  static_pointer_cast<const ThisType, const ParentType>(
	    EnableSafePtrFromThis<ParentType>::SafePtr_from_this()
	  )
	);
      rrstack->find(this_ptr);
      return true;
    }
  
  
  template <template <class> class ERI, class F, int part, FunctionPosition where>
    ITR_11_TwoPRep_11<ERI,F,part,where>::ITR_11_TwoPRep_11(const SafePtr<ERI<F> >& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir), nchildren_(0), nflops_(0), label_(generate_label(Tint))
    {
      /// InKet
      if (where == InKet)
        throw ProgrammingError("ITR_11_TwoPRep_11<ERI,F,part,where>::ITR_11_TwoPRep_11() -- where=InKet not implementd yet");

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
      children_[0] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(1);
      if (is_simple()) {
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.TwoPRepITR_pfac0[part][dir],children_[0]));
        expr_ = expr0_ptr;
      }

      // c+1
      bra_ref->operator[](p_c).inc(dir);
      children_[1] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
      bra_ref->operator[](p_c).dec(dir);
      const unsigned int ni_c = bra_ref->operator[](p_c).qn(dir);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(2);
      if (is_simple()) {
        SafePtr<ExprType> expr_ptr(new ExprType(ExprType::OperatorTypes::Times, prefactors.TwoPRepITR_pfac1[part], children_[1]));
        SafePtr<ExprType> exprsum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr_ptr,expr_));
        expr_ = exprsum_ptr;
      }

      // See if a-2 exists
      bool a_minus_2_exists = true;
      if (!bra_ref->operator[](p_a).dec(dir)) {
        a_minus_2_exists = false;
      }
      if (a_minus_2_exists) {
        children_[2] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
        bra_ref->operator[](p_a).inc(dir);
        const unsigned int ni_a = bra_ref->operator[](p_a).qn(dir);
        nchildren_ += 1;
        nflops_ += ConvertNumFlops<F>(3);
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_a], prefactors.one_o_2alpha12[part]));
          SafePtr<ExprType> expr_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, children_[2]));
          SafePtr<ExprType> exprsum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr_ptr,expr_));
          expr_ = exprsum_ptr;
        }
      }

      // See if c-1 exists
      if (!bra_ref->operator[](p_c).dec(dir)) {
        return;
      }
      children_[3] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m);
      bra_ref->operator[](p_c).inc(dir);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(3);
      if (is_simple()) {
        SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_c], prefactors.one_o_2alpha12[part]));
        SafePtr<ExprType> expr_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, children_[3]));
        SafePtr<ExprType> exprsum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr_ptr,expr_));
        expr_ = exprsum_ptr;
      }
    };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    ITR_11_TwoPRep_11<ERI,F,part,where>::~ITR_11_TwoPRep_11()
    {
      if (part < 0 || part >= 2) {
        assert(false);
      }
    };

  template <template <class> class ERI, class F, int part, FunctionPosition where>
    SafePtr< ERI<F> >
    ITR_11_TwoPRep_11<ERI,F,part,where>::child(unsigned int i) const
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
    std::string
    ITR_11_TwoPRep_11<ERI,F,part,where>::generate_label(const SafePtr<TargetType>& target) const
    {
      ostringstream os;
      
      os << "TwoPRep ITR Part" << part << " " <<
      (where == InBra ? "bra" : "ket") << " ( ";
      F sh_a(target->bra(0,0)); os << sh_a.label() << " ";
      F sh_b(target->ket(0,0)); os << sh_b.label() << " | ";
      F sh_c(target->bra(1,0)); os << sh_c.label() << " ";
      F sh_d(target->ket(1,0)); os << sh_d.label() << " )";
      
      return os.str();
    }
    
  typedef ITR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InBra> ITR_a_11_TwoPRep_11_sh;
  typedef ITR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InBra> ITR_c_11_TwoPRep_11_sh;
  typedef ITR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InKet> ITR_b_11_TwoPRep_11_sh;
  typedef ITR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InKet> ITR_d_11_TwoPRep_11_sh;


};

#endif
