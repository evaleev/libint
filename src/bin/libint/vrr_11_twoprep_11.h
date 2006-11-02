
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

#ifndef _libint2_src_bin_libint_vrr11twoprep11_h_
#define _libint2_src_bin_libint_vrr11twoprep11_h_

using namespace std;


namespace libint2 {

  /** VRR Recurrence Relation for 2-e ERI. part specifies for which particle
  the angular momentum is raised. bool bra specifies whether the angular momentum
  is raised in bra (true) or ket (false). Class ERI specifies which particular implementation
  of ERI to use.
  */
  template <template <class> class ERI, class BFSet, int part, FunctionPosition where>
    class VRR_11_TwoPRep_11 : public RecurrenceRelation
    {

  public:
    typedef RecurrenceRelation ParentType;
    typedef VRR_11_TwoPRep_11 ThisType;
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
    ~VRR_11_TwoPRep_11();

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
    /// Implementation of RecurrenceRelation::spfunction_call()
    std::string spfunction_call(const SafePtr<CodeContext>& context,
                                const SafePtr<ImplicitDimensions>& dims) const;
    
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
    VRR_11_TwoPRep_11(const SafePtr<TargetType>&, unsigned int dir);

    /// registers this RR with the stack, if needed
    bool register_with_rrstack() const;

    static const unsigned int max_nchildren_ = 8;
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
    SafePtr< VRR_11_TwoPRep_11<ERI,F,part,where> >
    VRR_11_TwoPRep_11<ERI,F,part,where>::Instance(const SafePtr<TargetType>& Tint,
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
    VRR_11_TwoPRep_11<ERI,F,part,where>::register_with_rrstack() const
    {
      // only register RRs for shell sets
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
    VRR_11_TwoPRep_11<ERI,F,part,where>::VRR_11_TwoPRep_11(const SafePtr<ERI<F> >& Tint,
                                                           unsigned int dir) :
    target_(Tint), dir_(dir), nchildren_(0), nflops_(0), label_(generate_label(Tint))
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
        SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.XY_X[part][where][dir],children_[0]));
        SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.W_XY[part][dir],children_[1]));
        SafePtr<ExprType> expr0p1_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr0_ptr,expr1_ptr));
        expr_ = expr0p1_ptr;
      }

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
        bra_ref->operator[](p_a).inc(dir);
        const unsigned int ni_a = bra_ref->operator[](p_a).qn(dir);
        nchildren_ += 2;
        nflops_ += ConvertNumFlops<F>(5);
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.rho_o_alpha12[part], children_[3]));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Minus, children_[2], expr_intmd0));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_a], prefactors.one_o_2alpha12[part]));
          SafePtr<ExprType> expr2_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd2, expr_intmd1));
          SafePtr<ExprType> expr012_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr2_ptr,expr_));
          expr_ = expr012_ptr;
        }
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
        ket_ref->operator[](p_a).inc(dir);
        const unsigned int ni_b = ket_ref->operator[](p_a).qn(dir);
        nchildren_ += 2;
        nflops_ += ConvertNumFlops<F>(5);
        if (is_simple()) {
          SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.rho_o_alpha12[part], children_[5]));
          SafePtr<ExprType> expr_intmd1(new ExprType(ExprType::OperatorTypes::Minus, children_[4], expr_intmd0));
          SafePtr<ExprType> expr_intmd2(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_b], prefactors.one_o_2alpha12[part]));
          SafePtr<ExprType> expr3_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd2, expr_intmd1));
          SafePtr<ExprType> expr0123_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr3_ptr,expr_));
          expr_ = expr0123_ptr;
        }
      }

      // See if c-1 exists
      try {
        bra_ref->operator[](p_c).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[6] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      bra_ref->operator[](p_c).inc(dir);
      const unsigned int ni_c = bra_ref->operator[](p_c).qn(dir);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(3);
      if (is_simple()) {
        SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_c], prefactors.one_o_2alphasum));
        SafePtr<ExprType> expr4_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, children_[6]));
        SafePtr<ExprType> exprsum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr4_ptr,expr_));
        expr_ = exprsum_ptr;
      }

      // See if d-1 exists
      try {
        ket_ref->operator[](p_c).dec(dir);
      }
      catch (InvalidDecrement) {
        return;
      }
      children_[7] = ERI<F>::Instance(bra[0],ket[0],bra[1],ket[1],m+1);
      ket_ref->operator[](p_c).inc(dir);
      const unsigned int ni_d = ket_ref->operator[](p_c).qn(dir);
      nchildren_ += 1;
      nflops_ += ConvertNumFlops<F>(3);
      if (is_simple()) {
        SafePtr<ExprType> expr_intmd0(new ExprType(ExprType::OperatorTypes::Times, prefactors.N_i[ni_d], prefactors.one_o_2alphasum));
        SafePtr<ExprType> expr5_ptr(new ExprType(ExprType::OperatorTypes::Times, expr_intmd0, children_[7]));
        SafePtr<ExprType> exprsum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr5_ptr,expr_));
        expr_ = exprsum_ptr;
      }
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
    std::string
    VRR_11_TwoPRep_11<ERI,F,part,where>::generate_label(const SafePtr<TargetType>& target) const
    {
      ostringstream os;
      
      os << "OS VRR Part" << part << " " <<
      (where == InBra ? "bra" : "ket") << " ( ";
      F sh_a(target->bra(0,0)); os << sh_a.label() << " ";
      F sh_b(target->ket(0,0)); os << sh_b.label() << " | ";
      F sh_c(target->bra(1,0)); os << sh_c.label() << " ";
      F sh_d(target->ket(1,0)); os << sh_d.label() << " )";
      
      return os.str();
    }
    
   template <template <class> class ERI, class F, int part, FunctionPosition where>
    std::string
    VRR_11_TwoPRep_11<ERI,F,part,where>::spfunction_call(
    const SafePtr<CodeContext>& context, const SafePtr<ImplicitDimensions>& dims) const
    {
      ostringstream os;
      os << context->label_to_name(label_to_funcname(context->cparams()->api_prefix() + label()))
         // First argument is the library object
         << "(inteval, "
         // Second is the target
         << context->value_to_pointer(rr_target()->symbol());
      // then come children
      const unsigned int nchildren = num_children();
      for(int c=0; c<nchildren; c++) {
        os << ", " << context->value_to_pointer(rr_child(c)->symbol());
      }
      os << ")" << context->end_of_stat() << endl;
      return os.str();
    }
    
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InBra> VRR_a_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InBra> VRR_c_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,0,InKet> VRR_b_11_TwoPRep_11_sh;
  typedef VRR_11_TwoPRep_11<TwoPRep_11_11,CGShell,1,InKet> VRR_d_11_TwoPRep_11_sh;


};

#endif
