
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <assert.h>
#include <dgvertex.h>
#include <rr.h>
#include <integral.h>
#include <tig12_11_11.h>
#include <algebra.h>
#include <flop.h>
#include <prefactors.h>
#include <context.h>
#include <default_params.h>

#ifndef _libint2_src_bin_libint_cr11tig1211_h_
#define _libint2_src_bin_libint_cr11tig1211_h_

using namespace std;


namespace libint2 {

  /** Compute relation for 2-e integrals of the Ti_G12 operators.
  I<BFSet,K> is the integral set specialization that describes the
  integrals of the Ti_G12 operator.
  */
  template <template <class,int> class I, class BFSet, int K>
  class CR_11_TiG12_11 : public RecurrenceRelation {

  public:
    typedef I<BFSet,K> TargetType;
    typedef R12kG12_11_11_base<BFSet> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef AlgebraicOperator<DGVertex> ExprType;

    /**
      dir specifies which quantum number is incremented.
      For example, dir can be 0 (x), 1(y), or 2(z) if BFSet is
      a Cartesian Gaussian.
     */
    CR_11_TiG12_11(const SafePtr<TargetType>&);
    ~CR_11_TiG12_11();

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
    SafePtr<DGVertex> rr_child(unsigned int i) const { return dynamic_pointer_cast<DGVertex,ChildType>(child(i)); }
    /// Implementation of RecurrenceRelation::rr_expr()
    SafePtr<DGVertex> rr_expr(unsigned int i) const { return static_pointer_cast<DGVertex,ExprType>(expr(i)); }
    /// Implementation of RecurrenceRelation::is_simple()
    bool is_simple() const {
      return TrivialBFSet<BFSet>::result;
    }
    /// Implementation of RecurrenceRelation::invariant_type()
    bool invariant_type() const {
      return true;
    }
    /// Implementation of RecurrenceRelation::label()
    std::string label() const { return label_; }

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
    static const unsigned int max_nchildren_ = 18;
    static const unsigned int max_nexpr_ = 1;
    
    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[max_nchildren_];
    SafePtr<ExprType> expr_[max_nexpr_];

    unsigned int nchildren_;
    unsigned int nexpr_;
    unsigned int nflops_;
 
    std::string label_;
    std::string generate_label(const SafePtr<TargetType>& target) const;

    void add_expr(const SafePtr<ExprType>&, int minus=1);
  };
  
  template <template <class,int> class I, class F, int K>
    CR_11_TiG12_11<I,F,K>::CR_11_TiG12_11(const SafePtr<I<F,K> >& Tint) :
    target_(Tint), nchildren_(0), nexpr_(0), nflops_(0), label_(generate_label(Tint))
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

      // On which particle to act
      int p_a = K;
      int p_c = (p_a == 0) ? 1 : 0;

      for(int braket=0; braket<=1; braket++) {
        FunctionPosition where = (FunctionPosition)braket;

        // Use indirection to choose bra or ket
        vector<F>* bra_ref = &bra;
        vector<F>* ket_ref = &ket;
        if (where == InKet) {
          bra_ref = &ket;
          ket_ref = &bra;
        }

        typedef R12kG12_11_11<F,0> child_type;

        const unsigned int ndirs = is_simple() ? 3 : 1;
        for(int xyz=0; xyz<ndirs; xyz++) {
          
          bool am1_exists = true;
          bool am2_exists = true;
          try {
            bra_ref->operator[](p_a).dec(xyz);
          }
          catch (InvalidDecrement) {
            am1_exists = false;
            am2_exists = false;
          }
          
          if (am1_exists) {
            try {
              bra_ref->operator[](p_a).dec(xyz);
            }
            catch (InvalidDecrement) {
              am2_exists = false;
              bra_ref->operator[](p_a).inc(xyz);
            }
          }
          
          if (am2_exists) {
            int next_child = nchildren_;
            children_[next_child] = child_type::Instance(bra[0],ket[0],bra[1],ket[1],0);
            nchildren_ += 1;
            if (is_simple()) {
              const unsigned int ni_a = bra_ref->operator[](p_a).qn(xyz) + 2;
              SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[ni_a * (ni_a-1)],rr_child(next_child)));
              add_expr(expr0_ptr);
            }
            nflops_ += ConvertNumFlops<F>(1);
            bra_ref->operator[](p_a).inc(xyz);
            bra_ref->operator[](p_a).inc(xyz);
          }

          // a+2
          {
            bra_ref->operator[](p_a).inc(xyz);
            bra_ref->operator[](p_a).inc(xyz);
            int next_child = nchildren_;
            children_[next_child] = child_type::Instance(bra[0],ket[0],bra[1],ket[1],0);
            nchildren_ += 1;
            if (is_simple()) {
              const unsigned int ni_a = bra_ref->operator[](p_a).qn(xyz);
              unsigned int pfac = 2*(ni_a + 1);
              if (am1_exists)
                pfac += 2*ni_a;
              SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.zeta[K][where],prefactors.zeta[K][where]));
              SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.Cdouble(4.0),expr0_ptr));
              SafePtr<ExprType> expr2_ptr(new ExprType(ExprType::OperatorTypes::Times,expr1_ptr,rr_child(next_child)));
              add_expr(expr0_ptr);
            }
            nflops_ += ConvertNumFlops<F>(3);
            bra_ref->operator[](p_a).dec(xyz);
            bra_ref->operator[](p_a).dec(xyz);
          }

          // a
          {
            int next_child = nchildren_;
            children_[next_child] = child_type::Instance(bra[0],ket[0],bra[1],ket[1],0);
            nchildren_ += 1;
            if (is_simple()) {
              const unsigned int ni_a = bra_ref->operator[](p_a).qn(xyz);
              unsigned int pfac = 2*(ni_a + 1);
              if (am1_exists)
                pfac += 2*ni_a;
              SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[pfac],prefactors.zeta[K][where]));
              SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,expr0_ptr,rr_child(next_child)));
              add_expr(expr0_ptr,-1);
            }
            nflops_ += ConvertNumFlops<F>(1);
          }
        }
      }
      // scale by -0.5
      SafePtr<ExprType> scaled(new ExprType(ExprType::OperatorTypes::Times,prefactors.Cdouble(-0.5),expr_[0]));
      expr_[0] = scaled;
      nflops_ += ConvertNumFlops<F>(1);
    }
  
  template <template <class,int> class I, class F, int K>
    CR_11_TiG12_11<I,F,K>::~CR_11_TiG12_11()
    {
      if (K < 0 || K >= 2) {
        assert(false);
      }
    };

  template <template <class,int> class I, class F, int K>
    SafePtr< typename CR_11_TiG12_11<I,F,K>::ChildType >
    CR_11_TiG12_11<I,F,K>::child(unsigned int i) const
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

  template <template <class,int> class I, class F, int K>
    SafePtr< AlgebraicOperator<DGVertex> >
    CR_11_TiG12_11<I,F,K>::expr(unsigned int i) const
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
    }

  template <template <class,int> class I, class F, int K>
    std::string
    CR_11_TiG12_11<I,F,K>::generate_label(const SafePtr<TargetType>& target) const
    {
      ostringstream os;
      
      os << "RR ( ";
      F sh_a(target->bra(0,0)); os << sh_a.label() << " ";
      F sh_b(target->ket(0,0)); os << sh_b.label() << " | [T_" << K << ", G12] | ";
      F sh_c(target->bra(1,0)); os << sh_c.label() << " ";
      F sh_d(target->ket(1,0)); os << sh_d.label() << " )";
      
      return os.str();
    }
    
   template <template <class,int> class I, class F, int K>
    std::string
    CR_11_TiG12_11<I,F,K>::spfunction_call(
    const SafePtr<CodeContext>& context, const SafePtr<ImplicitDimensions>& dims) const
    {
      ostringstream os;
      os << context->label_to_name(label_to_funcname(label()))
         // First argument is the library object
         << "(libint, "
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
    

  template <template <class,int> class I, class F, int K>
    void
    CR_11_TiG12_11<I,F,K>::add_expr(const SafePtr<ExprType>& expr, int minus)
    {
      if (nexpr_ == 0) {
        if (minus != -1) {
          expr_[0] = expr;
          nexpr_ = 1;
        }
        else {
          SafePtr<ExprType> negative(new ExprType(ExprType::OperatorTypes::Times,expr,prefactors.Cdouble(-1.0)));
          expr_[0] = negative;
          nflops_ += ConvertNumFlops<F>(1);
        }          
      }
      else {
        if (minus != -1) {
          SafePtr<ExprType> sum(new ExprType(ExprType::OperatorTypes::Plus,expr,expr_[0]));
          expr_[0] = sum;
          nflops_ += ConvertNumFlops<F>(1);
        }
        else {
          SafePtr<ExprType> sum(new ExprType(ExprType::OperatorTypes::Minus,expr,expr_[0]));
          expr_[0] = sum;
          nflops_ += ConvertNumFlops<F>(1);
        }
      }
    }
    
  /*
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,0,InBra> VRR_a_11_TwoPRep_11_sh;
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,1,InBra> VRR_c_11_TwoPRep_11_sh;
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,0,InKet> VRR_b_11_TwoPRep_11_sh;
  typedef VRR_11_R12kG12_11<R12kG12_11_11,CGShell,1,InKet> VRR_d_11_TwoPRep_11_sh;
  */

};

#endif
