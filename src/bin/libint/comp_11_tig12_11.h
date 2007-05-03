
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
  class CR_11_TiG12_11 : public RecurrenceRelation
    {

  public:
    typedef RecurrenceRelation ParentType;
    typedef BFSet BasisFunctionType;
    typedef CR_11_TiG12_11<I,BFSet,K> ThisType;
    typedef I<BFSet,K> TargetType;
    typedef R12kG12_11_11<BFSet,0> ChildType;
    /// The type of expressions in which RecurrenceRelations result.
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&);
    ~CR_11_TiG12_11() {
      if (K < 0 || K >= 2) {
        assert(false);
      }
    }

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
    CR_11_TiG12_11(const SafePtr<TargetType>&);

    SafePtr<TargetType> target_;
    static const unsigned int max_nchildren_ = 18;
    SafePtr<ChildType> children_[max_nchildren_];
    unsigned int nchildren_;

    std::string generate_label() const
    {
      ostringstream os;
      os << "RR ( " << rr_target()->label() << " )";
      return os.str();
    }
  };

  template <template <class,int> class I, class F, int K>
    SafePtr< CR_11_TiG12_11<I,F,K> >
    CR_11_TiG12_11<I,F,K>::Instance(const SafePtr<TargetType>& Tint)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack<ThisType>();
      }
      return this_ptr;
    }

  template <template <class,int> class I, class F, int K>
    CR_11_TiG12_11<I,F,K>::CR_11_TiG12_11(const SafePtr<I<F,K> >& Tint) :
    target_(Tint), nchildren_(0)
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

        const unsigned int ndirs = is_simple() ? 3 : 1;
        for(int xyz=0; xyz<ndirs; xyz++) {
          
          bool am1_exists = true;
          bool am2_exists = true;
	  if (!bra_ref->operator[](p_a).dec(xyz)) {
            am1_exists = false;
            am2_exists = false;
          }
          
          if (am1_exists) {
	    if (!bra_ref->operator[](p_a).dec(xyz)) {
              am2_exists = false;
	      // return to a
              bra_ref->operator[](p_a).inc(xyz);
            }
          }
          
          if (am2_exists) {
            int next_child = nchildren_;
            children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0);
            nchildren_ += 1;
            if (is_simple()) {
              const unsigned int ni_a = bra_ref->operator[](p_a).qn(xyz) + 2;
              SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[ni_a * (ni_a-1)],rr_child(next_child)));
              if (where == InBra)
                add_expr(expr0_ptr);
              else
                add_expr(expr0_ptr,-1);
	      nflops_ += (1);
            }
            bra_ref->operator[](p_a).inc(xyz);
            bra_ref->operator[](p_a).inc(xyz);
          }

          // a+2
          {
            bra_ref->operator[](p_a).inc(xyz);
            bra_ref->operator[](p_a).inc(xyz);
            int next_child = nchildren_;
            children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0);
            nchildren_ += 1;
            if (is_simple()) {
              SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.Cdouble(4.0),prefactors.zeta2[K][where]));
              SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,expr0_ptr,rr_child(next_child)));
              if (where == InBra)
                add_expr(expr1_ptr);
              else
                add_expr(expr1_ptr,-1);
	      nflops_ += (3);
            }
            bra_ref->operator[](p_a).dec(xyz);
            bra_ref->operator[](p_a).dec(xyz);
          }

        }
      }

      // now add a
      {
        int next_child = nchildren_;
        children_[next_child] = ChildType::Instance(bra[0],ket[0],bra[1],ket[1],0);
        nchildren_ += 1;
        if (is_simple()) {
          // prefactor in front of (ab|cd) is (4*l_a+6)*zeta_a - (4*l_b+6)*zeta_b
          unsigned int l_x[2];
          for(int braket=0; braket<=1; braket++) {
            FunctionPosition where = (FunctionPosition)braket;
            l_x[braket] = 0;
            for(int xyz=0; xyz<3; xyz++) {
              if (where == InBra)
                l_x[braket] += bra[p_a].qn(xyz);
              else
                l_x[braket] += ket[p_a].qn(xyz);
            }
          }
          SafePtr<ExprType> pfaca_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.Cdouble(4*l_x[InBra]+6),prefactors.zeta[K][InBra]));
          SafePtr<ExprType> pfacb_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.Cdouble(4*l_x[InKet]+6),prefactors.zeta[K][InKet]));
          SafePtr<ExprType> pfac_ptr(new ExprType(ExprType::OperatorTypes::Minus,pfaca_ptr,pfacb_ptr));
          SafePtr<ExprType> expr_ptr(new ExprType(ExprType::OperatorTypes::Times,pfac_ptr,rr_child(next_child)));
          add_expr(expr_ptr,-1);
	  nflops_ += (4);
        }
      }

      if (is_simple()) {
	// scale by -0.5
	SafePtr<ExprType> scaled(new ExprType(ExprType::OperatorTypes::Times,prefactors.Cdouble(-0.5),expr_));
	expr_ = scaled;
	nflops_ += (1);
      }
    }
  
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

};

#endif
