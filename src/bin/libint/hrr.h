
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
#include <dummyintegral.h>
#include <algebra.h>
#include <dgvertex.h>
#include <prefactors.h>
#include <default_params.h>
#include <dims.h>
#include <task.h>

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
  template <class IntType, class BFSet, int part,
  FunctionPosition loc_a, unsigned int pos_a,
  FunctionPosition loc_b, unsigned int pos_b>
  class HRR : public RecurrenceRelation
    {

  public:
    typedef RecurrenceRelation ParentType;
    typedef HRR<IntType,BFSet,part,loc_a,pos_a,loc_b,pos_b> ThisType;
    typedef IntType TargetType;
    typedef IntType ChildType;
    /// A short alias
    typedef RecurrenceRelation::ExprType ExprType;

    /** Use Instance() to obtain an instance of RR. This function is provided to avoid
        issues with getting a SafePtr from constructor (as needed for registry to work).

        dir specifies which quantum number of a and b is shifted.
        For example, dir can be 0 (x), 1(y), or 2(z) if F is
        a Cartesian Gaussian.
    */
    static SafePtr<ThisType> Instance(const SafePtr<TargetType>&, unsigned int dir = 0);
    ~HRR();

    /// Implementation of RecurrenceRelation::num_children()
    const unsigned int num_children() const { return nchildren_; };
    /// returns pointer to the target
    SafePtr<TargetType> target() const { return target_; };
    /// child(i) returns pointer to i-th child
    SafePtr<ChildType> child(unsigned int i) const;
    /// Implementation of RecurrenceRelation::target()
    SafePtr<DGVertex> rr_target() const { return static_pointer_cast<DGVertex,TargetType>(target()); }
    /// Implementation of RecurrenceRelation::child()
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
    /// Reimplementation of RecurrenceRelation::description()
    const std::string& description() const
      {
        if (descr_.empty()) {
          ostringstream oss;
          oss << label() << "   target = " << target_->label();
          descr_ = oss.str();
        }
        return descr_;
      }
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
      dir specifies which quantum number of a and b is shifted.
      For example, dir can be 0 (x), 1(y), or 2(z) if F is
      a Cartesian Gaussian.
      */
    HRR(const SafePtr<TargetType>&, unsigned int dir);

    /// registers this RR with the stack, if needed
    bool register_with_rrstack() const;

    static const unsigned int max_nchildren_ = 2;
    unsigned int dir_;

    SafePtr<TargetType> target_;
    SafePtr<ChildType> children_[max_nchildren_];
    SafePtr<ExprType> expr_;

    unsigned int nchildren_;
    unsigned int nflops_;
    void oper_checks() const;

    mutable std::string descr_;
    std::string label_;
    std::string generate_label(const SafePtr<TargetType>& target) const;
    /// Overload of RecurrenceRelation::adapt_dims_()
    SafePtr<ImplicitDimensions> adapt_dims_(const SafePtr<ImplicitDimensions>& dims) const;
    /** return true if the high dimension must be shown explicitly. For example,
        cd-HRR applied (ss|pp) has high dimension of rank 1 but since the code for such
        RR is not specific to ab=(ss|, the rank of high dimension must be shown explicitly.
      */
    bool expl_high_dim() const;
    bool expl_low_dim() const;
  };

  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    SafePtr< HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b> >
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::Instance(const SafePtr<TargetType>& Tint, unsigned int dir)
    {
      SafePtr<ThisType> this_ptr(new ThisType(Tint,dir));
      // Do post-construction duties
      if (this_ptr->num_children() != 0) {
        this_ptr->register_with_rrstack();
      }
      return this_ptr;
    }
  
  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::HRR(const SafePtr<TargetType>& Tint, unsigned int dir) :
    target_(Tint), dir_(dir), nchildren_(0), nflops_(0), label_(generate_label(Tint))
    {
      target_ = Tint;
      typename IntType::AuxQuantaType aux = Tint->aux();

      typedef typename IntType::BraType IBraType;
      typedef typename IntType::KetType IKetType;
      IBraType* bra = new IBraType(Tint->bra());
      IKetType* ket = new IKetType(Tint->ket());

      //
      // InBra and InKet cases have to treated explicitly since BraType and KetType don't have to match
      //
      if (loc_b == InBra) {
        // See if b-1 exists
        F sh_b(bra->member(part,pos_b));
	if (!sh_b.dec(dir_)) {
          delete bra;
          delete ket;
          return;
        }
        bra->set_member(sh_b,part,pos_b);
        children_[1] = IntType::Instance(*bra,*ket,aux);

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
        children_[0] = IntType::Instance(*bra,*ket,aux);
        nchildren_ += 2;

        if (is_simple()) {
          SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[1],children_[0]));
          SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.X_Y[part][dir],children_[1]));
          if (loc_a == InBra && loc_b == InKet) {
            SafePtr<ExprType> sum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr0_ptr,expr1_ptr));
            expr_ = sum_ptr;
          }
          else
            throw std::runtime_error("HRR::HRR() -- geometric prefactor is not general enough. Please, contact main developer.");
        }
      }
      else {
        // See if b-1 exists
        F sh_b(ket->member(part,pos_b));
	if (!sh_b.dec(dir_)) {
          delete bra;
          delete ket;
          return;
        }
        ket->set_member(sh_b,part,pos_b);
        children_[1] = IntType::Instance(*bra,*ket,aux);

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
        children_[0] = IntType::Instance(*bra,*ket,aux);
        nchildren_ += 2;

        if (is_simple()) {
          SafePtr<ExprType> expr0_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.N_i[1],children_[0]));
          SafePtr<ExprType> expr1_ptr(new ExprType(ExprType::OperatorTypes::Times,prefactors.X_Y[part][dir],children_[1]));
          if (loc_a == InBra && loc_b == InKet) {
            SafePtr<ExprType> sum_ptr(new ExprType(ExprType::OperatorTypes::Plus,expr0_ptr,expr1_ptr));
            expr_ = sum_ptr;
          }
          else
            throw std::runtime_error("HRR::HRR() -- geometric prefactor is not general enough. Please, contact main developer.");
        }
      }

      delete bra;
      delete ket;
    }

  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    bool
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::register_with_rrstack() const
    {
      // This is an ugly hack -- add HRR's to the RRStack preemptively for targets in which all functions not involved
      // in transfer have zero quanta. The reason is that for the HRR quartet level code to work correctly
      // I must use these particular instances of HRR to generate the source.
      
      // only register RRs with for shell sets
      if (TrivialBFSet<F>::result)
        return false;
      typedef typename IntType::BraType IBraType;
      typedef typename IntType::KetType IKetType;
      const IBraType& bra = target_->bra();
      const IKetType& ket = target_->ket();
      
      //check for nonzero quanta for all particles other than part
      bool nonzero_quanta = false;
      unsigned const int npart = IntType::OperatorType::Properties::np;
      for(int p=0; p<npart; p++) {
        if (p == part)
          continue;
        int nfbra = bra.num_members(p);
        for(int f=0; f<nfbra; f++)
#if USE_BRAKET_H
          if (!bra.member(p,f).zero())
#else
          if (!bra.member(p,f)->zero())
#endif
            nonzero_quanta = true;
        int nfket = ket.num_members(p);
        for(int f=0; f<nfket; f++)
#if USE_BRAKET_H
          if (!ket.member(p,f).zero())
#else
          if (!ket.member(p,f)->zero())
#endif
            nonzero_quanta = true;
      }
      // if all bfsets not involved in transfer have zero quanta then this instance needs to be added to the stack
      if (!nonzero_quanta) {
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
      
      //
      // else create the needed instance of HRR
      //
      
      // zero out unneeded bfs'
      IBraType bra_zero(bra);
      IKetType ket_zero(ket);
      for(int p=0; p<npart; p++) {
        if (p == part)
          continue;
        int nfbra = bra_zero.num_members(p);
        for(int f=0; f<nfbra; f++) {
          typedef typename IBraType::bfs_type bfs_type;
          typedef typename IBraType::bfs_ref bfs_ref;
          bfs_ref bfs = bra_zero.member(p,f);
#if USE_BRAKET_H
          if (!bfs.zero()) {
#else
          if (!bfs->zero()) {
#endif
            bfs_type null_bfs;
            swap(bfs,null_bfs);
          }
        }
        int nfket = ket_zero.num_members(p);
        for(int f=0; f<nfket; f++) {
          typedef typename IKetType::bfs_type bfs_type;
          typedef typename IKetType::bfs_ref bfs_ref;
          bfs_ref bfs = ket_zero.member(p,f);
#if USE_BRAKET_H
          if (!bfs.zero()) {
#else
          if (!bfs->zero()) {
#endif
            bfs_type null_bfs;
            swap(bfs,null_bfs);
          }
        }
      }
      // create a generic GenIntegralSet over a multiplicative operator
      typedef GenSymmOper< OperatorProperties<IntType::OperatorType::Properties::np,true,PermutationalSymmetry::symm> > DummyOper;
      typedef typename IBraType::bfs_type bfs_type;
      typedef EmptySet DummyQuanta;
      typedef GenIntegralSet<DummyOper, IncableBFSet, IBraType, IKetType, DummyQuanta> DummyIntegral;
      DummyOper dummy_oper;
      DummyQuanta dummy_quanta(std::vector<int>(0,0));
      SafePtr<DummyIntegral> dummy_integral = DummyIntegral::Instance(bra_zero,ket_zero,dummy_quanta,dummy_oper);
      
      // Construct generic HRR and add it to the stack instead of this HRR
      typedef HRR<DummyIntegral,F,part,loc_a,pos_a,loc_b,pos_b> DummyHRR;
      SafePtr<DummyHRR> dummy_hrr = DummyHRR::Instance(dummy_integral,dir_);
      SafePtr<RRStack> rrstack = RRStack::Instance();
      rrstack->find(dummy_hrr);
      return true;
      // not done coding -- throw an exception for now
      throw CodeDoesNotExist("Have not finished generalizing HRR::register with rrstack()");
    }


  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::~HRR()
    {
      oper_checks();
    }

  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    void
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::oper_checks() const
    {
      //
      // Here we check basic HRR applicability requirements on the integral class
      //

      // part is within the range
      typedef typename IntType::OperatorType Oper;
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
          
  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    SafePtr<typename HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::ChildType>
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::child(unsigned int i) const
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

  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    std::string
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::generate_label(const SafePtr<TargetType>& target) const
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
    
  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    std::string
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::spfunction_call(
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
      // then dimensions of basis function sets not involved in the transfer
      unsigned int hsr = 1;
      // a cleaner way to count the number of function sets referring
      // to some particles is to construct a dummy integral and
      // use subiterator policy
      // WARNING !!!
      for(int p=0; p<part; p++) {
        unsigned int nbra = target_->bra().num_members(p);
        for(int i=0; i<nbra; i++) {
          SubIterator* iter = target_->bra().member_subiter(p,i);
          hsr *= iter->num_iter();
	  delete iter;
        }
        unsigned int nket = target_->ket().num_members(p);
        for(int i=0; i<nket; i++) {
          SubIterator* iter = target_->ket().member_subiter(p,i);
          hsr *= iter->num_iter();
	  delete iter;
        }
      }
      // Use TaskParameters to keep track of maximum hsr
      LibraryTaskManager& taskmgr = LibraryTaskManager::Instance();
      taskmgr.current().params()->max_hrr_hsrank(hsr);
      
      // can only do a simple bra->ket or ket->bra transfer so far
      unsigned int isr = 1;
      if (loc_a == loc_b && pos_a != 0 && pos_b != 0)
        throw CodeDoesNotExist("HRR::spfunction_call -- has not been generalized yet");
      
      /// WARNING !!!
      unsigned int lsr = 1;
      unsigned int np = IntType::OperType::Properties::np;
      for(int p=part+1; p<np; p++) {
        unsigned int nbra = target_->bra().num_members(p);
        for(int i=0; i<nbra; i++) {
          SubIterator* iter = target_->bra().member_subiter(p,i);
          lsr *= iter->num_iter();
	  delete iter;
        }
        unsigned int nket = target_->ket().num_members(p);
        for(int i=0; i<nket; i++) {
          SubIterator* iter = target_->ket().member_subiter(p,i);
          lsr *= iter->num_iter();
	  delete iter;
        }
      }
      // Use TaskParameters to keep track of maximum hsr
      taskmgr.current().params()->max_hrr_hsrank(hsr);
      
      if (expl_high_dim())
        os << "," << hsr;
      if (expl_low_dim())
        os << "," << lsr;
      os << ")" << context->end_of_stat() << endl;
      return os.str();
    }
  
  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    bool
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::expl_high_dim() const
    {
      unsigned int np = IntType::OperType::Properties::np;
      bool high = true;
      if (part == 0)
        high = false;
      return high;
    }
  
  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    bool
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::expl_low_dim() const
    {
      unsigned int np = IntType::OperType::Properties::np;
      bool low = true;
      if (part == np -1)
        low = false;
      return low;
    }
  
  template <class IntType, class F, int part,
    FunctionPosition loc_a, unsigned int pos_a,
    FunctionPosition loc_b, unsigned int pos_b>
    SafePtr<ImplicitDimensions>
    HRR<IntType,F,part,loc_a,pos_a,loc_b,pos_b>::adapt_dims_(const SafePtr<ImplicitDimensions>& dims) const
    {
      bool high_rank = expl_high_dim();
      bool low_rank = expl_low_dim();

      SafePtr<Entity> high_dim, low_dim;
      if (high_rank) {
        high_dim = SafePtr<Entity>(new RTimeEntity<EntityTypes::Int>("highdim"));
      }
      else {
        high_dim = dims->high();
      }
      if (low_rank) {
        low_dim = SafePtr<Entity>(new RTimeEntity<EntityTypes::Int>("lowdim"));
      }
      else {
        low_dim = dims->low();
      }

      SafePtr<ImplicitDimensions> localdims(new ImplicitDimensions(high_dim,low_dim,dims->vecdim()));
      return localdims;
    }
  
  typedef HRR<TwoPRep_11_11<CGShell>,CGShell,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11<CGShell>,CGShell,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11<CGShell>,CGShell,0,InKet,0,InBra,0> HRR_ba_11_TwoPRep_11_sh;
  typedef HRR<TwoPRep_11_11<CGShell>,CGShell,1,InKet,0,InBra,0> HRR_dc_11_TwoPRep_11_sh;

  typedef HRR<TwoPRep_11_11<CGF>,CGF,0,InBra,0,InKet,0> HRR_ab_11_TwoPRep_11_int;
  typedef HRR<TwoPRep_11_11<CGF>,CGF,1,InBra,0,InKet,0> HRR_cd_11_TwoPRep_11_int;

  typedef HRR<DummySymmIntegral_11_11_sq,CGShell,0,InBra,0,InKet,0> HRR_ab_11_Dummy_11_sh;
  typedef HRR<DummySymmIntegral_11_11_sq,CGShell,1,InBra,0,InKet,0> HRR_cd_11_Dummy_11_sh;
  typedef HRR<DummySymmIntegral_11_11_int,CGF,0,InBra,0,InKet,0> HRR_ab_11_Dummy_11_int;
  typedef HRR<DummySymmIntegral_11_11_int,CGF,1,InBra,0,InKet,0> HRR_cd_11_Dummy_11_int;

};

#endif
