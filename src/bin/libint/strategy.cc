
#include <vector>
#include <algorithm>
#include <strategy.h>
#include <dg.h>
#include <rr.h>
#include <rr.templ.h>
#include <comp_11_tig12_11.h>
#include <graph_registry.h>

#include <master.h>

#if 1
#include <boost/mpl/list.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>
#endif

#define USE_HRR 1
#define LOCAL_DEBUG 0

using namespace std;
using namespace libint2;
namespace libint2 {
  
  template <class T> struct MasterStrategy;
  template <> struct MasterStrategy<TwoPRep_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_TwoPRep_11_sh,
    HRR_cd_11_TwoPRep_11_sh,
    VRR_a_11_TwoPRep_11_sh,
    VRR_c_11_TwoPRep_11_sh
    > value;
  };
  template <> struct MasterStrategy<TwoPRep_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_TwoPRep_11_int,
    HRR_cd_11_TwoPRep_11_int,
    VRR_a_11_TwoPRep_11_int,
    VRR_c_11_TwoPRep_11_int
    > value;
  };
  template <> struct MasterStrategy<R12kG12_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_R12kG12_11_sh,
    HRR_cd_11_R12kG12_11_sh,
    VRR_a_11_R12kG12_11_sh,
    VRR_c_11_R12kG12_11_sh
    > value;
  };
  template <> struct MasterStrategy<R12kG12_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_R12kG12_11_int,
    HRR_cd_11_R12kG12_11_int,
    VRR_a_11_R12kG12_11_int,
    VRR_c_11_R12kG12_11_int
    > value;
  };
  template <> struct MasterStrategy<DummySymmIntegral_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_Dummy_11_sh,
    HRR_cd_11_Dummy_11_sh
    > value;
  };
  template <> struct MasterStrategy<DummySymmIntegral_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_Dummy_11_int,
    HRR_cd_11_Dummy_11_int
    > value;
  };

  /// transform<RRType> encapsulates RRType and the action associated with RRType
  /// it's used by apply_strategy::operator()<RRType> via mpl::for_each
  template <class RRType>
  class apply_strategy_transform {
    public:
    typedef typename RRType::TargetType IntType;
    // return true if no more visitations necessary
    static bool visit(const SafePtr<DirectedGraph>& dg,
                      const SafePtr<IntType>& integral,
                      const SafePtr<Tactic>& tactic,
                      SafePtr<RecurrenceRelation>& rr,
                      Tactic::rr_stack& rrstack) {
      if (boost::is_same<typename IntType::BasisFunctionType,CGF>::value)
        return _visit_cgf(dg,integral,tactic,rr,rrstack);
      else
        return _visit_cgshell(dg,integral,tactic,rr,rrstack);
    }
    private:
      static bool _visit_cgshell(const SafePtr<DirectedGraph>& dg,
                          const SafePtr<IntType>& integral,
                          const SafePtr<Tactic>& tactic,
                          SafePtr<RecurrenceRelation>& rr,
                          Tactic::rr_stack& rrstack) {
        SafePtr<RRType> rr_ptr = RRType::Instance(integral, 0);
        // only first RR is needed in CGShell case
        if (rr_ptr->num_children()) {
          rr = static_pointer_cast<RecurrenceRelation, RRType>(rr_ptr);
          return true;
        }
        return false;
      }
      static bool _visit_cgf(const SafePtr<DirectedGraph>& dg,
                      const SafePtr<IntType>& integral,
                      const SafePtr<Tactic>& tactic,
                      SafePtr<RecurrenceRelation>& rr,
                      Tactic::rr_stack& rrstack) {
        // in CGF case collect all rrs on rrstack
        for (int xyz = 2; xyz >= 0; xyz--) {
          SafePtr<RRType> rr_ptr = RRType::Instance(integral, xyz);
          if (rr_ptr->num_children())
            rrstack.push_back(static_pointer_cast<RecurrenceRelation, RRType>(rr_ptr));
        }
        return false;
      }

    
  };

  /** This type helps with processing lists of integral types via mpl::for_each.
      It will attemp to cast integral to each T. For the first such match it will
      use Strategy<T> (a list of types also) to determine the optimal recurrence relation.

      This design follows section 9.1.2 of "C++ Template Metaprogramming" by Abrahams and Gurtovoy.
    */
  template <class IntType> class apply_strategy {
    public:

      struct Impl {
        Impl(const SafePtr<DirectedGraph>& dg,
             const SafePtr<IntType>& integral,
             const SafePtr<Tactic>& tactic) :
               dg_(dg), integral_(integral), tactic_(tactic), done_(false) {
            }
            
            const SafePtr<RecurrenceRelation>& rr() {
              // if rr() is called then we should no longer do any processing
              done_ = true;
              // some postprocessing may be needed to determine optimal rr
              if (boost::is_same<typename IntType::BasisFunctionType, CGF>::value)
                postprocess_rr_cgf(tactic_, rr_, rrstack_);
              return rr_;
            }
            
            SafePtr<DirectedGraph> dg_;
            SafePtr<IntType> integral_;
            SafePtr<Tactic> tactic_;
            SafePtr<RecurrenceRelation> rr_;
            bool done_;
            Tactic::rr_stack rrstack_;

            // for shell quartets the relation has been determined, for CGF use rrstack and tactic
            void postprocess_rr_cgf(const SafePtr<Tactic>& tactic,
                                    SafePtr<RecurrenceRelation>& rr,
                                    const Tactic::rr_stack& rrstack) {
              rr = tactic->optimal_rr(rrstack);
            }
            
        };
      
    apply_strategy(const SafePtr<Impl>& impl) : impl_(impl) {}
    apply_strategy(const apply_strategy& app) : impl_(app.impl_) {}
    const apply_strategy& operator=(const apply_strategy& app) {
      impl_ = app.impl_;
      return *this;
    }
    
    template <class Visitor>
    void operator()(const Visitor&) {
      if (!impl_->done_)
        impl_->done_ = Visitor::visit(impl_->dg_, impl_->integral_, impl_->tactic_, impl_->rr_, impl_->rrstack_);
    }
    
    const SafePtr<Impl>& impl() const { return impl_; }
    
    private:
      SafePtr<Impl> impl_;
  };

  /// transform<T> encapsulates T and the action associated with T
  /// it's used by operator()<T> via mpl::for_each
  template <class T>
  struct match_first_inttype_transform {
    
    static bool visit(const SafePtr<DirectedGraph>& dg,
                      const SafePtr<DGVertex>& integral,
                      const SafePtr<Tactic>& tactic,
                      SafePtr<RecurrenceRelation>& rr) {
      SafePtr<T> tptr = dynamic_pointer_cast<T,DGVertex>(integral);
      //std::cout << "In match_inttype_first::visit() integral=" << integral->label() << endl;
      if (tptr != 0) {
        using namespace boost;
        using namespace boost::mpl::placeholders;
        // Try to unroll first, if this is a shell set
        const unsigned int size = integral->size();
        const bool can_unroll = boost::is_same<typename T::BasisFunctionType,CGShell>::value &&
                                (size <= dg->registry()->unroll_threshold());
        if (can_unroll) {
          typedef IntegralSet_to_Integrals<T> ISet2I;
          SafePtr<ISet2I> x(new ISet2I(tptr));
          rr = static_pointer_cast<RecurrenceRelation,ISet2I>(x);
        }
        // else: apply the known strategy
        else {
          typedef apply_strategy<T> apply_strategy_t;
          typedef typename apply_strategy_t::Impl apply_strategy_t_impl;
          SafePtr<apply_strategy_t_impl> applier_impl(new apply_strategy_t_impl(dg,tptr,tactic));
          apply_strategy_t applier(applier_impl);
          mpl::for_each<typename MasterStrategy<T>::value, apply_strategy_transform<_1>, apply_strategy_t& >(applier);
          rr = applier_impl->rr();
        }
        return true;
      }
    }
    
  };

  /** This type helps with processing lists of integral types via mpl::for_each.
      It will attemp to cast integral to each T. For the first such match it will
      use Strategy<T> (a list of types also) to determine the optimal recurrence relation.

      This design follows section 9.1.2 of "C++ Template Metaprogramming" by Abrahams and Gurtovoy.
    */
  class match_first_inttype {
  public:
    struct Impl {
      Impl(const SafePtr<DirectedGraph>& dg,
           const SafePtr<DGVertex>& integral,
           const SafePtr<Tactic>& tactic) :
                   dg_(dg),
                   integral_(integral),
                   tactic_(tactic),
                   found_this_type_(false)
                   {
                   }

      const SafePtr<RecurrenceRelation>& rr() const { return rr_; }
      
      SafePtr<DirectedGraph> dg_;
      SafePtr<DGVertex> integral_;
      SafePtr<Tactic> tactic_;
      SafePtr<RecurrenceRelation> rr_;
      bool found_this_type_;
    };

    match_first_inttype(const SafePtr<Impl>& impl) : impl_(impl) {}
    match_first_inttype(const match_first_inttype& x) : impl_(x.impl_) {}
    const match_first_inttype& operator=(const match_first_inttype& x) {
      impl_ = x.impl_;
      return *this;
    }
    
    template <class Visitor>
    void operator()(const Visitor&) {
      if (!impl_->found_this_type_)
        impl_->found_this_type_ = Visitor::visit(impl_->dg_,impl_->integral_,impl_->tactic_,impl_->rr_);
    }

    const SafePtr<Impl>& impl() const { return impl_; }
    
    private:
      SafePtr<Impl> impl_;
};
  

}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr(const SafePtr<DirectedGraph>& graph,
                     const SafePtr<DGVertex>& integral,
                     const SafePtr<Tactic>& tactic)
{
#if 0
  {
    std::cout << "Strategy::optimal_rr() -- integral:" << std::endl;
    integral->print(std::cout);
  }
#endif

#if 1
  using namespace boost;
  using namespace boost::mpl::placeholders;

  // by iterate over the master typelist determine to the type of this integral
  // use the type to search through the type-specific strategy (see apply_strategy<T>)
  SafePtr<match_first_inttype::Impl> matcher_impl(new match_first_inttype::Impl(graph,integral,tactic));
  match_first_inttype matcher(matcher_impl);
  mpl::for_each<MasterIntegralTypeList, match_first_inttype_transform<_1>, match_first_inttype& >(matcher);
  return matcher_impl->rr();
#endif

#if 0
  //
  // We must first determine the type of the integral
  //
  {
    SafePtr<TwoPRep_11_11_sq> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_sq,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_sq(graph,eri_ptr,tactic);
  }
  {
    SafePtr<TwoPRep_11_11_int> eri_ptr = dynamic_pointer_cast<TwoPRep_11_11_int,DGVertex>(integral);
    if (eri_ptr != 0)
      return optimal_rr_twoprep1111_int(graph,eri_ptr,tactic);
  }
  {
    typedef R12kG12_11_11_sq IType;
    SafePtr<IType> bptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (bptr != 0) {
      return optimal_rr_R12kG121111_sq(graph,bptr,tactic);
    }
  }
  {
    typedef R12kG12_11_11_int IType;
    SafePtr<IType> bptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (bptr != 0) {
      return optimal_rr_R12kG121111_int(graph,bptr,tactic);
    }
  }
#if 0
  {
    typedef TiG12_11_11_base<CGShell> base_type;
    SafePtr<base_type> bptr = dynamic_pointer_cast<base_type,DGVertex>(integral);
    if (bptr != 0) {
      int k = TiG12_11_11_Util::i<CGShell>(bptr);
      switch (k) {
        case 0:
          return optimal_rr_TiG121111_sq<0>(graph,bptr,tactic);
        case 1:
          return optimal_rr_TiG121111_sq<1>(graph,bptr,tactic);
        default:
          throw logic_error("Strategy::optimal_rr() unable to determine K for TiG12_11_11<CGShell,K> class");
      };
    }
  }
  {
    typedef TiG12_11_11_base<CGF> base_type;
    SafePtr<base_type> bptr = dynamic_pointer_cast<base_type,DGVertex>(integral);
    if (bptr != 0) {
      int k = TiG12_11_11_Util::i<CGF>(bptr);
      switch (k) {
        case 0:
          return optimal_rr_TiG121111_int<0>(graph,bptr,tactic);
        case 1:
          return optimal_rr_TiG121111_int<1>(graph,bptr,tactic);
        default:
          throw logic_error("Strategy::optimal_rr() unable to determine K for TiG12_11_11<CGF,K> class");
      };
    }
  }
  {
    typedef R1dotR1G12_11_11_sq IType;
    SafePtr<IType> iptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (iptr != 0)
      return optimal_rr_R1dotR1G121111_sq(graph,iptr,tactic);
  }
  {
    typedef R1dotR1G12_11_11_int IType;
    SafePtr<IType> iptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (iptr != 0)
      return optimal_rr_R1dotR1G121111_int(graph,iptr,tactic);
  }
  {
    typedef R2dotR2G12_11_11_sq IType;
    SafePtr<IType> iptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (iptr != 0)
      return optimal_rr_R2dotR2G121111_sq(graph,iptr,tactic);
  }
  {
    typedef R2dotR2G12_11_11_int IType;
    SafePtr<IType> iptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (iptr != 0)
      return optimal_rr_R2dotR2G121111_int(graph,iptr,tactic);
  }
  {
    typedef R1dotR2G12_11_11_sq IType;
    SafePtr<IType> iptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (iptr != 0)
      return optimal_rr_R1dotR2G121111_sq(graph,iptr,tactic);
  }
  {
    typedef R1dotR2G12_11_11_int IType;
    SafePtr<IType> iptr = dynamic_pointer_cast<IType,DGVertex>(integral);
    if (iptr != 0)
      return optimal_rr_R1dotR2G121111_int(graph,iptr,tactic);
  }
#endif
  
  // Type insensitive RR can be applied to almost any integral
  {
    typedef DummySymmIntegral_11_11_sq int_type;
    SafePtr<int_type> int_ptr = dynamic_pointer_cast<int_type,DGVertex>(integral);
    if (int_ptr != 0)
      return optimal_rr_Dummy1111_sq(graph,int_ptr,tactic);
  }
  {
    typedef DummySymmIntegral_11_11_int int_type;
    SafePtr<int_type> int_ptr = dynamic_pointer_cast<int_type,DGVertex>(integral);
    if (int_ptr != 0)
      return optimal_rr_Dummy1111_int(graph,int_ptr,tactic);
  }

  // Don't know how to apply any RR
  return SafePtr<RecurrenceRelation>();
#endif
}

#if 0
SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_sq(const SafePtr<DirectedGraph>& graph,
                                    const SafePtr<TwoPRep_11_11_sq>& integral,
                                    const SafePtr<Tactic>& tactic)
{
  static SafePtr<RecurrenceRelation> nullptr;
#if DEBUG
  std::cout << "Strategy::optimal_rr_twoprep1111_sq: called for " << integral->label() << std::endl;
#endif

  //
  // This is a basic strategy for computing integral
  // 1) first see if should convert the set to infividual integrals
  // 2) if possible apply HRR
  // 3) else apply VRR
  //
  const unsigned int size = integral->size();
  const bool can_unroll = graph->registry()->can_unroll();

  if (size == 1 || (can_unroll && size <= max_size_to_unroll_)) {
#if DEBUG
    std::cout << "Strategy::optimal_rr_twoprep1111_sq: " << integral->label() << " to be unrolled" << std::endl;
#endif
    return unroll_intset<TwoPRep_11_11_sq>(integral);
  }

#if LIBINT_ERI_STRATEGY == 1 || LIBINT_ERI_STRATEGY == 2
  //
  // Particle 0 is most significant for storage, hence want to perform HRR on it last,
  // when functions of particle 1 are ready. This will maximize the length of the loops
  // in HRRPart0... code.
  //

  // shift from B to A
  {
    typedef HRR_ab_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children()) {
#if LOCAL_DEBUG
      std::cout << "Applying HRR(ab) to ";  integral->print(std::cout);  std::cout << std::endl;
#endif
      return rr_cast(rr_ptr);
    }
  }

  // shift from D to C
  {
    typedef HRR_cd_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children()) {
#if LOCAL_DEBUG
      std::cout << "Applying HRR(cd) to ";  integral->print(std::cout);  std::cout << std::endl;
#endif
      return rr_cast(rr_ptr);
    }
  }
#endif
  
#if LIBINT_ERI_STRATEGY == 2
  {
    typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11,CGShell,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
#endif
  
  {
    typedef VRR_a_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children()) {
#if LOCAL_DEBUG
      std::cout << "Applying VRR(a) to ";  integral->print(std::cout);  std::cout << std::endl;
#endif
      return rr_cast(rr_ptr);
    }
  }

#if !USE_HRR
  {
    typedef VRR_b_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children()) {
#if LOCAL_DEBUG
      std::cout << "Applying VRR(b) to ";  integral->print(std::cout);  std::cout << std::endl;
#endif
      return rr_cast(rr_ptr);
    }
  }
#endif
  
  {
    typedef VRR_c_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children()) {
#if LOCAL_DEBUG
      std::cout << "Applying VRR(c) to ";  integral->print(std::cout);  std::cout << std::endl;
#endif
      return rr_cast(rr_ptr);
    }
  }
  
#if !USE_HRR
  {
    typedef VRR_d_11_TwoPRep_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children()) {
#if LOCAL_DEBUG
      std::cout << "Applying VRR(d) to ";  integral->print(std::cout);  std::cout << std::endl;
#endif
      return rr_cast(rr_ptr);
    }
  }
#endif
  
  return SafePtr<RecurrenceRelation>();
}



// Generate all possible recurrence relations and then
// use a Tactic object to decide which to use

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_twoprep1111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr<TwoPRep_11_11_int>& integral,
                                     const SafePtr<Tactic>& tactic)
{
  vector<RR> rrstack;  // stack of all recurrence relations
  
#if LIBINT_ERI_STRATEGY == 1 || LIBINT_ERI_STRATEGY == 2
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_ab_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_cd_11_TwoPRep_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
#endif

#if LIBINT_ERI_STRATEGY == 2
  // shift from A to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef ITR_11_TwoPRep_11<GenIntegralSet_11_11,CGF,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
#endif

  // decrease A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<GenIntegralSet_11_11,CGF,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  // decrease B
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<GenIntegralSet_11_11,CGF,0,InKet> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  // Decrease C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<GenIntegralSet_11_11,CGF,1,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // Decrease D
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef VRR_11_TwoPRep_11<GenIntegralSet_11_11,CGF,1,InKet> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
  
  return tactic->optimal_rr(rrstack);
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R12kG121111_sq(const SafePtr<DirectedGraph>& graph,
                                    const SafePtr< R12kG12_11_11_sq >& integral,
                                    const SafePtr<Tactic>& tactic)
{
  typedef R12kG12_11_11_sq inttype;
      
  //
  // This is a basic strategy for computing integral
  // 1) first see if should convert the SafePtr<RecurrenceRelation>
  //    set to infividual integrals
  // 2) if possible apply HRR
  // 3) else apply VRR
  //
  const unsigned int size = integral->size();
  if (size == 1 || (size <= max_size_to_unroll_ && graph->registry()->can_unroll()))
    return unroll_intset<inttype>(integral);
    
  {
    // AB HRR
    typedef HRR<inttype,CGShell,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
      
  {
    // CD HRR
    typedef HRR<inttype,CGShell,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
      
  {
    // A VRR
    typedef VRR_11_R12kG12_11<GenIntegralSet_11_11,CGShell,0,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }
      
  {
    // C VRR
    typedef VRR_11_R12kG12_11<GenIntegralSet_11_11,CGShell,1,InBra> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
   }
      
  return SafePtr<RecurrenceRelation>();
}

// Generate all possible recurrence relations and then
// use a Tactic object to decide which to use
SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R12kG121111_int(const SafePtr<DirectedGraph>& graph,
                                     const SafePtr< R12kG12_11_11_int >& integral,
                                     const SafePtr<Tactic>& tactic)
{
  typedef R12kG12_11_11_int inttype;
  vector<RR> rrstack;  // stack of all recurrence relations
  
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<inttype,CGF,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
      
  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<inttype,CGF,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // only apply VRR is AM on B and D is zero
#if USE_BRAKET_H
  if (integral->ket(0,0).zero() && integral->ket(1,0).zero()) {
#else
  if (integral->ket(0,0)->zero() && integral->ket(1,0)->zero()) {
#endif
    // decrease A
    for(int xyz = 2; xyz >= 0; xyz--) {
      typedef VRR_11_R12kG12_11<GenIntegralSet_11_11,CGF,0,InBra> rr_type;
      SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
      if (rr_ptr->num_children())
        rrstack.push_back(rr_cast(rr_ptr));
    }
        
    // Decrease C
    for(int xyz = 2; xyz >= 0; xyz--) {
      typedef VRR_11_R12kG12_11<GenIntegralSet_11_11,CGF,1,InBra> rr_type;
      SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
      if (rr_ptr->num_children())
        rrstack.push_back(rr_cast(rr_ptr));
    }
  }
      
  return tactic->optimal_rr(rrstack);
}

#if 0
SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R1dotR1G121111_sq(const SafePtr<DirectedGraph>& graph,
				       const SafePtr<R1dotR1G12_11_11_sq>& integral,
				       const SafePtr<Tactic>& tactic)
{
#define __IType_tmpl1 R1dotR1G12_11_11
  typedef R1dotR1G12_11_11_sq IType;

  const unsigned int size = integral->size();
  if (size == 1 || (graph->registry()->can_unroll() && size <= max_size_to_unroll_)) {
#if DEBUG
    std::cout << "Strategy::optimal_rr_r1dotr1g121111_sq: " << integral->label() << " to be unrolled" << std::endl;
#endif
    return unroll_intset<IType>(integral);
  }

  {
    // AB HRR
    typedef HRR<IType,IType::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    // CD HRR
    typedef HRR<IType,IType::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    // Build relation
    typedef CR_11_R1dotR1G12_11<__IType_tmpl1,IType::BasisFunctionType> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R1dotR1G121111_int(const SafePtr<DirectedGraph>& graph,
					const SafePtr<R1dotR1G12_11_11_int>& integral,
					const SafePtr<Tactic>& tactic)
{
#define __IType_tmpl2 R1dotR1G12_11_11
  typedef R1dotR1G12_11_11_int IType;
  vector<RR> rrstack;  // stack of all recurrence relations

  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<IType,IType::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<IType,IType::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  {
    // Build relation
    typedef CR_11_R1dotR1G12_11<__IType_tmpl2,IType::BasisFunctionType> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  return tactic->optimal_rr(rrstack);
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R2dotR2G121111_sq(const SafePtr<DirectedGraph>& graph,
				       const SafePtr<R2dotR2G12_11_11_sq>& integral,
				       const SafePtr<Tactic>& tactic)
{
#define __IType_tmpl5 R2dotR2G12_11_11
  typedef R2dotR2G12_11_11_sq IType;

  const unsigned int size = integral->size();
  if (size == 1 || (graph->registry()->can_unroll() && size <= max_size_to_unroll_)) {
#if DEBUG
    std::cout << "Strategy::optimal_rr_r2dotr2g121111_sq: " << integral->label() << " to be unrolled" << std::endl;
#endif
    return unroll_intset<IType>(integral);
  }

  {
    // AB HRR
    typedef HRR<IType,IType::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    // CD HRR
    typedef HRR<IType,IType::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    // Build relation
    typedef CR_11_R2dotR2G12_11<__IType_tmpl5,IType::BasisFunctionType> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R2dotR2G121111_int(const SafePtr<DirectedGraph>& graph,
					const SafePtr<R2dotR2G12_11_11_int>& integral,
					const SafePtr<Tactic>& tactic)
{
#define __IType_tmpl6 R2dotR2G12_11_11
  typedef R2dotR2G12_11_11_int IType;
  vector<RR> rrstack;  // stack of all recurrence relations

  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<IType,IType::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<IType,IType::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  {
    // Build relation
    typedef CR_11_R2dotR2G12_11<__IType_tmpl6,IType::BasisFunctionType> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  return tactic->optimal_rr(rrstack);
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R1dotR2G121111_sq(const SafePtr<DirectedGraph>& graph,
				       const SafePtr<R1dotR2G12_11_11_sq>& integral,
				       const SafePtr<Tactic>& tactic)
{
#define __IType_tmpl3 R1dotR2G12_11_11
  typedef R1dotR2G12_11_11_sq IType;

  const unsigned int size = integral->size();
  if (size == 1 || (graph->registry()->can_unroll() && size <= max_size_to_unroll_)) {
#if DEBUG
    std::cout << "Strategy::optimal_rr_r1dotr2g121111_sq: " << integral->label() << " to be unrolled" << std::endl;
#endif
    return unroll_intset<IType>(integral);
  }

  {
    // AB HRR
    typedef HRR<IType,IType::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    // CD HRR
    typedef HRR<IType,IType::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    // Build relation
    typedef CR_11_R1dotR2G12_11<__IType_tmpl3,IType::BasisFunctionType> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_R1dotR2G121111_int(const SafePtr<DirectedGraph>& graph,
					const SafePtr<R1dotR2G12_11_11_int>& integral,
					const SafePtr<Tactic>& tactic)
{
#define __IType_tmpl4 R1dotR2G12_11_11
  typedef R1dotR2G12_11_11_int IType;
  vector<RR> rrstack;  // stack of all recurrence relations

  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<IType,IType::BasisFunctionType,0,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR<IType,IType::BasisFunctionType,1,InBra,0,InKet,0> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  {
    // Build relation
    typedef CR_11_R1dotR2G12_11<__IType_tmpl4,IType::BasisFunctionType> rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  return tactic->optimal_rr(rrstack);
}
#endif

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_Dummy1111_sq(const SafePtr<DirectedGraph>& graph,
				  const SafePtr<DummySymmIntegral_11_11_sq>& integral,
				  const SafePtr<Tactic>& tactic)
{
  const unsigned int size = integral->size();
  if (size == 1 || (size <= max_size_to_unroll_ && graph->registry()->can_unroll()))
    return unroll_intset<DummySymmIntegral_11_11_sq>(integral);

  {
    typedef HRR_ab_11_Dummy_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  {
    typedef HRR_cd_11_Dummy_11_sh rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,0);
    if (rr_ptr->num_children())
      return rr_cast(rr_ptr);
  }

  return SafePtr<RecurrenceRelation>();
}

SafePtr<RecurrenceRelation>
Strategy::optimal_rr_Dummy1111_int(const SafePtr<DirectedGraph>& graph,
				   const SafePtr<DummySymmIntegral_11_11_int>& integral,
				   const SafePtr<Tactic>& tactic)
{
  vector<RR> rrstack;  // stack of all recurrence relations
  
#if USE_HRR
  // shift from B to A
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_ab_11_Dummy_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }

  // shift from D to C
  for(int xyz = 2; xyz >= 0; xyz--) {
    typedef HRR_cd_11_Dummy_11_int rr_type;
    SafePtr<rr_type> rr_ptr = rr_type::Instance(integral,xyz);
    if (rr_ptr->num_children())
      rrstack.push_back(rr_cast(rr_ptr));
  }
#endif

  return tactic->optimal_rr(rrstack);
}

#endif
