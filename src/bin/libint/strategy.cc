
//
// ... likely the most impressive piece of C++ code I've ever written ...
//                    ------ ALL HAIL MPL ------

#include <vector>
#include <algorithm>
#include <strategy.h>
#include <dg.h>
#include <rr.h>
#include <rr.templ.h>
#include <graph_registry.h>
#include <intset_to_ints.h>
#include <singl_stack.timpl.h>

#include <master_ints_list.h>
#include <master_rrs_list.h>

#include <boost/mpl/list.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>

using namespace std;
using namespace libint2;

namespace libint2 {

  //
  // Particle 0 is most significant for storage, hence want to perform HRR on it last,
  // when functions of particle 1 are ready. This will maximize the length of the loops
  // in HRRPart0... code.
  //

#if LIBINT_ERI_STRATEGY == 0
# error "Not all recurrence relations are implemented yet for pure OS scheme (have 5 minutes to fix this?)"
#endif
  template <class T> struct MasterStrategy;
  template <> struct MasterStrategy<TwoPRep_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_TwoPRep_11_sh,
    HRR_cd_11_TwoPRep_11_sh,
#if LIBINT_ERI_STRATEGY == 2
    ITR_a_11_TwoPRep_11_sh,
#endif
    VRR_a_11_TwoPRep_11_sh,
    VRR_c_11_TwoPRep_11_sh
    > value;
  };
  template <> struct MasterStrategy<TwoPRep_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_TwoPRep_11_int,
    HRR_cd_11_TwoPRep_11_int,
#if LIBINT_ERI_STRATEGY == 2
    ITR_a_11_TwoPRep_11_int,
#endif
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
  template <> struct MasterStrategy<R12kR12lG12_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_R12kR12lG12_11_sh,
    HRR_cd_11_R12kR12lG12_11_sh,
    CR_11_R12kR12lG12_11_sh
    > value;
  };
  template <> struct MasterStrategy<R12kR12lG12_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_R12kR12lG12_11_int,
    HRR_cd_11_R12kR12lG12_11_int,
    CR_11_R12kR12lG12_11_int
    > value;
  };
  template <> struct MasterStrategy<TiG12_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_TiG12_11_sh,
    HRR_cd_11_TiG12_11_sh,
    CR_11_TiG12_11_sh
    > value;
  };
  template <> struct MasterStrategy<TiG12_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_TiG12_11_int,
    HRR_cd_11_TiG12_11_int,
    CR_11_TiG12_11_int
    > value;
  };
  template <> struct MasterStrategy<DivG12prime_xTx_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_DivG12prime_xTx_sh,
    HRR_cd_11_DivG12prime_xTx_sh,
    CR_11_DivG12prime_xTx_11_sh
    > value;
  };
  template <> struct MasterStrategy<DivG12prime_xTx_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_DivG12prime_xTx_int,
    HRR_cd_11_DivG12prime_xTx_int,
    CR_11_DivG12prime_xTx_11_int
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
        if (rr_ptr != 0) {
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
          // TODO: can I use the knowledge of Tactic behavior to skip some iteration?
          if (rr_ptr != 0)
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
#if DEBUG
          std::cout << "Unrolled " << tptr->label() << std::endl;
#endif
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

  using namespace boost;
  using namespace boost::mpl::placeholders;

  // by iterate over the master typelist determine to the type of this integral
  // matcher then uses the type to search through the type-specific strategy (see apply_strategy<T>)
  SafePtr<match_first_inttype::Impl> matcher_impl(new match_first_inttype::Impl(graph,integral,tactic));
  match_first_inttype matcher(matcher_impl);
  mpl::for_each<MasterIntegralTypeList, match_first_inttype_transform<_1>, match_first_inttype& >(matcher);
  return matcher_impl->rr();
}
