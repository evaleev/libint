/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#include <vector>
#include <algorithm>

#include <strategy.h>
#include <dg.h>
#include <rr.h>
#include <graph_registry.h>
#include <intset_to_ints.h>
#include <uncontract.h>
#include <singl_stack.h>

#include <master_ints_list.h>
#include <master_rrs_list.h>

// MPL is painful
#include <boost/mpl/list.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>

using namespace std;
using namespace libint2;

namespace libint2 {

  // Depending on the set of shell quartets, may need to adjust the strategy
#define LIBINT_SHELLQUARTET_STRATEGY_A0C0 1
#define LIBINT_SHELLQUARTET_STRATEGY_0B0D 2
#define LIBINT_SHELLQUARTET_STRATEGY LIBINT_SHELLQUARTET_STRATEGY_A0C0
#if LIBINT_SHELL_SET == LIBINT_SHELL_SET_ORCA
# undef LIBINT_SHELLQUARTET_STRATEGY
# define LIBINT_SHELLQUARTET_STRATEGY LIBINT_SHELLQUARTET_STRATEGY_0B0D
#endif

  //
  // Particle 0 is most significant for storage, hence want to perform HRR on it last,
  // when functions of particle 1 are ready. This will maximize the length of the loops
  // in HRRPart0... code.
  //

#if LIBINT_ERI_STRATEGY == 0
# error "Not all recurrence relations are implemented yet for pure OS scheme (have 5 minutes to fix this?)"
#endif
  template <class T> struct MasterStrategy;
#if LIBINT_SHELLQUARTET_STRATEGY == LIBINT_SHELLQUARTET_STRATEGY_A0C0
  template <> struct MasterStrategy<TwoPRep_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_TwoPRep_11_sh,
    HRR_cd_11_TwoPRep_11_sh,
    Deriv_a_11_TwoPRep_11_sh,
    Deriv_b_11_TwoPRep_11_sh,
    Deriv_c_11_TwoPRep_11_sh,
    Deriv_d_11_TwoPRep_11_sh,
#if LIBINT_ERI_STRATEGY == 2
    ITR_a_11_TwoPRep_11_sh,
    ITR_c_11_TwoPRep_11_sh,
#endif
    VRR_a_11_TwoPRep_11_sh,
    VRR_c_11_TwoPRep_11_sh
    > value;
  };
  template <> struct MasterStrategy<TwoPRep_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_TwoPRep_11_int,
    HRR_cd_11_TwoPRep_11_int,
    Deriv_a_11_TwoPRep_11_int,
    Deriv_b_11_TwoPRep_11_int,
    Deriv_c_11_TwoPRep_11_int,
    Deriv_d_11_TwoPRep_11_int,
#if LIBINT_ERI_STRATEGY == 2
    ITR_a_11_TwoPRep_11_int,
    ITR_c_11_TwoPRep_11_int,
#endif
    VRR_a_11_TwoPRep_11_int,
    VRR_c_11_TwoPRep_11_int
    > value;
  };
#else  // 0B0D strategy
  template <> struct MasterStrategy<TwoPRep_11_11_sq> {
    typedef mpl::list<
    HRR_ba_11_TwoPRep_11_sh,
    HRR_dc_11_TwoPRep_11_sh,
    Deriv_a_11_TwoPRep_11_sh,
    Deriv_b_11_TwoPRep_11_sh,
    Deriv_c_11_TwoPRep_11_sh,
    Deriv_d_11_TwoPRep_11_sh,
#if LIBINT_ERI_STRATEGY == 2
    ITR_b_11_TwoPRep_11_sh,
    ITR_d_11_TwoPRep_11_sh,
#endif
    VRR_b_11_TwoPRep_11_sh,
    VRR_d_11_TwoPRep_11_sh
    > value;
  };
  template <> struct MasterStrategy<TwoPRep_11_11_int> {
    typedef mpl::list<
    HRR_ba_11_TwoPRep_11_int,
    HRR_dc_11_TwoPRep_11_int,
    Deriv_a_11_TwoPRep_11_int,
    Deriv_b_11_TwoPRep_11_int,
    Deriv_c_11_TwoPRep_11_int,
    Deriv_d_11_TwoPRep_11_int,
#if LIBINT_ERI_STRATEGY == 2
    ITR_b_11_TwoPRep_11_int,
    ITR_d_11_TwoPRep_11_int,
#endif
    VRR_b_11_TwoPRep_11_int,
    VRR_d_11_TwoPRep_11_int
    > value;
  };
#endif

#if LIBINT_SHELLQUARTET_STRATEGY == LIBINT_SHELLQUARTET_STRATEGY_A0C0
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
#else // 0B0D strategy
  template <> struct MasterStrategy<R12kG12_11_11_sq> {
    typedef mpl::list<
    HRR_ba_11_R12kG12_11_sh,
    HRR_dc_11_R12kG12_11_sh,
    VRR_b_11_R12kG12_11_sh,
    VRR_d_11_R12kG12_11_sh
    > value;
  };
  template <> struct MasterStrategy<R12kG12_11_11_int> {
    typedef mpl::list<
    HRR_ba_11_R12kG12_11_int,
    HRR_dc_11_R12kG12_11_int,
    VRR_b_11_R12kG12_11_int,
    VRR_d_11_R12kG12_11_int
    > value;
  };
#endif

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
  template <> struct MasterStrategy<G12TiG12_11_11_sq> {
    typedef mpl::list<
    HRR_ab_11_G12TiG12_11_sh,
    HRR_cd_11_G12TiG12_11_sh,
    CR_11_G12TiG12_11_sh
    > value;
  };
  template <> struct MasterStrategy<G12TiG12_11_11_int> {
    typedef mpl::list<
    HRR_ab_11_G12TiG12_11_int,
    HRR_cd_11_G12TiG12_11_int,
    CR_11_G12TiG12_11_int
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
#if LIBINT_SHELLQUARTET_STRATEGY == LIBINT_SHELLQUARTET_STRATEGY_A0C0
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
#else // 0B0D strategy
  template <> struct MasterStrategy<DummySymmIntegral_11_11_sq> {
    typedef mpl::list<
    HRR_ba_11_Dummy_11_sh,
    HRR_dc_11_Dummy_11_sh
    > value;
  };
  template <> struct MasterStrategy<DummySymmIntegral_11_11_int> {
    typedef mpl::list<
    HRR_ba_11_Dummy_11_int,
    HRR_dc_11_Dummy_11_int
    > value;
  };
#endif

#if LIBINT_SUPPORT_ONEBODYINTS
  template <> struct MasterStrategy<Overlap_1_1_sh> {
      typedef mpl::list<
        CR_XYZ_1_1<CGShell,OverlapOper>
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_int> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF,OverlapOper>
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_sh_x> {
      typedef mpl::list<
        CR_XYZ_1_1<CGShell1d<CartesianAxis_X>,OverlapOper>
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_sh_y> {
      typedef mpl::list<
        CR_XYZ_1_1<CGShell1d<CartesianAxis_Y>,OverlapOper>
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_sh_z> {
      typedef mpl::list<
        CR_XYZ_1_1<CGShell1d<CartesianAxis_Z>,OverlapOper>
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_int_x> {
      typedef mpl::list<
      VRR_a_1_Overlap_1_int_x,
      VRR_b_1_Overlap_1_int_x
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_int_y> {
      typedef mpl::list<
      VRR_a_1_Overlap_1_int_y,
      VRR_b_1_Overlap_1_int_y
      > value;
    };
  template <> struct MasterStrategy<Overlap_1_1_int_z> {
      typedef mpl::list<
      VRR_a_1_Overlap_1_int_z,
      VRR_b_1_Overlap_1_int_z
      > value;
    };
# if LIBINT_SHELLQUARTET_STRATEGY == LIBINT_SHELLQUARTET_STRATEGY_A0C0
  template <> struct MasterStrategy<ElecPot_1_1_sh> {
      typedef mpl::list<
      HRR_ab_1_ElecPot_1_sh,
      VRR_a_1_ElecPot_1_sh
      > value;
    };
  template <> struct MasterStrategy<ElecPot_1_1_int> {
      typedef mpl::list<
      HRR_ab_1_ElecPot_1_int,
      VRR_a_1_ElecPot_1_int
      > value;
    };
# else //  // 0B0D strategy
  template <> struct MasterStrategy<ElecPot_1_1_sh> {
      typedef mpl::list<
      HRR_ba_1_ElecPot_1_sh,
      VRR_b_1_ElecPot_1_sh
      > value;
    };
  template <> struct MasterStrategy<ElecPot_1_1_int> {
      typedef mpl::list<
      HRR_ba_1_ElecPot_1_int,
      VRR_b_1_ElecPot_1_int
      > value;
    };
# endif // strategy: A0 or 0B
  template <> struct MasterStrategy<Kinetic_1_1_sh> {
      typedef mpl::list<
        CR_XYZ_1_1<CGShell,KineticOper>
      > value;
    };
  template <> struct MasterStrategy<Kinetic_1_1_int> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF,KineticOper>
      > value;
    };
  template <> struct MasterStrategy<Kinetic_1_1_int_x> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF1d<CartesianAxis_X>,KineticOper>
      > value;
    };
  template <> struct MasterStrategy<Kinetic_1_1_int_y> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF1d<CartesianAxis_Y>,KineticOper>
      > value;
    };
  template <> struct MasterStrategy<Kinetic_1_1_int_z> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF1d<CartesianAxis_Z>,KineticOper>
      > value;
    };
  template <> struct MasterStrategy<CMultipole_1_1_sh> {
      typedef mpl::list<
        CR_XYZ_1_1<CGShell,CartesianMultipoleOper<3u>>
      > value;
    };
  template <> struct MasterStrategy<CMultipole_1_1_int> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF,CartesianMultipoleOper<3u>>
      > value;
    };
  template <> struct MasterStrategy<CMultipole_1_1_int_x> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF1d<CartesianAxis_X>,CartesianMultipoleOper<1u>>
      > value;
    };
  template <> struct MasterStrategy<CMultipole_1_1_int_y> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF1d<CartesianAxis_Y>,CartesianMultipoleOper<1u>>
      > value;
    };
  template <> struct MasterStrategy<CMultipole_1_1_int_z> {
      typedef mpl::list<
        CR_XYZ_1_1<CGF1d<CartesianAxis_Z>,CartesianMultipoleOper<1u>>
      > value;
    };
#endif // LIBINT_SUPPORT_ONEBODYINTS

  template <typename T> struct MasterStrategy {
      typedef mpl::list<> value;
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
        const bool rr_is_directional = RRType::directional();
        for (int xyz = (rr_is_directional ? 2 : 0); xyz >= 0; xyz--) {
          SafePtr<RRType> rr_ptr = RRType::Instance(integral, xyz);
          if (rr_ptr != 0)
            rrstack.push_back(static_pointer_cast<RecurrenceRelation, RRType>(rr_ptr));
        }
        return false;
      }
      static bool _visit_cgf(const SafePtr<DirectedGraph>& dg,
                      const SafePtr<IntType>& integral,
                      const SafePtr<Tactic>& tactic,
                      SafePtr<RecurrenceRelation>& rr,
                      Tactic::rr_stack& rrstack) {
        // If given NullTactic -- skip
        SafePtr<NullTactic> ntactic = dynamic_pointer_cast<NullTactic,Tactic>(tactic);
        if (ntactic)
          return false;

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
              // determine optimal rr
              postprocess_rr(tactic_, rr_, rrstack_);
              return rr_;
            }

            SafePtr<DirectedGraph> dg_;
            SafePtr<IntType> integral_;
            SafePtr<Tactic> tactic_;
            SafePtr<RecurrenceRelation> rr_;
            bool done_;
            Tactic::rr_stack rrstack_;

            // determine the optimal RR to use given the rrstack and tactic
            void postprocess_rr(const SafePtr<Tactic>& tactic,
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
      if (tptr != 0) {
        //std::cout << "In match_inttype_first::visit() integral=" << integral->label() << std::endl;
        using namespace boost;
        using namespace boost::mpl::placeholders;
        // Try to unroll first, if this is a shell set
        const unsigned int size = integral->size();
        const bool can_unroll = not TrivialBFSet<typename T::BasisFunctionType>::result &&
                                (size <= dg->registry()->unroll_threshold());
//        std::cout << "  size=" << size << " can unroll? " << (can_unroll ? "yes" : "no") << std::endl;
//        if (not can_unroll) {
//          std::cout << "    unroll_threshold=" << dg->registry()->unroll_threshold()
//                    << " shell set? " << (not TrivialBFSet<typename T::BasisFunctionType>::result ? "yes" : "no")
//                    << std::endl;
//        }
#if 0
        // for now only allow unrolling in primitive-basis code
        // TODO solve the problem with allowing unrolling in contracted code:
        // (ss|ps) top(HRR)-level code unrolls the quartet to integrals, but these integrals
        // are contracted, hence their evaluation is deferred to the prerequsite step
        // when constructing prereq graph these integrals are added and assigned addresses in
        // arbitrary order; to avoid this "for now" do this hack
        // for a more sound solution see PrerequisitesExtractor in dg.cc, unfortunately it doesn't seem to fully work
        // right now I don't have time to mess with this anymore
        const bool can_uncontract = dg->registry()->uncontract();
        if (can_unroll && can_uncontract) {
#endif
        if (can_unroll) {
          typedef IntegralSet_to_Integrals<T> ISet2I;
          SafePtr<ISet2I> x(new ISet2I(tptr));
          rr = static_pointer_cast<RecurrenceRelation,ISet2I>(x);
#if DEBUG
          std::cout << "Unrolled " << tptr->label() << std::endl;
#endif
        }
        else {
          // if allowed to uncontract -- try that first
          const bool can_uncontract = dg->registry()->uncontract();
          if (can_uncontract) {
            typedef Uncontract_Integral<T> UncI;
            SafePtr<UncI> x(new UncI(tptr));
            rr = static_pointer_cast<RecurrenceRelation,UncI>(x);
            if (rr != 0) {
              if (rr->num_children() != 0) {
#if DEBUG
                std::cout << "Uncontracted " << tptr->label() << std::endl;
#endif
                return true;
              }
            }
          }

          // if uncontraction failed -- apply the known strategy
          typedef apply_strategy<T> apply_strategy_t;
          typedef typename apply_strategy_t::Impl apply_strategy_t_impl;
          SafePtr<apply_strategy_t_impl> applier_impl(new apply_strategy_t_impl(dg,tptr,tactic));
          apply_strategy_t applier(applier_impl);
          mpl::for_each<typename MasterStrategy<T>::value, apply_strategy_transform<_1>, apply_strategy_t& >(applier);
          rr = applier_impl->rr();
#if DEBUG
          if (rr != 0)
            std::cout << "Selected the following RR: " << rr->label() << std::endl;
#endif
        }
        return true;
      }
      return false;
    }

  };

  /** This type helps with processing lists of integral types via mpl::for_each.
      It will attempt to cast integral to each T. For the first such match it will
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

  // iterate over the master typelist to determine the type of this integral
  // matcher then uses the type to search through the type-specific strategy (see apply_strategy<T>)
  SafePtr<match_first_inttype::Impl> matcher_impl(new match_first_inttype::Impl(graph,integral,tactic));
  match_first_inttype matcher(matcher_impl);
  mpl::for_each<MasterIntegralTypeList, match_first_inttype_transform<_1>, match_first_inttype& >(matcher);
  return matcher_impl->rr();
}
