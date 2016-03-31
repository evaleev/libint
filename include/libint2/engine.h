/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_lib_libint_engine_h_
#define _libint2_src_lib_libint_engine_h_

#include <libint2/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "libint2/engine.h requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>
#include <map>
#include <memory>

#include <libint2.h>
#include <libint2/boys.h>
#include <libint2/shell.h>
#include <libint2/timer.h>
#include <libint2/solidharmonics.h>
#include <libint2/any.h>
#include <libint2/util/compressed_pair.h>
#include <libint2/cxxapi.h>

#include <libint2/boost/preprocessor.hpp>

// extra PP macros

#define BOOST_PP_MAKE_TUPLE_INTERNAL(z, i, last) i BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i,last))
/// BOOST_PP_MAKE_TUPLE(n) returns (0,1,....n-1,n)
#define BOOST_PP_MAKE_TUPLE(n) ( BOOST_PP_REPEAT( n , BOOST_PP_MAKE_TUPLE_INTERNAL, BOOST_PP_DEC(n) ) )

#include <Eigen/Core>

// the engine will be profiled by default if library was configured with --enable-profile
#ifdef LIBINT2_PROFILE
#  define LIBINT2_ENGINE_TIMERS
// uncomment if want to profile each integral class
#  define LIBINT2_ENGINE_PROFILE_CLASS
#endif
// uncomment if want to profile the engine even if library was configured without --enable-profile
//#  define LIBINT2_ENGINE_TIMERS

namespace libint2 {

  constexpr size_t num_geometrical_derivatives(size_t ncenter,
                                               size_t deriv_order) {
    return (deriv_order > 0) ? (num_geometrical_derivatives(ncenter, deriv_order-1) * (3*ncenter+deriv_order-1))/deriv_order : 1;
  }

  template <typename T, unsigned N>
  typename std::remove_all_extents<T>::type *
  to_ptr1(T (&a)[N])  { return reinterpret_cast<
                                 typename std::remove_all_extents<T>::type *
                               >(&a); }

#if defined(LIBINT2_SUPPORT_ONEBODY)

  /**
   * OneBodyEngine computes integrals of operators (or operator sets) given by OneBodyEngine::operator_type
   */
  class OneBodyEngine {

    private:
      typedef struct {} empty_pod;

    public:

      /// types of operators (operator sets) supported by OneBodyEngine.
      /// \warning These must start with 0 and appear in same order as elements of BOOST_PP_ONEBODY_OPERATOR_LIST preprocessor macro.
      enum operator_type {
        overlap=0,        //!< overlap
        kinetic=1,        //!< electronic kinetic energy, i.e. \f$ -\frac{1}{2} \Nabla^2 \f$
        nuclear=2,        //!< Coulomb potential due to point charges
        emultipole1=3,    //!< overlap + (Cartesian) electric dipole moment, \f$ x_O, y_O, z_O \f$, where \f$ x_O \equiv x - O_x \f$ is relative to origin \f$ \vec{O} \f$
        emultipole2=4,    //!< emultipole1 + (Cartesian) electric quadrupole moment, \f$ x^2, xy, xz, y^2, yz, z^2 \f$
        emultipole3=5,    //!< emultipole2 + (Cartesian) electric octupole moment, \f$ x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3 \f$
        _invalid=-1
      };

/// list of libint task names for each operator type. These MUST appear in the same order as in operator_type
#define BOOST_PP_ONEBODY_OPERATOR_LIST  (overlap,         \
                                        (kinetic,         \
                                        (elecpot,         \
                                        (1emultipole,     \
                                        (2emultipole,     \
                                        (3emultipole,     \
                                         BOOST_PP_NIL))))))

#define BOOST_PP_ONEBODY_OPERATOR_INDEX_TUPLE BOOST_PP_MAKE_TUPLE( BOOST_PP_LIST_SIZE(BOOST_PP_ONEBODY_OPERATOR_LIST) )
#define BOOST_PP_ONEBODY_OPERATOR_INDEX_LIST BOOST_PP_TUPLE_TO_LIST( BOOST_PP_ONEBODY_OPERATOR_INDEX_TUPLE )

// make list of derivative orders for 1-body ints
#define BOOST_PP_ONEBODY_DERIV_ORDER_TUPLE BOOST_PP_MAKE_TUPLE( BOOST_PP_INC(LIBINT2_DERIV_ONEBODY_ORDER) )
#define BOOST_PP_ONEBODY_DERIV_ORDER_LIST BOOST_PP_TUPLE_TO_LIST( BOOST_PP_ONEBODY_DERIV_ORDER_TUPLE )

      /// alias to operator_type for backward compatibility to pre-05/13/2015 code
      /// \deprecated use operator_type instead
      DEPRECATED typedef operator_type integral_type;

      /// describes operator sets given by OneBodyOperator
      /// \note needs to be specialized for some operator types
      template <operator_type O> struct operator_traits;

      typedef libint2::FmEval_Taylor<real_t, 7> coulomb_core_eval_t;


      /// creates a default (unusable) OneBodyEngine; to be used as placeholder for copying a usable engine
      OneBodyEngine() : type_(_invalid), primdata_(), stack_size_(0), lmax_(-1) {}

      /// Constructs a (usable) OneBodyEngine

      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_order if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_order
      /// \param params a value of type OneBodyEngine::operator_traits<type>::oper_params_type specifying the parameters of
      ///               the operator set, e.g. position and magnitude of the charges creating the Coulomb potential
      ///               for type == nuclear. For most values of type this is not needed.
      ///               \sa OneBodyEngine::operator_traits
      template <typename Params = empty_pod>
      OneBodyEngine(operator_type type,
                    size_t max_nprim,
                    int max_l,
                    int deriv_order = 0,
                    Params params = empty_pod()) :
        type_(type),
        primdata_(max_nprim * max_nprim),
        lmax_(max_l),
        deriv_order_(deriv_order),
        params_(enforce_params_type(type,params)),
        fm_eval_(type == nuclear ? coulomb_core_eval_t::instance(2*max_l+deriv_order, 1e-25) : decltype(fm_eval_){})
      {
        initialize();
      }

      /// move constructor
      // intel does not accept "move ctor = default"
      OneBodyEngine(OneBodyEngine&& other) :
        type_(other.type_),
        primdata_(std::move(other.primdata_)),
        stack_size_(other.stack_size_),
        lmax_(other.lmax_),
        hard_lmax_(other.hard_lmax_),
        deriv_order_(other.deriv_order_),
        params_(std::move(other.params_)),
        fm_eval_(std::move(other.fm_eval_)),
        scratch_(std::move(other.scratch_)),
        scratch2_(other.scratch2_),
        buildfnptrs_(other.buildfnptrs_) {
      }

      /// (deep) copy constructor
      OneBodyEngine(const OneBodyEngine& other) :
        type_(other.type_),
        primdata_(other.primdata_.size()),
        stack_size_(other.stack_size_),
        lmax_(other.lmax_),
        deriv_order_(other.deriv_order_),
        params_(other.params_),
        fm_eval_(other.fm_eval_) {
        initialize();
      }

      ~OneBodyEngine() {
        finalize();
      }

      /// move assignment is default
      OneBodyEngine& operator=(OneBodyEngine&& other) {
        type_ = other.type_;
        primdata_ = std::move(other.primdata_);
        stack_size_ = other.stack_size_;
        lmax_ = other.lmax_;
        hard_lmax_ = other.hard_lmax_;
        deriv_order_ = other.deriv_order_;
        params_ = std::move(other.params_);
        fm_eval_ = std::move(other.fm_eval_);
        scratch_ = std::move(other.scratch_);
        scratch2_ = other.scratch2_;
        buildfnptrs_ = other.buildfnptrs_;
        return *this;
      }

      /// (deep) copy assignment
      OneBodyEngine& operator=(const OneBodyEngine& other) {
        type_ = other.type_;
        primdata_.resize(other.primdata_.size());
        stack_size_ = other.stack_size_;
        lmax_ = other.lmax_;
        deriv_order_ = other.deriv_order_;
        params_ = other.params_;
        fm_eval_ = other.fm_eval_;
        initialize();
        return *this;
      }

      /// resets operator parameters; this may be useful if need to compute Coulomb potential
      /// integrals over batches of charges for the sake of parallelism.
      template <typename Params>
      void set_params(const Params& params) {
        params_ = params;
        reset_scratch();
      }
      /// alias to set_params() for backward compatibility with pre-05/13/2015 code
      /// \deprecated use set_params() instead
      template <typename Params>
      DEPRECATED void set_q(const Params& params) {
        set_params(params);
      }

      /// reports the number of shell sets that each call to compute() produces.
      /// this depends on the order of geometrical derivatives requested and
      /// on the operator set. \sa compute()
      /// \note need to specialize for some operator types
      unsigned int nshellsets() const {
        const unsigned int num_operator_geometrical_derivatives = (type_ == nuclear) ? this->nparams() : 0;
        const auto ncenters = 2 + num_operator_geometrical_derivatives;
        return nopers() * num_geometrical_derivatives(ncenters,deriv_order_);
      }

      /// Given the Cartesian geometric derivative index that refers to center set
      /// (0...n-1) with one center omitted compute the derivative index referring
      /// to the full center set
      /// \param deriv_idx index of the derivative referring to the set with center \c omitted_center omitted
      /// \param deriv_order order of the geometric derivative
      /// \param ncenters number of centers in the full set
      /// \param omitted_center the omitted center, must be less than \c ncenters.
      /// \return the index of the derivative referring to full set
      static unsigned int to_target_deriv_index(unsigned int deriv_idx,
                                                unsigned int deriv_order,
                                                unsigned int ncenters,
                                                unsigned int omitted_center) {
        auto ncenters_reduced = ncenters - 1;
        auto nderiv_1d = 3 * ncenters_reduced;
        assert(deriv_idx < num_geometrical_derivatives(ncenters_reduced,deriv_order));
        switch (deriv_order) {
          case 1: return deriv_idx >= omitted_center*3 ? deriv_idx+3 : deriv_idx;
          default: assert(deriv_order<=1); // not implemented, won't be needed when Libint computes all derivatives
        }
        assert(false); // unreachable
        return 0;
      }

      /// computes shell set of integrals
      /// \note result is stored in row-major order
      const real_t* compute(const libint2::Shell& s1,
                            const libint2::Shell& s2) {

        // can only handle 1 contraction at a time
        assert(s1.ncontr() == 1 && s2.ncontr() == 1);

        const auto l1 = s1.contr[0].l;
        const auto l2 = s2.contr[0].l;

        // if want nuclear, make sure there is at least one nucleus .. otherwise the user likely forgot to call set_params
        if (type_ == nuclear and nparams() == 0)
          throw std::runtime_error("libint2::OneBodyEngine<nuclear>, but no charges found; forgot to call set_params()?");

        const auto n1 = s1.size();
        const auto n2 = s2.size();
        const auto n12 = n1 * n2;
        const auto ncart1 = s1.cartesian_size();
        const auto ncart2 = s2.cartesian_size();
        const auto ncart12 = ncart1 * ncart2;
        const auto nops = nopers();

        const auto tform_to_solids = s1.contr[0].pure || s2.contr[0].pure;

        // assert # of primitive pairs
        const auto nprim1 = s1.nprim();
        const auto nprim2 = s2.nprim();
        const auto nprimpairs = nprim1 * nprim2;
        assert(nprimpairs <= primdata_.size());

        auto nparam_sets = nparams();

        // how many shell sets will be returned?
        auto num_shellsets = nshellsets();
        // Libint computes derivatives with respect to one center fewer, will use translational invariance to recover
        const auto geometry_independent_operator = type_ == overlap || type_ == kinetic;
        const auto num_deriv_centers_computed = geometry_independent_operator ? 1 : 2;
        auto num_shellsets_computed = nopers() *
                                      num_geometrical_derivatives(num_deriv_centers_computed,
                                                                  deriv_order_);
        // size of ints block computed by Libint
        const auto target_buf_size = num_shellsets_computed * ncart12;

        // will use scratch_ if:
        // - Coulomb ints are computed 1 charge at a time, contributions are accumulated in scratch_ (unless la==lb==0)
        // - derivatives on the missing center need to be reconstructed (no need to accumulate into scratch though)
        // will only use scratch to accumulate ints when
        const auto accumulate_ints_in_scratch = (type_ == nuclear);

        // adjust max angular momentum, if needed
        const auto lmax = std::max(l1, l2);
        assert (lmax <= lmax_);

        // where cartesian ints are located varies, sometimes we compute them in scratch, etc.
        // this is the most likely location
        auto cartesian_ints = primdata_[0].stack;

        // simple (s|s) ints will be computed directly and accumulated in the first element of stack
        const auto compute_directly = lmax == 0 && deriv_order_ == 0 && (type_ == overlap || type_ == nuclear);
        if (compute_directly) {
          primdata_[0].stack[0] = 0;
        }
        else if (accumulate_ints_in_scratch)
          memset(static_cast<void*>(&scratch_[0]), 0, sizeof(real_t)*target_buf_size);

        // loop over accumulation batches
        for(auto pset=0u; pset!=nparam_sets; ++pset) {

          if (type_!=nuclear) assert(nparam_sets == 1);

          auto p12 = 0;
          for(auto p1=0; p1!=nprim1; ++p1) {
            for(auto p2=0; p2!=nprim2; ++p2, ++p12) {
              compute_primdata(primdata_[p12],s1,s2,p1,p2,pset);
            }
          }
          primdata_[0].contrdepth = p12;

          if (compute_directly) {
            auto& result = cartesian_ints[0];
            switch (type_) {
              case overlap:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12]._0_Overlap_0_x[0]
                          * primdata_[p12]._0_Overlap_0_y[0]
                          * primdata_[p12]._0_Overlap_0_z[0];
                  break;
              case nuclear:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12].LIBINT_T_S_ELECPOT_S(0)[0];
                  break;
              default:
                assert(false);
            }
            primdata_[0].targets[0] = cartesian_ints;
          }
          else {

            buildfnptrs_[s1.contr[0].l*hard_lmax_ + s2.contr[0].l](&primdata_[0]);
            cartesian_ints = primdata_[0].targets[0];

            if (accumulate_ints_in_scratch) {
              cartesian_ints = &scratch_[0];
              std::transform(primdata_[0].targets[0], primdata_[0].targets[0] + target_buf_size,
                             &scratch_[0],
                             &scratch_[0], std::plus<real_t>());

              // need to reconstruct derivatives of nuclear ints for each nucleus
              if (deriv_order_ > 0){
                const auto nints_per_center = target_buf_size/2;
                // first two blocks are derivatives with respect to Gaussian positions
                // rest are derivs with respect to nuclear coordinates
                auto dest = &scratch_[0] + (2+pset)*nints_per_center;
                auto src = primdata_[0].targets[0];
                for(auto i=0; i!=nints_per_center; ++i) {
                  dest[i] = -src[i];
                }
                src = primdata_[0].targets[0] + nints_per_center;
                for(auto i=0; i!=nints_per_center; ++i) {
                  dest[i] -= src[i];
                }
                num_shellsets_computed+=3; // we just added 3 shell sets
              } // reconstruct derivatives

            }
          } // ltot != 0

        } // pset (accumulation batches)

        auto result = cartesian_ints; // will be adjusted as we proceed
        if (tform_to_solids) {
          // where do spherical ints go?
          auto* spherical_ints = (cartesian_ints == &scratch_[0]) ? scratch2_ : &scratch_[0];
          result = spherical_ints;

          // transform to solid harmonics, one shell set at a time:
          // for each computed shell set ...
          for(auto s=0ul; s!=num_shellsets_computed; ++s, cartesian_ints+=ncart12) {
            // ... find its index in the target shell set:
            // 1. if regular ints do nothing
            // 2. for derivatives the target set includes derivatives w.r.t omitted centers,
            //    to be computed later (for all ints) or already computed (for nuclear);
            //    in the former case the "omitted" set of derivatives always comes at the end
            //    hence the index of the current shell set does not change (for 2-body ints
            //    the rules are different, but Libint will eliminate the need to reconstruct via
            //    translational invariance soon, so this logic will be unnecessary).
            auto s_target = s;
            // .. and compute the destination
            spherical_ints = result + n12 * s_target;
            if (s1.contr[0].pure && s2.contr[0].pure) {
              libint2::solidharmonics::tform(l1, l2, cartesian_ints, spherical_ints);
            }
            else {
              if (s1.contr[0].pure)
                libint2::solidharmonics::tform_rows(l1, n2, cartesian_ints, spherical_ints);
              else
                libint2::solidharmonics::tform_cols(n1, l2, cartesian_ints, spherical_ints);
            }
          } // loop cartesian shell set

        } // tform to solids

        // if computing derivatives of ints of geometry-independent operators
        // compute the omitted derivatives using translational invariance
        if (deriv_order_ > 0 && geometry_independent_operator) {
          assert(deriv_order_ == 1); // assuming 1st-order derivs here, arbitrary order later

          const auto nints_computed = n12*num_shellsets_computed; // target # of ints is twice this

          // make sure there is enough room left in libint stack
          // if not, copy into scratch2_
          if (not tform_to_solids) {
            const auto stack_size_remaining = stack_size_ - (result-primdata_[0].stack) - nints_computed;
            const auto copy_to_scratch2 = stack_size_remaining < nints_computed;
            if (copy_to_scratch2) {
              // this is tricky ... copy does not allow scratch2_ in [result, result + nints_computed)
              // but this would only happen in scratch2_ == result, but definition of scratch2_ ensures this
              std::copy(result, result + nints_computed, scratch2_);
              result = scratch2_;
            }
          }

          const auto src = result;
          const auto dest = result + nints_computed;
          for(auto f=0ul; f!=nints_computed; ++f) {
            dest[f] = -src[f];
          }
        } // rebuild omitted derivatives of Cartesian ints

        return result;
      }

      inline void compute_primdata(Libint_t& primdata,
                                   const Shell& s1, const Shell& s2,
                                   size_t p1, size_t p2,
                                   size_t oset);

    private:
      operator_type type_;
      std::vector<Libint_t> primdata_;
      size_t stack_size_;  // amount allocated by libint2_init_xxx in primdata_[0].stack
      int lmax_;
      int hard_lmax_;      // max L supported by library for this operator type + 1
      size_t deriv_order_;
      any params_;
      std::shared_ptr<coulomb_core_eval_t> fm_eval_; // this is for Coulomb only
      std::vector<real_t> scratch_; // for transposes and/or transforming to solid harmonics
      real_t* scratch2_;            // &scratch_[0] points to the first block large enough to hold all target ints
                                    // scratch2_ points to second such block. It could point into scratch_ or at primdata_[0].stack
      typedef void (*buildfnptr_t)(const Libint_t*);
      buildfnptr_t* buildfnptrs_;

      void reset_scratch() {
        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;
        const auto target_shellset_size = nshellsets() * ncart_max * ncart_max;
        // need to be able to hold 2 sets of target shellsets: the worst case occurs when dealing with
        // 1-body Coulomb ints derivatives ... have 2+natom derivative sets that are stored in scratch
        // then need to transform to solids. To avoid copying back and forth make sure that there is enough
        // room to transform all ints and save them in correct order in single pass
        const auto need_extra_large_scratch = stack_size_ < target_shellset_size;
        scratch_.resize(need_extra_large_scratch ? 2*target_shellset_size : target_shellset_size);
        scratch2_ = need_extra_large_scratch ? &scratch_[target_shellset_size] : primdata_[0].stack;
      }

      void initialize() {
        assert(libint2::initialized());
        assert(deriv_order_ <= LIBINT2_DERIV_ONEBODY_ORDER);

#define BOOST_PP_ONEBODYENGINE_MCR2(r,product)                                                                \
         if (type_ == BOOST_PP_TUPLE_ELEM(2,0,product) && deriv_order_ == BOOST_PP_TUPLE_ELEM(2,1,product) ) {\
           assert(lmax_ <= BOOST_PP_CAT(LIBINT2_MAX_AM_ ,                                                     \
                                        BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_OPERATOR_LIST,                      \
                                                         BOOST_PP_TUPLE_ELEM(2,0,product) )                   \
                                       ) );                                                                   \
           stack_size_ =                                                                                      \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_need_memory_ ,                                       \
                                      BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_OPERATOR_LIST,                        \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                 )(lmax_);                                                                    \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_init_ ,                                              \
                                      BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_OPERATOR_LIST,                        \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                )(&primdata_[0], lmax_, 0);                                                   \
           buildfnptrs_ = to_ptr1(                                                                            \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_build_ ,                                             \
                                      BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_OPERATOR_LIST,                        \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                )                                                                             \
                                 );                                                                           \
           hard_lmax_ =           BOOST_PP_CAT(                                                               \
                                    LIBINT2_MAX_AM_ ,                                                         \
                                    BOOST_PP_CAT(                                                             \
                                      BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_OPERATOR_LIST,                        \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) ),                    \
                                      BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),     \
                                                    BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()        \
                                                )                                                             \
                                    )                                                                         \
                                  ) + 1;                                                                      \
           reset_scratch();                                                                                   \
           return;                                                                                            \
         }

BOOST_PP_LIST_FOR_EACH_PRODUCT ( BOOST_PP_ONEBODYENGINE_MCR2, 2, (BOOST_PP_ONEBODY_OPERATOR_INDEX_LIST, BOOST_PP_ONEBODY_DERIV_ORDER_LIST) )

        assert(false); // either deriv_order_ or type_ is wrong
      } // initialize()

      void finalize() {
        if (primdata_.size() != 0) {

#define BOOST_PP_ONEBODYENGINE_MCR3(r,product)                                                                \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_cleanup_ ,                                           \
                                      BOOST_PP_LIST_AT(BOOST_PP_ONEBODY_OPERATOR_LIST,                        \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                )(&primdata_[0]);                                                             \
           return;

BOOST_PP_LIST_FOR_EACH_PRODUCT ( BOOST_PP_ONEBODYENGINE_MCR3, 2, (BOOST_PP_ONEBODY_OPERATOR_INDEX_LIST, BOOST_PP_ONEBODY_DERIV_ORDER_LIST) )

        }
      } // finalize()

      //-------
      // utils
      //-------
      inline unsigned int nparams() const;
      inline unsigned int nopers() const;
      /// if Params == operator_traits<type>::oper_params_type, will return any(params)
      /// else will set return any initialized with default value for operator_traits<type>::oper_params_type
      /// @param throw_if_wrong_type if true, and Params != operator_traits<type>::oper_params_type, will throw std::bad_cast
      template <typename Params>
      static any enforce_params_type(operator_type type,
                                     const Params& params,
                                     bool throw_if_wrong_type = not std::is_same<Params,empty_pod>::value);

  }; // struct OneBodyEngine

  template <OneBodyEngine::operator_type Op> struct OneBodyEngine::operator_traits {
      typedef struct {} oper_params_type;
      static oper_params_type default_params() { return oper_params_type{}; }
      static constexpr unsigned int nopers = 1;
  };

  template <> struct OneBodyEngine::operator_traits<OneBodyEngine::nuclear> {
      /// point charges and their positions
      typedef std::vector<std::pair<double, std::array<double, 3>>> oper_params_type;
      static oper_params_type default_params() { return oper_params_type{}; }
      static constexpr unsigned int nopers = 1;
  };

  template <> struct OneBodyEngine::operator_traits<OneBodyEngine::emultipole1> {
      /// Cartesian coordinates of the origin with respect to which the dipole moment is defined
      typedef std::array<double, 3> oper_params_type;
      static oper_params_type default_params() { return oper_params_type{{0.0,0.0,0.0}}; }
      static constexpr unsigned int nopers = 4; //!< overlap + 3 dipole components
  };
  template <> struct OneBodyEngine::operator_traits<OneBodyEngine::emultipole2> {
      /// Cartesian coordinates of the origin with respect to which the multipole moments are defined
      typedef OneBodyEngine::operator_traits<OneBodyEngine::emultipole1>::oper_params_type oper_params_type;
      static oper_params_type default_params() { return OneBodyEngine::operator_traits<OneBodyEngine::emultipole1>::default_params(); }
      static constexpr unsigned int nopers = OneBodyEngine::operator_traits<OneBodyEngine::emultipole1>::nopers + 6; //!< overlap + 3 dipoles + 6 quadrupoles
  };
  template <> struct OneBodyEngine::operator_traits<OneBodyEngine::emultipole3> {
      /// Cartesian coordinates of the origin with respect to which the multipole moments are defined
      typedef OneBodyEngine::operator_traits<OneBodyEngine::emultipole1>::oper_params_type oper_params_type;
      static oper_params_type default_params() { return OneBodyEngine::operator_traits<OneBodyEngine::emultipole1>::default_params(); }
      static constexpr unsigned int nopers = OneBodyEngine::operator_traits<OneBodyEngine::emultipole2>::nopers + 10;
  };

  inline unsigned int OneBodyEngine::nparams() const {
    switch (type_) {
      case nuclear:
        return params_.as<operator_traits<nuclear>::oper_params_type>().size();
      default:
        return 1;
    }
    return 1;
  }
  inline unsigned int OneBodyEngine::nopers() const {
    switch (type_) {
#define BOOST_PP_ONEBODYENGINE_MCR4(r,data,i,elem)  case i : return operator_traits< static_cast<OneBodyEngine::operator_type> ( i ) >::nopers;
BOOST_PP_LIST_FOR_EACH_I ( BOOST_PP_ONEBODYENGINE_MCR4, _, BOOST_PP_ONEBODY_OPERATOR_LIST)
      default: break;
    }
    assert(false); // unreachable
    return 0;
  }
  template <typename Params>
  any OneBodyEngine::enforce_params_type(operator_type type,
                                         const Params& params,
                                         bool throw_if_wrong_type) {
    any result;
    switch(type) {

#define BOOST_PP_ONEBODYENGINE_MCR5(r,data,i,elem)                                                           \
      case i :                                                                                               \
      if (std::is_same<Params,operator_traits< static_cast<operator_type> ( i ) >::oper_params_type>::value) \
        result = params;                                                                                     \
      else {                                                                                                 \
        if (throw_if_wrong_type) throw std::bad_cast();                                                      \
        result = operator_traits<static_cast<operator_type> ( i ) >::default_params();                       \
      }                                                                                                      \
      break;

BOOST_PP_LIST_FOR_EACH_I ( BOOST_PP_ONEBODYENGINE_MCR5, _, BOOST_PP_ONEBODY_OPERATOR_LIST)

      default:
        assert(false); // missed a case?
    }
    return result;
  }

  inline void OneBodyEngine::compute_primdata(Libint_t& primdata,
                                              const Shell& s1, const Shell& s2,
                                              size_t p1, size_t p2,
                                              size_t oset) {

    const auto& A = s1.O;
    const auto& B = s2.O;

    const auto alpha1 = s1.alpha[p1];
    const auto alpha2 = s2.alpha[p2];

    const auto c1 = s1.contr[0].coeff[p1];
    const auto c2 = s2.contr[0].coeff[p2];

    const auto gammap = alpha1 + alpha2;
    const auto oogammap = 1.0 / gammap;
    const auto rhop_over_alpha1 = alpha2 * oogammap;
    const auto rhop = alpha1 * rhop_over_alpha1;
    const auto Px = (alpha1 * A[0] + alpha2 * B[0]) * oogammap;
    const auto Py = (alpha1 * A[1] + alpha2 * B[1]) * oogammap;
    const auto Pz = (alpha1 * A[2] + alpha2 * B[2]) * oogammap;
    const auto AB_x = A[0] - B[0];
    const auto AB_y = A[1] - B[1];
    const auto AB_z = A[2] - B[2];
    const auto AB2_x = AB_x * AB_x;
    const auto AB2_y = AB_y * AB_y;
    const auto AB2_z = AB_z * AB_z;
    const auto AB2 = AB2_x + AB2_y + AB2_z;

    assert (LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD);

    // overlap and kinetic energy ints don't use HRR, hence VRR on both centers
    // Coulomb potential do HRR on center 1 only
#if LIBINT2_DEFINED(eri,PA_x)
      primdata.PA_x[0] = Px - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
      primdata.PA_y[0] = Py - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
      primdata.PA_z[0] = Pz - A[2];
#endif

    if (type_ != nuclear) {

#if LIBINT2_DEFINED(eri,PB_x)
      primdata.PB_x[0] = Px - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
      primdata.PB_y[0] = Py - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
      primdata.PB_z[0] = Pz - B[2];
#endif
    }

    if (type_ == emultipole1 || type_ == emultipole2 || type_ == emultipole3) {
      auto& O = params_.as<operator_traits<emultipole1>::oper_params_type>(); // same as emultipoleX
#if LIBINT2_DEFINED(eri,BO_x)
      primdata.BO_x[0] = B[0] - O[0];
#endif
#if LIBINT2_DEFINED(eri,BO_y)
      primdata.BO_y[0] = B[1] - O[1];
#endif
#if LIBINT2_DEFINED(eri,BO_z)
      primdata.BO_z[0] = B[2] - O[2];
#endif
    }

#if LIBINT2_DEFINED(eri,oo2z)
    primdata.oo2z[0] = 0.5*oogammap;
#endif

    if (type_ == nuclear) { // additional factor for electrostatic potential
      auto& params = params_.as<operator_traits<nuclear>::oper_params_type>();
      const auto& C = params[oset].second;
#if LIBINT2_DEFINED(eri,PC_x)
      primdata.PC_x[0] = Px - C[0];
#endif
#if LIBINT2_DEFINED(eri,PC_y)
      primdata.PC_y[0] = Py - C[1];
#endif
#if LIBINT2_DEFINED(eri,PC_z)
      primdata.PC_z[0] = Pz - C[2];
#endif
      // elecpot uses HRR
#if LIBINT2_DEFINED(eri,AB_x)
      primdata.AB_x[0] = A[0] - B[0];
#endif
#if LIBINT2_DEFINED(eri,AB_y)
      primdata.AB_y[0] = A[1] - B[1];
#endif
#if LIBINT2_DEFINED(eri,AB_z)
      primdata.AB_z[0] = A[2] - B[2];
#endif

    }

    decltype(c1) sqrt_PI(1.77245385090551602729816748334);
    const auto xyz_pfac = sqrt_PI * sqrt(oogammap);
    const auto ovlp_ss_x = exp(- rhop * AB2_x) * xyz_pfac * c1 * c2;
    const auto ovlp_ss_y = exp(- rhop * AB2_y) * xyz_pfac;
    const auto ovlp_ss_z = exp(- rhop * AB2_z) * xyz_pfac;

    primdata._0_Overlap_0_x[0] = ovlp_ss_x;
    primdata._0_Overlap_0_y[0] = ovlp_ss_y;
    primdata._0_Overlap_0_z[0] = ovlp_ss_z;

    if (type_ == kinetic || (deriv_order_ > 0)) {
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
      primdata.two_alpha0_bra[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
      primdata.two_alpha0_ket[0] = 2.0 * alpha2;
#endif
    }

    if (type_ == nuclear) {
#if LIBINT2_DEFINED(eri,rho12_over_alpha1) || LIBINT2_DEFINED(eri,rho12_over_alpha2)
      if (deriv_order_ > 0) {
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
        primdata.rho12_over_alpha1[0] = rhop_over_alpha1;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
        primdata.rho12_over_alpha2[0] = alpha1 * oogammap;
#endif
      }
#endif
#if LIBINT2_DEFINED(eri,PC_x) && LIBINT2_DEFINED(eri,PC_y) && LIBINT2_DEFINED(eri,PC_z)
      const auto PC2 = primdata.PC_x[0] * primdata.PC_x[0] +
                       primdata.PC_y[0] * primdata.PC_y[0] +
                       primdata.PC_z[0] * primdata.PC_z[0];
      const auto U = gammap * PC2;
      const auto mmax = s1.contr[0].l + s2.contr[0].l + deriv_order_;
      auto* fm_ptr = &(primdata.LIBINT_T_S_ELECPOT_S(0)[0]);
      fm_eval_->eval(fm_ptr, U, mmax);

      decltype(U) two_o_sqrt_PI(1.12837916709551257389615890312);
      const auto q = params_.as<operator_traits<nuclear>::oper_params_type>()[oset].first;
      const auto pfac = - q * sqrt(gammap) * two_o_sqrt_PI * ovlp_ss_x * ovlp_ss_y * ovlp_ss_z;
      const auto m_fence = mmax + 1;
      for(auto m=0; m!=m_fence; ++m) {
        fm_ptr[m] *= pfac;
      }
#endif
    }

  } // OneBodyEngine::compute_primdata()

#undef BOOST_PP_ONEBODY_OPERATOR_LIST
#undef BOOST_PP_ONEBODY_OPERATOR_INDEX_TUPLE
#undef BOOST_PP_ONEBODY_OPERATOR_INDEX_LIST
#undef BOOST_PP_ONEBODY_DERIV_ORDER_TUPLE
#undef BOOST_PP_ONEBODY_DERIV_ORDER_LIST
#undef BOOST_PP_ONEBODYENGINE_MCR2
#undef BOOST_PP_ONEBODYENGINE_MCR3
#undef BOOST_PP_ONEBODYENGINE_MCR4
#undef BOOST_PP_ONEBODYENGINE_MCR5

#endif // LIBINT2_SUPPORT_ONEBODY

  /// types of multiplicative spherically-symmetric two-body kernels known by TwoBodyEngine
  enum MultiplicativeSphericalTwoBodyKernel {
    Coulomb,            //!< \f$ 1/r_{12} \f$
    cGTG,               //!< contracted Gaussian geminal = \f$ \sum_i c_i \exp(- \alpha r_{12}^2) \f$
    cGTG_times_Coulomb, //!< contracted Gaussian geminal times Coulomb
    DelcGTG_square,     //!< (\f$ \nabla \f$ cGTG) \f$ \cdot \f$ (\f$ \nabla \f$ cGTG)
    Delta               //!< \f$ \delta\left({\bf r}_1 - {\bf r}_2\right) \f$
  };

  /// contracted Gaussian geminal = \f$ \sum_i c_i \exp(- \alpha r_{12}^2) \f$, represented as a vector of
  /// {\f$ \alpha_i \f$, \f$ c_i \f$ } pairs
  typedef std::vector<std::pair<double,double>> ContractedGaussianGeminal;

  namespace detail {
    template <int K> struct R12_K_G12_to_Kernel;
    template <> struct R12_K_G12_to_Kernel<-1> {
        static const MultiplicativeSphericalTwoBodyKernel value = cGTG_times_Coulomb;
    };
    template <> struct R12_K_G12_to_Kernel<0> {
        static const MultiplicativeSphericalTwoBodyKernel value = cGTG;
    };
    template <> struct R12_K_G12_to_Kernel<2> {
        static const MultiplicativeSphericalTwoBodyKernel value = DelcGTG_square;
    };

    template <MultiplicativeSphericalTwoBodyKernel Kernel> struct TwoBodyEngineDispatcher;

  } // namespace detail

  template <MultiplicativeSphericalTwoBodyKernel Kernel> struct TwoBodyEngineTraits;
  template <> struct TwoBodyEngineTraits<Coulomb> {
      typedef libint2::FmEval_Chebyshev3<double> core_eval_type;
      //typedef libint2::FmEval_Chebyshev7<double> core_eval_type;
      //typedef libint2::FmEval_Taylor<double, 7> core_eval_type;
      typedef struct {} oper_params_type;
  };
  template <> struct TwoBodyEngineTraits<cGTG> {
      typedef libint2::GaussianGmEval<real_t, 0> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
  };
  template <> struct TwoBodyEngineTraits<cGTG_times_Coulomb> {
      typedef libint2::GaussianGmEval<real_t, -1> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
  };
  template <> struct TwoBodyEngineTraits<DelcGTG_square> {
      typedef libint2::GaussianGmEval<real_t, 2> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
  };

  struct delta_gm_eval {
      void operator()(double* Gm,
                      double rho,
                      double T,
                      int mmax) {
        const static auto one_over_two_pi = 1.0 / (2.0 * M_PI);
        const auto G0 = exp(-T) * rho * one_over_two_pi;
        std::fill(Gm, Gm+mmax+1, G0);
      }
  };
  template <> struct TwoBodyEngineTraits<Delta> {
      typedef libint2::GenericGmEval<delta_gm_eval> core_eval_type; // core ints are too trivial to bother
      typedef struct {} oper_params_type;
  };

#ifdef LIBINT2_SUPPORT_ERI
  /**
   * TwoBodyEngine computes (ab|O|cd) (i.e. <em>four-center</em>) integrals over
   * a two-body kernel of type MultiplicativeSphericalTwoBodyKernel using Obara-Saika-Ahlrichs relations.
   *
   * \tparam Kernel kernel type, the supported values are enumerated by MultiplicativeSphericalTwoBodyKernel
   */
  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  class TwoBodyEngine {
    public:
      typedef typename libint2::TwoBodyEngineTraits<Kernel>::oper_params_type oper_params_type;

      /// creates a default (unusable) TwoBodyEngine
      TwoBodyEngine() : primdata_(), lmax_(-1), core_eval_pack_() {
        set_precision(std::numeric_limits<real_t>::epsilon());
      }

      /// Constructs a (usable) TwoBodyEngine

      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_order if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_order
      /// \param precision specifies the target precision with which the integrals will be computed; the default is the "epsilon"
      ///        of \c real_t type, given by \c std::numeric_limits<real_t>::epsilon(). The precision
      ///        control is somewhat empirical, hence be conservative. \sa set_precision()
      /// \param oper_params specifies the operator parameters. The type of \c oper_params depends on \c Kernel as follows:
      ///        <ol>
      ///        <li> \c Coulomb : empty type (does not need to be provided) </li>
      ///        <li> \c cGTG : ContractedGaussianGeminal </li>
      ///        <li> \c cGTG_times_Coulomb : ContractedGaussianGeminal </li>
      ///        <li> \c DelcGTG_square : ContractedGaussianGeminal </li>
      ///        </ol>
      /// \warning currently derivative integrals are not supported
      /// \warning currently only one-contraction Shell objects are supported; i.e. generally-contracted Shells are not yet supported
      TwoBodyEngine(size_t max_nprim, int max_l,
                    int deriv_order = 0,
                    real_t precision = std::numeric_limits<real_t>::epsilon(),
                    const oper_params_type& oparams = oper_params_type()) :
        primdata_(max_nprim * max_nprim * max_nprim * max_nprim),
        spbra_(max_nprim), spket_(max_nprim),
        lmax_(max_l), deriv_order_(deriv_order),
        core_eval_pack_(core_eval_type::instance(4*lmax_ + deriv_order,
                                                 std::numeric_limits<real_t>::epsilon()),
                        core_eval_scratch_type(4*lmax_ + deriv_order)
                       )
      {
        set_precision(precision);
        initialize();
        init_core_ints_params(oparams);
      }

      /// move constructor
      // intel does not support "move ctor = default"
      TwoBodyEngine(TwoBodyEngine&& other) :
        primdata_(std::move(other.primdata_)),
        spbra_(std::move(other.spbra_)), spket_(std::move(other.spket_)),
        lmax_(other.lmax_), deriv_order_(other.deriv_order_),
        precision_(other.precision_), ln_precision_(other.ln_precision_),
        core_eval_pack_(std::move(other.core_eval_pack_)),
        core_ints_params_(std::move(other.core_ints_params_)),
        scratch_(std::move(other.scratch_)) {
      }

      /// (deep) copy constructor
      TwoBodyEngine(const TwoBodyEngine& other) :
        primdata_(other.primdata_),
        spbra_(other.spbra_), spket_(other.spket_),
        lmax_(other.lmax_), deriv_order_(other.deriv_order_),
        precision_(other.precision_), ln_precision_(other.ln_precision_),
        core_eval_pack_(other.core_eval_pack_), core_ints_params_(other.core_ints_params_)
      {
        initialize();
      }

      ~TwoBodyEngine() {
        finalize();
      }

      /// move assignment
      // intel does not support "move asgnmt = default"
      TwoBodyEngine& operator=(TwoBodyEngine&& other) {
        primdata_ = std::move(other.primdata_);
        lmax_ = other.lmax_;
        deriv_order_ = other.deriv_order_;
        precision_ = other.precision_;
        ln_precision_ = other.ln_precision_;
        core_eval_pack_ = std::move(other.core_eval_pack_);
        core_ints_params_ = std::move(other.core_ints_params_);
        scratch_ = std::move(other.scratch_);
        return *this;
      }

      /// (deep) copy assignment
      TwoBodyEngine& operator=(const TwoBodyEngine& other)
      {
        primdata_ = other.primdata_;
        lmax_ = other.lmax_;
        deriv_order_ = other.deriv_order_;
        precision_ = other.precision_;
        ln_precision_ = other.ln_precision_;
        core_eval_pack_ = other.core_eval_pack_;
        core_ints_params_ = other.core_ints_params_;
        initialize();
        return *this;
      }

      static bool skip_core_ints;

#ifdef LIBINT2_ENGINE_TIMERS
      Timers<3> timers; // timers[0] -> prereqs
                        // timers[1] -> build (only meaningful if LIBINT2_PROFILE is not defined
                        // timers[2] -> tform
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
      struct class_id {
          size_t l[4];
          template <typename I> class_id(I l0, I l1, I l2, I l3) {
            l[0] = l0;
            l[1] = l1;
            l[2] = l2;
            l[3] = l3;
          }
          bool operator<(const class_id& other) const {
            return ordinal(l) < ordinal(other.l);
          }
          static size_t ordinal(const size_t (&l)[4]) {
            return ((l[0]*LIBINT2_MAX_AM_eri + l[1])*LIBINT2_MAX_AM_eri + l[2])*LIBINT2_MAX_AM_eri + l[3];
          }
          std::string to_string() const {
            std::ostringstream oss;
            oss << "(" << Shell::am_symbol(l[0]) << Shell::am_symbol(l[1])
                << "|" << Shell::am_symbol(l[2]) << Shell::am_symbol(l[3]) << ")";
            return oss.str();
          }
      };
      struct class_profile {
          double prereqs;
          double build_hrr;
          double build_vrr;
          double tform;
          size_t nshellset; // total number of shell sets
          size_t nprimset;  // total number of primitive sets
          class_profile() { clear(); }
          class_profile(const class_profile& other) = default;
          void clear() {
            prereqs = build_hrr = build_vrr = tform = 0.;
            nprimset = nshellset = 0;
          }
      };
      std::map<class_id, class_profile> class_profiles;
#  endif
#endif

      /// computes shell set of integrals
      /// \note result is stored in the "chemists" form, i.e. (tbra1 tbra2 |tket1 tket2), in row-major order
      const real_t* compute(const libint2::Shell& tbra1,
                            const libint2::Shell& tbra2,
                            const libint2::Shell& tket1,
                            const libint2::Shell& tket2) {

        //
        // i.e. bra and ket refer to chemists bra and ket
        //

        // can only handle 1 contraction at a time
        assert(tbra1.ncontr() == 1 && tbra2.ncontr() == 1 &&
               tket1.ncontr() == 1 && tket2.ncontr() == 1);

        // derivatives not supported for now
        assert(deriv_order_ == 0);

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // standard angular momentum ordering
        auto swap_bra = (tbra1.contr[0].l < tbra2.contr[0].l);
        auto swap_ket = (tket1.contr[0].l < tket2.contr[0].l);
        auto swap_braket = (tbra1.contr[0].l + tbra2.contr[0].l > tket1.contr[0].l + tket2.contr[0].l);
#else // orca angular momentum ordering
        auto swap_bra = (tbra1.contr[0].l > tbra2.contr[0].l);
        auto swap_ket = (tket1.contr[0].l > tket2.contr[0].l);
        auto swap_braket = (tbra1.contr[0].l + tbra2.contr[0].l < tket1.contr[0].l + tket2.contr[0].l);
#endif
        const auto& bra1 = swap_braket ? (swap_ket ? tket2 : tket1) : (swap_bra ? tbra2 : tbra1);
        const auto& bra2 = swap_braket ? (swap_ket ? tket1 : tket2) : (swap_bra ? tbra1 : tbra2);
        const auto& ket1 = swap_braket ? (swap_bra ? tbra2 : tbra1) : (swap_ket ? tket2 : tket1);
        const auto& ket2 = swap_braket ? (swap_bra ? tbra1 : tbra2) : (swap_ket ? tket1 : tket2);

        const bool tform = tbra1.contr[0].pure || tbra2.contr[0].pure || tket1.contr[0].pure || tket2.contr[0].pure;
        const bool use_scratch = (swap_braket || swap_bra || swap_ket || tform);

        // assert # of primitive pairs
        auto nprim_bra1 = bra1.nprim();
        auto nprim_bra2 = bra2.nprim();
        auto nprim_ket1 = ket1.nprim();
        auto nprim_ket2 = ket2.nprim();

        // adjust max angular momentum, if needed
        auto lmax = std::max(std::max(bra1.contr[0].l, bra2.contr[0].l), std::max(ket1.contr[0].l, ket2.contr[0].l));
        assert (lmax <= lmax_);
        if (lmax == 0) // (ss|ss) ints will be accumulated in the first element of stack
          primdata_[0].stack[0] = 0;

#ifdef LIBINT2_ENGINE_PROFILE_CLASS
        class_id id(bra1.contr[0].l, bra2.contr[0].l, ket1.contr[0].l, ket2.contr[0].l);
        if (class_profiles.find(id) == class_profiles.end()) {
          class_profile dummy;
          class_profiles[id] = dummy;
        }
#endif

        // compute primitive data
#ifdef LIBINT2_ENGINE_TIMERS
        timers.start(0);
#endif
        {
          auto p = 0;
          // this is far less aggressive than should be, but proper analysis
          // involves both bra and ket *bases* and thus cannot be done on shell-set basis
          // probably ln_precision_/2 - 10 is enough
          spbra_.init(bra1, bra2, ln_precision_);
          spket_.init(ket1, ket2, ln_precision_);
          const auto npbra = spbra_.primpairs.size();
          const auto npket = spket_.primpairs.size();
          for(auto pb=0; pb!=npbra; ++pb) {
            for(auto pk=0; pk!=npket; ++pk) {
              if (spbra_.primpairs[pb].scr + spket_.primpairs[pk].scr > ln_precision_) {
                if (compute_primdata(primdata_[p],bra1,bra2,ket1,ket2,spbra_,pb,spket_,pk)) {
                  ++p;
                }
              }
            }
          }
          primdata_[0].contrdepth = p;
        }

#ifdef LIBINT2_ENGINE_TIMERS
        const auto t0 = timers.stop(0);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
        class_profiles[id].prereqs += t0.count();
        if (primdata_[0].contrdepth != 0) {
          class_profiles[id].nshellset += 1;
          class_profiles[id].nprimset += primdata_[0].contrdepth;
        }
#  endif
#endif

        // all primitive combinations screened out? return zeroes
        if (primdata_[0].contrdepth == 0) {
          const size_t n = bra1.size() * bra2.size() * ket1.size() * ket2.size();
          memset(primdata_[0].stack, 0, sizeof(real_t)*n);
          return primdata_[0].stack;
        }

        real_t* result = nullptr;

        if (lmax == 0) { // (ss|ss)
#ifdef LIBINT2_ENGINE_TIMERS
          timers.start(1);
#endif
          auto& stack = primdata_[0].stack[0];
          stack = 0;
          for(auto p=0; p != primdata_[0].contrdepth; ++p)
            stack += primdata_[p].LIBINT_T_SS_EREP_SS(0)[0];
          primdata_[0].targets[0] = primdata_[0].stack;
#ifdef LIBINT2_ENGINE_TIMERS
          const auto t1 = timers.stop(1);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
          class_profiles[id].build_vrr += t1.count();
#  endif
#endif

          result = primdata_[0].targets[0];
        }
        else { // not (ss|ss)
#ifdef LIBINT2_ENGINE_TIMERS
#    ifdef LIBINT2_PROFILE
          const auto t1_hrr_start = primdata_[0].timers->read(0);
          const auto t1_vrr_start = primdata_[0].timers->read(1);
#    endif
          timers.start(1);
#endif
          LIBINT2_PREFIXED_NAME(libint2_build_eri)[bra1.contr[0].l][bra2.contr[0].l][ket1.contr[0].l][ket2.contr[0].l](&primdata_[0]);
#ifdef LIBINT2_ENGINE_TIMERS
          const auto t1 = timers.stop(1);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
#    ifndef LIBINT2_PROFILE
          class_profiles[id].build_vrr += t1.count();
#    else
          class_profiles[id].build_hrr += primdata_[0].timers->read(0) - t1_hrr_start;
          class_profiles[id].build_vrr += primdata_[0].timers->read(1) - t1_vrr_start;
#    endif
#  endif
#endif

          result = primdata_[0].targets[0];

#ifdef LIBINT2_ENGINE_TIMERS
          timers.start(2);
#endif

          // if needed, permute and transform
          if (use_scratch) {

            constexpr auto using_scalar_real = std::is_same<double,real_t>::value || std::is_same<float,real_t>::value;
            static_assert(using_scalar_real, "Libint2 C++11 API only supports fundamental real types");
            typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;

            // a 2-d view of the 4-d source tensor
            const auto nr1_cart = bra1.cartesian_size();
            const auto nr2_cart = bra2.cartesian_size();
            const auto nc1_cart = ket1.cartesian_size();
            const auto nc2_cart = ket2.cartesian_size();
            const auto ncol_cart = nc1_cart * nc2_cart;
            const auto nr1 = bra1.size();
            const auto nr2 = bra2.size();
            const auto nc1 = ket1.size();
            const auto nc2 = ket2.size();
            const auto nrow = nr1 * nr2;
            const auto ncol = nc1 * nc2;

            // a 2-d view of the 4-d target tensor
            const auto nr1_tgt = tbra1.size();
            const auto nr2_tgt = tbra2.size();
            const auto nc1_tgt = tket1.size();
            const auto nc2_tgt = tket2.size();
            const auto ncol_tgt = nc1_tgt * nc2_tgt;

            // transform to solid harmonics first, then unpermute, if necessary
            auto mainbuf = result;
            auto scratchbuf = &scratch_[0];
            if (bra1.contr[0].pure) {
              libint2::solidharmonics::transform_first(bra1.contr[0].l, nr2_cart*ncol_cart,
                                                       mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }
            if (bra2.contr[0].pure) {
              libint2::solidharmonics::transform_inner(bra1.size(), bra2.contr[0].l, ncol_cart,
                                                       mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }
            if (ket1.contr[0].pure) {
              libint2::solidharmonics::transform_inner(nrow, ket1.contr[0].l, nc2_cart,
                                                       mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }
            if (ket2.contr[0].pure) {
              libint2::solidharmonics::transform_last(bra1.size()*bra2.size()*ket1.size(), ket2.contr[0].l,
                                                      mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }

            // loop over rows of the source matrix
            const auto* src_row_ptr = mainbuf;
            auto tgt_ptr = scratchbuf;
            for(auto r1=0; r1!=nr1; ++r1) {
              for(auto r2=0; r2!=nr2; ++r2, src_row_ptr+=ncol) {

                typedef Eigen::Map<const Matrix> ConstMap;
                typedef Eigen::Map<Matrix> Map;
                typedef Eigen::Map<Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > StridedMap;

                // represent this source row as a matrix
                ConstMap src_blk_mat(src_row_ptr, nc1, nc2);

                // and copy to the block of the target matrix
                if (swap_braket) {
                  // if swapped bra and ket, a row of source becomes a column of target
                  // source row {r1,r2} is mapped to target column {r1,r2} if !swap_ket, else to {r2,r1}
                  const auto tgt_col_idx = !swap_ket ? r1 * nr2 + r2 : r2 * nr1 + r1;
                  StridedMap tgt_blk_mat(tgt_ptr + tgt_col_idx, nr1_tgt, nr2_tgt,
                                         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(nr2_tgt*ncol_tgt,ncol_tgt));
                  if (swap_bra)
                    tgt_blk_mat = src_blk_mat.transpose();
                  else
                    tgt_blk_mat = src_blk_mat;
                }
                else {
                  // source row {r1,r2} is mapped to target row {r1,r2} if !swap_bra, else to {r2,r1}
                  const auto tgt_row_idx = !swap_bra ? r1 * nr2 + r2 : r2 * nr1 + r1;
                  Map tgt_blk_mat(tgt_ptr + tgt_row_idx*ncol, nc1_tgt, nc2_tgt);
                  if (swap_ket)
                    tgt_blk_mat = src_blk_mat.transpose();
                  else
                    tgt_blk_mat = src_blk_mat;
                }

              } // end of loop
            }   // over rows of source

            result = scratchbuf;

          } // if need_scratch => needed to transpose

#ifdef LIBINT2_ENGINE_TIMERS
          const auto t2 = timers.stop(2);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
          class_profiles[id].tform += t2.count();
#  endif
#endif

        } // not (ss|ss)

        return result;
      }

      /** this specifies target precision for computing the integrals.
       *  target precision \f$ \epsilon \f$ is used in 3 ways:
       *  (1) to screen out primitive pairs in ShellPair object for which
       *      \f$ {\rm scr}_{12} = \max|c_1| \max|c_2| \exp(-\rho_{12} |AB|^2)/\gamma_{12} < \epsilon \f$ ;
       *  (2) to screen out primitive quartets outside compute_primdata() for which \f$ {\rm scr}_{12} {\rm scr}_{34} <  \epsilon \f$;
       *  (3) to screen out primitive quartets inside compute_primdata() for which the prefactor of \f$ F_m(\rho, T) \f$ is smaller
       *      than \f$ \epsilon \f$ .
       */
      void set_precision(real_t prec) {
        if (prec <= 0.) {
          precision_ = 0.;
          ln_precision_ = std::numeric_limits<real_t>::lowest();
        }
        else {
          precision_ = prec;
          ln_precision_ = std::log(precision_);
        }
      }
      /// @return the target precision for computing the integrals
      /// @sa set_precision(real_t)
      real_t precision() const {
        return precision_;
      }

      void print_timers() {
#ifdef LIBINT2_ENGINE_TIMERS
        std::cout << "timers: prereq = " << timers.read(0);
#  ifndef LIBINT2_PROFILE // if libint's profiling was on, engine's build timer will include its overhead
                          // do not report it, detailed profiling from libint will be printed below
        std::cout << " build = " << timers.read(1);
#  endif
        std::cout << " tform = " << timers.read(2) << std::endl;
#endif
#ifdef LIBINT2_PROFILE
        std::cout << "build timers: hrr = " << primdata_[0].timers->read(0)
                << " vrr = " << primdata_[0].timers->read(1) << std::endl;
#endif
#ifdef LIBINT2_ENGINE_TIMERS
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
        for(const auto& p: class_profiles) {
          printf("{\"%s\", %10.5lf, %10.5lf, %10.5lf, %10.5lf, %ld, %ld},\n",
                 p.first.to_string().c_str(),
                 p.second.prereqs,
                 p.second.build_vrr,
                 p.second.build_hrr,
                 p.second.tform,
                 p.second.nshellset,
                 p.second.nprimset);
        }
#  endif
#endif
      }

    private:

      inline bool compute_primdata(Libint_t& primdata,
                                   const Shell& sbra1,
                                   const Shell& sbra2,
                                   const Shell& sket1,
                                   const Shell& sket2,
                                   const ShellPair& spbra, size_t pbra,
                                   const ShellPair& spket, size_t pket);

      std::vector<Libint_t> primdata_;
      ShellPair spbra_, spket_;
      int lmax_;
      size_t deriv_order_;
      real_t precision_;
      real_t ln_precision_;

      typedef typename libint2::TwoBodyEngineTraits<Kernel>::core_eval_type core_eval_type;
      typedef typename libint2::detail::CoreEvalScratch<core_eval_type> core_eval_scratch_type;
      typedef libint2::detail::compressed_pair<std::shared_ptr<core_eval_type>,
                                               core_eval_scratch_type>
              core_eval_pack_type;
      core_eval_pack_type core_eval_pack_;

      typedef oper_params_type core_ints_params_type; // currently core ints params are always same type as operator params
      core_ints_params_type core_ints_params_;
      /// converts operator parameters to core ints params
      void init_core_ints_params(const oper_params_type& oparams);

      std::vector<real_t> scratch_; // for transposes and/or transforming to solid harmonics

      friend struct detail::TwoBodyEngineDispatcher<Kernel>;

      void initialize() {
        assert(libint2::initialized());
        assert(lmax_ <= LIBINT2_MAX_AM_eri);
        assert(deriv_order_ <= LIBINT2_DERIV_ERI_ORDER);
	
        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;
        const auto max_shellpair_size = ncart_max * ncart_max;
        const auto max_shellset_size = max_shellpair_size * max_shellpair_size;

        switch(deriv_order_) {

            case 0:
              libint2_init_eri(&primdata_[0], lmax_, 0);
              scratch_.resize(max_shellset_size);
              break;
            case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0
              libint2_init_eri1(&primdata_[0], lmax_, 0);
              scratch_.resize(9 * max_shellset_size);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1
              libint2_init_eri2(&primdata_[0], lmax_, 0);
              scratch_.resize(45 * max_shellset_size);
#endif
              break;
            default: assert(deriv_order_ < 3);
        }

#ifdef LIBINT2_ENGINE_TIMERS
        timers.set_now_overhead(25);
#endif
#ifdef LIBINT2_PROFILE
        primdata_[0].timers->set_now_overhead(25);
#endif
      }

      void finalize() {
        if (primdata_.size() != 0) {
          switch(deriv_order_) {

            case 0:
              libint2_cleanup_eri(&primdata_[0]);
              break;
            case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0
              libint2_cleanup_eri1(&primdata_[0]);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1
              libint2_cleanup_eri2(&primdata_[0]);
#endif
              break;
          }
        }
      }

  }; // struct TwoBodyEngine

  namespace detail {
    template <> struct TwoBodyEngineDispatcher<Coulomb> {
        static void core_eval(TwoBodyEngine<Coulomb>* engine,
                              real_t* Fm,
                              int mmax,
                              real_t T,
                              real_t) {
          engine->core_eval_pack_.first()->eval(Fm, T, mmax);
        }
    };

    template <>
    struct TwoBodyEngineDispatcher<cGTG_times_Coulomb> {
        static void core_eval(TwoBodyEngine<cGTG_times_Coulomb>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_pack_.first()->eval(Gm, rho, T, mmax, engine->core_ints_params_, &engine->core_eval_pack_.second());
        }
    };
    template <>
    struct TwoBodyEngineDispatcher<cGTG> {
        static void core_eval(TwoBodyEngine<cGTG>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_pack_.first()->eval(Gm, rho, T, mmax, engine->core_ints_params_);
        }
    };
    template <>
    struct TwoBodyEngineDispatcher<DelcGTG_square> {
        static void core_eval(TwoBodyEngine<DelcGTG_square>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_pack_.first()->eval(Gm, rho, T, mmax, engine->core_ints_params_);
        }
    };
    template <>
    struct TwoBodyEngineDispatcher<Delta> {
        static void core_eval(TwoBodyEngine<Delta>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_pack_.first()->eval(Gm, rho, T, mmax);
        }
    };
  }

  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  bool TwoBodyEngine<Kernel>::skip_core_ints = false;

  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  inline bool TwoBodyEngine<Kernel>::compute_primdata(Libint_t& primdata,
                                                      const Shell& sbra1,
                                                      const Shell& sbra2,
                                                      const Shell& sket1,
                                                      const Shell& sket2,
                                                      const ShellPair& spbra, size_t pbra,
                                                      const ShellPair& spket, size_t pket) {

    const auto& A = sbra1.O;
    const auto& B = sbra2.O;
    const auto& C = sket1.O;
    const auto& D = sket2.O;
    const auto& AB = spbra.AB;
    const auto& CD = spket.AB;

    const auto& spbrapp = spbra.primpairs[pbra];
    const auto& spketpp = spket.primpairs[pket];
    const auto& pbra1 = spbrapp.p1;
    const auto& pbra2 = spbrapp.p2;
    const auto& pket1 = spketpp.p1;
    const auto& pket2 = spketpp.p2;

    const auto alpha0 = sbra1.alpha[pbra1];
    const auto alpha1 = sbra2.alpha[pbra2];
    const auto alpha2 = sket1.alpha[pket1];
    const auto alpha3 = sket2.alpha[pket2];

    const auto c0 = sbra1.contr[0].coeff[pbra1];
    const auto c1 = sbra2.contr[0].coeff[pbra2];
    const auto c2 = sket1.contr[0].coeff[pket1];
    const auto c3 = sket2.contr[0].coeff[pket2];

    const auto amtot = sbra1.contr[0].l + sket1.contr[0].l +
                       sbra2.contr[0].l + sket2.contr[0].l;

    const auto gammap = alpha0 + alpha1;
    const auto oogammap = spbrapp.one_over_gamma;
    const auto rhop = alpha0 * alpha1 * oogammap;

    const auto gammaq = alpha2 + alpha3;
    const auto oogammaq = spketpp.one_over_gamma;
    const auto rhoq = alpha2 * alpha3 * oogammaq;

    const auto& P = spbrapp.P;
    const auto& Q = spketpp.P;
    const auto PQx = P[0] - Q[0];
    const auto PQy = P[1] - Q[1];
    const auto PQz = P[2] - Q[2];
    const auto PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

    const auto K12 = spbrapp.K * spketpp.K;
    decltype(K12) two_times_M_PI_to_25(34.986836655249725693); // (2 \pi)^{5/2}
    const auto gammapq = gammap + gammaq;
    const auto sqrt_gammapq = sqrt(gammapq);
    const auto oogammapq = 1.0 / (gammapq);
    auto pfac = two_times_M_PI_to_25 * K12 * sqrt_gammapq * oogammapq;
    pfac *= c0 * c1 * c2 * c3;

    if (std::abs(pfac) < precision_)
      return false;

    const auto rho = gammap * gammaq * oogammapq;
    const auto T = PQ2*rho;
    auto* fm_ptr = &(primdata.LIBINT_T_SS_EREP_SS(0)[0]);
    const auto mmax = amtot + deriv_order_;

    if (!skip_core_ints)
      detail::TwoBodyEngineDispatcher<Kernel>::core_eval(this, fm_ptr, mmax, T, rho);

    for(auto m=0; m!=mmax+1; ++m) {
      fm_ptr[m] *= pfac;
    }

    if (mmax == 0)
      return true;

#if LIBINT2_DEFINED(eri,PA_x)
    primdata.PA_x[0] = P[0] - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
    primdata.PA_y[0] = P[1] - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
    primdata.PA_z[0] = P[2] - A[2];
#endif
#if LIBINT2_DEFINED(eri,PB_x)
    primdata.PB_x[0] = P[0] - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
    primdata.PB_y[0] = P[1] - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
    primdata.PB_z[0] = P[2] - B[2];
#endif

#if LIBINT2_DEFINED(eri,QC_x)
    primdata.QC_x[0] = Q[0] - C[0];
#endif
#if LIBINT2_DEFINED(eri,QC_y)
    primdata.QC_y[0] = Q[1] - C[1];
#endif
#if LIBINT2_DEFINED(eri,QC_z)
    primdata.QC_z[0] = Q[2] - C[2];
#endif
#if LIBINT2_DEFINED(eri,QD_x)
    primdata.QD_x[0] = Q[0] - D[0];
#endif
#if LIBINT2_DEFINED(eri,QD_y)
    primdata.QD_y[0] = Q[1] - D[1];
#endif
#if LIBINT2_DEFINED(eri,QD_z)
    primdata.QD_z[0] = Q[2] - D[2];
#endif

#if LIBINT2_DEFINED(eri,AB_x)
    primdata.AB_x[0] = AB[0];
#endif
#if LIBINT2_DEFINED(eri,AB_y)
    primdata.AB_y[0] = AB[1];
#endif
#if LIBINT2_DEFINED(eri,AB_z)
    primdata.AB_z[0] = AB[2];
#endif
#if LIBINT2_DEFINED(eri,BA_x)
    primdata.BA_x[0] = -AB[0];
#endif
#if LIBINT2_DEFINED(eri,BA_y)
    primdata.BA_y[0] = -AB[1];
#endif
#if LIBINT2_DEFINED(eri,BA_z)
    primdata.BA_z[0] = -AB[2];
#endif

#if LIBINT2_DEFINED(eri,CD_x)
    primdata.CD_x[0] = CD[0];
#endif
#if LIBINT2_DEFINED(eri,CD_y)
    primdata.CD_y[0] = CD[1];
#endif
#if LIBINT2_DEFINED(eri,CD_z)
    primdata.CD_z[0] = CD[2];
#endif
#if LIBINT2_DEFINED(eri,DC_x)
    primdata.DC_x[0] = -CD[0];
#endif
#if LIBINT2_DEFINED(eri,DC_y)
    primdata.DC_y[0] = -CD[1];
#endif
#if LIBINT2_DEFINED(eri,DC_z)
    primdata.DC_z[0] = -CD[2];
#endif

    const auto gammap_o_gammapgammaq = oogammapq * gammap;
    const auto gammaq_o_gammapgammaq = oogammapq * gammaq;

    const auto Wx = (gammap_o_gammapgammaq * P[0] + gammaq_o_gammapgammaq * Q[0]);
    const auto Wy = (gammap_o_gammapgammaq * P[1] + gammaq_o_gammapgammaq * Q[1]);
    const auto Wz = (gammap_o_gammapgammaq * P[2] + gammaq_o_gammapgammaq * Q[2]);

#if LIBINT2_DEFINED(eri,WP_x)
    primdata.WP_x[0] = Wx - P[0];
#endif
#if LIBINT2_DEFINED(eri,WP_y)
    primdata.WP_y[0] = Wy - P[1];
#endif
#if LIBINT2_DEFINED(eri,WP_z)
    primdata.WP_z[0] = Wz - P[2];
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
    primdata.WQ_x[0] = Wx - Q[0];
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
    primdata.WQ_y[0] = Wy - Q[1];
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
    primdata.WQ_z[0] = Wz - Q[2];
#endif
#if LIBINT2_DEFINED(eri,oo2z)
    primdata.oo2z[0] = 0.5*oogammap;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
    primdata.oo2e[0] = 0.5*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
    primdata.oo2ze[0] = 0.5*oogammapq;
#endif
#if LIBINT2_DEFINED(eri,roz)
    primdata.roz[0] = rho*oogammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
    primdata.roe[0] = rho*oogammaq;
#endif

    // using ITR?
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
    primdata.TwoPRepITR_pfac0_0_0_x[0] = - (alpha1*AB[0] + alpha3*CD[0]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
    primdata.TwoPRepITR_pfac0_0_0_y[0] = - (alpha1*AB[1] + alpha3*CD[1]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
    primdata.TwoPRepITR_pfac0_0_0_z[0] = - (alpha1*AB[2] + alpha3*CD[2]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
    primdata.TwoPRepITR_pfac0_1_0_x[0] = - (alpha1*AB[0] + alpha3*CD[0]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
    primdata.TwoPRepITR_pfac0_1_0_y[0] = - (alpha1*AB[1] + alpha3*CD[1]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
    primdata.TwoPRepITR_pfac0_1_0_z[0] = - (alpha1*AB[2] + alpha3*CD[2]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
    primdata.TwoPRepITR_pfac0_0_1_x[0] = (alpha0*AB[0] + alpha2*CD[0]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
    primdata.TwoPRepITR_pfac0_0_1_y[0] = (alpha0*AB[1] + alpha2*CD[1]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
    primdata.TwoPRepITR_pfac0_0_1_z[0] = (alpha0*AB[2] + alpha2*CD[2]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
    primdata.TwoPRepITR_pfac0_1_1_x[0] = (alpha0*AB[0] + alpha2*CD[0]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
    primdata.TwoPRepITR_pfac0_1_1_y[0] = (alpha0*AB[1] + alpha2*CD[1]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
    primdata.TwoPRepITR_pfac0_1_1_z[0] = (alpha0*AB[2] + alpha2*CD[2]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,eoz)
    primdata.eoz[0] = gammaq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,zoe)
    primdata.zoe[0] = gammap*oogammaq;
#endif

    // prefactors for derivative ERI relations
    if (deriv_order_ > 0) {
#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
      primdata.alpha1_rho_over_zeta2[0] = alpha0 * (oogammap * gammaq_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
      primdata.alpha2_rho_over_zeta2[0] = alpha1 * (oogammap * gammaq_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
      primdata.alpha3_rho_over_eta2[0] = alpha2 * (oogammaq * gammap_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
      primdata.alpha4_rho_over_eta2[0] = alpha3 * (oogammaq * gammap_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
      primdata.alpha1_over_zetapluseta[0] = alpha0 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
      primdata.alpha2_over_zetapluseta[0] = alpha1 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
      primdata.alpha3_over_zetapluseta[0] = alpha2 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
      primdata.alpha4_over_zetapluseta[0] = alpha3 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
      primdata.rho12_over_alpha1[0] = alpha1 * oogammap;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
      primdata.rho12_over_alpha2[0] = alpha0 * oogammap;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha3)
      primdata.rho34_over_alpha3[0] = alpha3 * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha4)
      primdata.rho34_over_alpha4[0] = alpha2 * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
      primdata.two_alpha0_bra[0] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
      primdata.two_alpha0_ket[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_bra)
      primdata.two_alpha1_bra[0] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_ket)
      primdata.two_alpha1_ket[0] = 2.0 * alpha3;
#endif
    }

    return true;
  }

  template <>
  inline void TwoBodyEngine<DelcGTG_square>::init_core_ints_params(
      const oper_params_type& oparams) {
    // [g12,[- \Del^2, g12] = 2 (\Del g12) \cdot (\Del g12)
    // (\Del exp(-a r_12^2) \cdot (\Del exp(-b r_12^2) = 4 a b (r_{12}^2 exp(- (a+b) r_{12}^2) )
    // i.e. need to scale each coefficient by 4 a b
    const auto ng = oparams.size();
    core_ints_params_.reserve(ng*(ng+1)/2);
    for(size_t b=0; b<ng; ++b)
      for(size_t k=0; k<=b; ++k) {
        const auto gexp = oparams[b].first + oparams[k].first;
        const auto gcoeff = oparams[b].second * oparams[k].second * (b == k ? 1 : 2); // if a != b include ab and ba
        const auto gcoeff_rescaled = 4 * oparams[b].first * oparams[k].first * gcoeff;
        core_ints_params_.push_back(std::make_pair(gexp, gcoeff_rescaled));
      }
  }

  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  inline void TwoBodyEngine<Kernel>::init_core_ints_params(
      const oper_params_type& oparams) {
    core_ints_params_ = oparams;
  }


#endif // LIBINT2_SUPPORT_ERI

#if 1

  /// types of operators (operator sets) supported by Engine.
  /// \warning These must start with 0 and appear in same order as elements of BOOST_PP_ONEBODY_OPERATOR_LIST preprocessor macro.
  /// \warning for the sake of nbody() order operators by # of particles
  enum class Operator {
    overlap=0,        //!< overlap
    kinetic,          //!< electronic kinetic energy, i.e. \f$ -\frac{1}{2} \Nabla^2 \f$
    nuclear,          //!< Coulomb potential due to point charges
    emultipole1,      //!< overlap + (Cartesian) electric dipole moment, \f$ x_O, y_O, z_O \f$, where \f$ x_O \equiv x - O_x \f$ is relative to origin \f$ \vec{O} \f$
    emultipole2,      //!< emultipole1 + (Cartesian) electric quadrupole moment, \f$ x^2, xy, xz, y^2, yz, z^2 \f$
    emultipole3,      //!< emultipole2 + (Cartesian) electric octupole moment, \f$ x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3 \f$
    delta,            //!< \f$ \delta(\vec{r}_1 - \vec{r}_2) \f$
    coulomb,          //!< (2-body) Coulomb operator
    cgtg,             //!< contracted Gaussian geminal
    cgtg_x_coulomb,   //!< contracted Gaussian geminal times Coulomb
    delcgtg2,         //!< |Delta . contracted Gaussian geminal|^2
    invalid=-1, // do not modify this
    // keep this updated
    first_1body_oper=overlap,
    last_1body_oper=emultipole3,
    first_2body_oper=delta,
    last_2body_oper=delcgtg2,
    first_oper=first_1body_oper,
    last_oper=last_2body_oper
  };

  /// @return the particle rank of \c oper
  inline int rank(Operator oper) {
    int n = 0;
    if (oper >= Operator::first_1body_oper && oper <= Operator::last_1body_oper) n = 1;
    else if (oper >= Operator::first_2body_oper && oper <= Operator::last_2body_oper) n = 2;
    return n;
  }

  namespace detail {
    struct default_operator_traits {
        typedef struct {} oper_params_type;
        static oper_params_type default_params() { return oper_params_type{}; }
        static constexpr unsigned int nopers = 1;
        typedef struct _core_eval_type {
            template <typename ...params> static std::shared_ptr<_core_eval_type> instance(params...) { return nullptr; }
        } core_eval_type;
    };
  }
  /// describes operator set \c Op
  /// @tparam Op a value of type Operator
  /// \note default describes operator set of size 1 that takes trivial \c oper_params_type and \c core_eval_type;
  ///       needs to be specialized for some operator types
  template <Operator Op> struct operator_traits : public detail::default_operator_traits {};

  template <> struct operator_traits<Operator::nuclear> : public detail::default_operator_traits {
      /// point charges and their positions
      typedef std::vector<std::pair<double, std::array<double, 3>>> oper_params_type;
      static oper_params_type default_params() { return oper_params_type{}; }
      typedef libint2::FmEval_Taylor<double, 7> core_eval_type;
  };

  template <> struct operator_traits<Operator::emultipole1> : public detail::default_operator_traits {
      /// Cartesian coordinates of the origin with respect to which the dipole moment is defined
      typedef std::array<double, 3> oper_params_type;
      static oper_params_type default_params() { return oper_params_type{{0.0,0.0,0.0}}; }
      static constexpr unsigned int nopers = 4; //!< overlap + 3 dipole components
  };
  template <> struct operator_traits<Operator::emultipole2> : public operator_traits<Operator::emultipole1> {
      static constexpr unsigned int nopers = operator_traits<Operator::emultipole1>::nopers + 6; //!< overlap + 3 dipoles + 6 quadrupoles
  };
  template <> struct operator_traits<Operator::emultipole3> : public operator_traits<Operator::emultipole1> {
      static constexpr unsigned int nopers = operator_traits<Operator::emultipole2>::nopers + 10;
  };

  template <> struct operator_traits<Operator::coulomb> : public detail::default_operator_traits {
      typedef libint2::FmEval_Chebyshev3<double> core_eval_type;
  };
  namespace detail {
    template <int K> struct cgtg_operator_traits : public detail::default_operator_traits {
      typedef libint2::GaussianGmEval<real_t, K> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
    };
  } // namespace detail
  template <> struct operator_traits<Operator::cgtg> : public detail::cgtg_operator_traits<0> {};
  template <> struct operator_traits<Operator::cgtg_x_coulomb> : public detail::cgtg_operator_traits<-1> {};
  template <> struct operator_traits<Operator::delcgtg2> : public detail::cgtg_operator_traits<2> {};

  template <> struct operator_traits<Operator::delta> : public detail::default_operator_traits {
      typedef libint2::GenericGmEval<delta_gm_eval> core_eval_type; // core ints are too trivial to bother
  };

  namespace detail {
    template <typename core_eval_type> using __core_eval_pack_type =
        compressed_pair<std::shared_ptr<core_eval_type>,
                        libint2::detail::CoreEvalScratch<core_eval_type>
        >;
    template <Operator Op> using core_eval_pack_type =
        __core_eval_pack_type<typename operator_traits<Op>::core_eval_type>;
  }

  /// types of shell sets supported by Engine, in chemist notation (i.e. '_' separates particles)
  /// \warning macro \c LIBINT2_MAX_BRAKET_INDEX must equal the maximum value of this enum
  enum class BraKet {
    x_x=0,
    xx_xx,
    xs_xx,
    xs_xs,
    invalid=-1,
    first_1body_braket=x_x,
    last_1body_braket=x_x,
    first_2body_braket=xx_xx,
    last_2body_braket=xs_xs,
    first_braket=first_1body_braket,
    last_braket=last_2body_braket
  };
#define LIBINT2_MAX_BRAKET_INDEX 3

  /// @return rank of \c braket
  inline int rank(BraKet braket) {
    int n = 0;
    switch(braket) {
      case BraKet::x_x:
      case BraKet::xs_xs:
        n = 2; break;
      case BraKet::xs_xx:
        n = 3; break;
      case BraKet::xx_xx:
        n = 4; break;
      default:
        assert(false);
    }
    return n;
  }

  constexpr size_t nopers_2body = (int)Operator::last_2body_oper - (int)Operator::first_2body_oper + 1;
  constexpr size_t nbrakets_2body = (int)BraKet::last_2body_braket - (int)BraKet::first_2body_braket + 1;
  constexpr size_t nderivorders_2body = LIBINT2_MAX_DERIV_ORDER+1;

  template <typename T> struct print_type;

  /**
   * Engine computes integrals of operators (or operator sets) specified by combination of Operator and BraKet.
   * This class deprecates OneBodyEngine and TwoBodyEngine.
   */
  class Engine {

    private:
      typedef struct {} empty_pod;

    public:

/// list of libint task names for each Operator type. These MUST appear in the same order as in Operator
#define BOOST_PP_NBODY_OPERATOR_LIST  (overlap,         \
                                      (kinetic,         \
                                      (elecpot,         \
                                      (1emultipole,     \
                                      (2emultipole,     \
                                      (3emultipole,     \
                                      (eri,             \
                                      (eri,             \
                                      (eri,             \
                                      (eri,             \
                                      (eri,             \
                                       BOOST_PP_NIL)))))))))))

#define BOOST_PP_NBODY_OPERATOR_INDEX_TUPLE BOOST_PP_MAKE_TUPLE( BOOST_PP_LIST_SIZE(BOOST_PP_NBODY_OPERATOR_LIST) )
#define BOOST_PP_NBODY_OPERATOR_INDEX_LIST BOOST_PP_TUPLE_TO_LIST( BOOST_PP_NBODY_OPERATOR_INDEX_TUPLE )

// make list of braket indices for n-body ints
#define BOOST_PP_NBODY_BRAKET_INDEX_TUPLE BOOST_PP_MAKE_TUPLE( BOOST_PP_INC(LIBINT2_MAX_BRAKET_INDEX) )
#define BOOST_PP_NBODY_BRAKET_INDEX_LIST BOOST_PP_TUPLE_TO_LIST( BOOST_PP_NBODY_BRAKET_INDEX_TUPLE )

// make list of derivative orders for n-body ints
#define BOOST_PP_NBODY_DERIV_ORDER_TUPLE BOOST_PP_MAKE_TUPLE( BOOST_PP_INC(LIBINT2_MAX_DERIV_ORDER) )
#define BOOST_PP_NBODY_DERIV_ORDER_LIST BOOST_PP_TUPLE_TO_LIST( BOOST_PP_NBODY_DERIV_ORDER_TUPLE )

      /// creates a default Engine that cannot be used for computing integrals;
      /// to be used as placeholder for copying a usable engine, OR for cleanup of thread-local data
      Engine() :
        oper_(Operator::invalid),
        braket_(BraKet::invalid),
        primdata_(),
        stack_size_(0),
        lmax_(-1) {
        set_precision(std::numeric_limits<real_t>::epsilon());
      }

      /// Constructs a (usable) Engine

      /// \param oper a value of Operator type
      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_order if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_order ,
      ///                    (default=0)
      /// \param precision specifies the target precision with which the integrals will be computed; the default is the "epsilon"
      ///        of \c real_t type, given by \c std::numeric_limits<real_t>::epsilon(). Currently precision control is implemented
      ///        for two-body integrals only. The precision control is somewhat empirical,
      ///        hence be conservative. \sa set_precision()
      /// \param params a value of type Engine::operator_traits<oper>::oper_params_type specifying the parameters of
      ///               the operator set, e.g. position and magnitude of the charges creating the Coulomb potential
      ///               for oper == Operator::nuclear. For most values of \c oper this is not needed.
      ///               \sa Engine::operator_traits
      /// \param braket a value of BraKet type
      /// \warning currently derivative integrals are not supported
      /// \warning currently only one-contraction Shell objects are supported; i.e. generally-contracted Shells are not yet supported
      template <typename Params = empty_pod>
      Engine(Operator oper,
          size_t max_nprim,
          int max_l,
          int deriv_order = 0,
          real_t precision = std::numeric_limits<real_t>::epsilon(),
          Params params = empty_pod(),
          BraKet braket = BraKet::invalid) :
        oper_(oper),
        braket_(braket),
        spbra_(max_nprim), spket_(max_nprim),
        lmax_(max_l),
        deriv_order_(deriv_order),
        params_(enforce_params_type(oper,params))
      {
        set_precision(precision);
        initialize(max_nprim);
        core_eval_pack_ = make_core_eval_pack(oper); // must follow initialize() to ensure default braket_ has been set
        init_core_ints_params(params_);
      }

      /// move constructor
      // intel does not accept "move ctor = default"
      Engine(Engine&& other) :
        oper_(other.oper_),
        braket_(other.braket_),
        primdata_(std::move(other.primdata_)),
        spbra_(std::move(other.spbra_)), spket_(std::move(other.spket_)),
        stack_size_(other.stack_size_),
        lmax_(other.lmax_),
        hard_lmax_(other.hard_lmax_),
        deriv_order_(other.deriv_order_),
        precision_(other.precision_), ln_precision_(other.ln_precision_),
        core_eval_pack_(other.core_eval_pack_),
        params_(std::move(other.params_)),
        core_ints_params_(std::move(other.core_ints_params_)),
        scratch_(std::move(other.scratch_)),
        scratch2_(other.scratch2_),
        buildfnptrs_(other.buildfnptrs_) {
      }

      /// (deep) copy constructor
      Engine(const Engine& other) :
        oper_(other.oper_),
        braket_(other.braket_),
        primdata_(other.primdata_.size()),
        spbra_(other.spbra_), spket_(other.spket_),
        stack_size_(other.stack_size_),
        lmax_(other.lmax_),
        deriv_order_(other.deriv_order_),
        precision_(other.precision_), ln_precision_(other.ln_precision_),
        core_eval_pack_(other.core_eval_pack_),
        params_(other.params_),
        core_ints_params_(other.core_ints_params_)
      {
        initialize();
      }

      ~Engine() {
        finalize();
      }

      /// move assignment is default
      Engine& operator=(Engine&& other) {
        oper_ = other.oper_;
        braket_ = other.braket_;
        primdata_ = std::move(other.primdata_);
        spbra_ = std::move(other.spbra_);
        spket_ = std::move(other.spket_);
        stack_size_ = other.stack_size_;
        lmax_ = other.lmax_;
        hard_lmax_ = other.hard_lmax_;
        deriv_order_ = other.deriv_order_;
        precision_ = other.precision_;
        ln_precision_ = other.ln_precision_;
        core_eval_pack_ = std::move(other.core_eval_pack_);
        params_ = std::move(other.params_);
        core_ints_params_ = std::move(other.core_ints_params_);
        scratch_ = std::move(other.scratch_);
        scratch2_ = other.scratch2_;
        buildfnptrs_ = other.buildfnptrs_;
        return *this;
      }

      /// (deep) copy assignment
      Engine& operator=(const Engine& other) {
        oper_ = other.oper_;
        braket_ = other.braket_;
        primdata_.resize(other.primdata_.size());
        spbra_ = other.spbra_;
        spket_ = other.spket_;
        stack_size_ = other.stack_size_;
        lmax_ = other.lmax_;
        deriv_order_ = other.deriv_order_;
        precision_ = other.precision_;
        ln_precision_ = other.ln_precision_;
        core_eval_pack_ = other.core_eval_pack_;
        params_ = other.params_;
        core_ints_params_ = other.core_ints_params_;
        initialize();
        return *this;
      }

      /// returns the particle rank of the operator
      int operator_rank() const {
        return rank(oper_);
      }

      /// rank of the braket
      int braket_rank() const {
        return rank(braket_);
      }

      /// resets operator type
      void set_oper(Operator new_oper) {
        if (rank(new_oper) != operator_rank())
          braket_ = BraKet::invalid;
        oper_ = new_oper;
        initialize();
      }

      /// resets braket type
      void set_braket(BraKet new_braket) {
        braket_ = new_braket;
        initialize();
      }

      /// resets operator parameters; this may be useful e.g. if need to compute Coulomb potential
      /// integrals over batches of charges for the sake of parallelism.
      template <typename Params>
      void set_params(const Params& params) {
        params_ = params;
        init_core_ints_params(params_);
        reset_scratch();
      }

      /// reports the number of shell sets that each call to compute() produces.
      /// this depends on the order of geometrical derivatives requested and
      /// on the operator set. \sa compute()
      /// \note need to specialize for some operator types
      unsigned int nshellsets() const {
        const unsigned int num_operator_geometrical_derivatives = (oper_ == Operator::nuclear) ? this->nparams() : 0;
        const auto ncenters = braket_rank() + num_operator_geometrical_derivatives;
        return nopers() * num_geometrical_derivatives(ncenters,deriv_order_);
      }

      /// Given the Cartesian geometric derivative index that refers to center set
      /// (0...n-1) with one center omitted compute the derivative index referring
      /// to the full center set
      /// \param deriv_idx index of the derivative referring to the set with center \c omitted_center omitted
      /// \param deriv_order order of the geometric derivative
      /// \param ncenters number of centers in the full set
      /// \param omitted_center the omitted center, must be less than \c ncenters.
      /// \return the index of the derivative referring to full set
      static unsigned int to_target_deriv_index(unsigned int deriv_idx,
                                                unsigned int deriv_order,
                                                unsigned int ncenters,
                                                unsigned int omitted_center) {
        auto ncenters_reduced = ncenters - 1;
        auto nderiv_1d = 3 * ncenters_reduced;
        assert(deriv_idx < num_geometrical_derivatives(ncenters_reduced,deriv_order));
        switch (deriv_order) {
          case 1: return deriv_idx >= omitted_center*3 ? deriv_idx+3 : deriv_idx;
          default: assert(deriv_order<=1); // not implemented, won't be needed when Libint computes all derivatives
        }
        assert(false); // unreachable
        return 0;
      }

      /// computes shell set of integrals
      /// \note result is stored in row-major order
      template <typename ... ShellPack>
      const real_t* compute(const libint2::Shell& first_shell, const ShellPack&... rest_of_shells) {
        constexpr size_t nargs = 1 + sizeof...(rest_of_shells);
        assert(nargs == braket_rank());

        std::array<std::reference_wrapper<const Shell>, nargs> shells{first_shell, rest_of_shells...};

        if (operator_rank() == 1) {
          if (nargs == 2)
            return compute1(shells[0], shells[1]);
        } else if (operator_rank() == 2) {
          auto compute_ptr_idx = (((int)oper_ - (int)Operator::first_2body_oper) * nbrakets_2body + ((int)braket_ - (int)BraKet::first_2body_braket)) * nderivorders_2body + deriv_order_;
          auto compute_ptr = compute2_ptrs()[compute_ptr_idx];
          assert(compute_ptr != nullptr);
          if (nargs == 2)
            return (this->*compute_ptr)(shells[0], Shell::unit(), shells[1], Shell::unit());
          if (nargs == 3)
            return (this->*compute_ptr)(shells[0], Shell::unit(), shells[1], shells[2]);
          if (nargs == 4)
            return (this->*compute_ptr)(shells[0], shells[1], shells[2], shells[3]);
        }

        assert(false); // only reached if missing a feature
        return nullptr;
      }

      /// computes shell set of integrals
      /// \note result is stored in row-major order
      const real_t* compute1(const libint2::Shell& s1,
                             const libint2::Shell& s2) {

        // can only handle 1 contraction at a time
        assert(s1.ncontr() == 1 && s2.ncontr() == 1);

        const auto l1 = s1.contr[0].l;
        const auto l2 = s2.contr[0].l;

        // if want nuclear, make sure there is at least one nucleus .. otherwise the user likely forgot to call set_params
        if (oper_ == Operator::nuclear and nparams() == 0)
          throw std::runtime_error("Engine<nuclear>, but no charges found; forgot to call set_params()?");

        const auto n1 = s1.size();
        const auto n2 = s2.size();
        const auto n12 = n1 * n2;
        const auto ncart1 = s1.cartesian_size();
        const auto ncart2 = s2.cartesian_size();
        const auto ncart12 = ncart1 * ncart2;
        const auto nops = nopers();

        const auto tform_to_solids = s1.contr[0].pure || s2.contr[0].pure;

        // assert # of primitive pairs
        const auto nprim1 = s1.nprim();
        const auto nprim2 = s2.nprim();
        const auto nprimpairs = nprim1 * nprim2;
        assert(nprimpairs <= primdata_.size());

        auto nparam_sets = nparams();

        // how many shell sets will be returned?
        auto num_shellsets = nshellsets();
        // Libint computes derivatives with respect to one center fewer, will use translational invariance to recover
        const auto geometry_independent_operator = not (oper_ == Operator::nuclear);
        const auto num_deriv_centers_computed = geometry_independent_operator ? braket_rank()-1 : braket_rank();
        auto num_shellsets_computed = nopers() *
                                      num_geometrical_derivatives(num_deriv_centers_computed,
                                                                  deriv_order_);
        // size of ints block computed by Libint
        const auto target_buf_size = num_shellsets_computed * ncart12;

        // will use scratch_ if:
        // - Coulomb ints are computed 1 charge at a time, contributions are accumulated in scratch_ (unless la==lb==0)
        // - derivatives on the missing center need to be reconstructed (no need to accumulate into scratch though)
        // will only use scratch to accumulate ints when
        const auto accumulate_ints_in_scratch = (oper_ == Operator::nuclear);

        // adjust max angular momentum, if needed
        const auto lmax = std::max(l1, l2);
        assert (lmax <= lmax_);

        // where cartesian ints are located varies, sometimes we compute them in scratch, etc.
        // this is the most likely location
        auto cartesian_ints = primdata_[0].stack;

        // simple (s|s) ints will be computed directly and accumulated in the first element of stack
        const auto compute_directly = lmax == 0 && deriv_order_ == 0 && (oper_ == Operator::overlap || oper_ == Operator::nuclear);
        if (compute_directly) {
          primdata_[0].stack[0] = 0;
        }
        else if (accumulate_ints_in_scratch)
          memset(static_cast<void*>(&scratch_[0]), 0, sizeof(real_t)*target_buf_size);

        // loop over accumulation batches
        for(auto pset=0u; pset!=nparam_sets; ++pset) {

          if (oper_!=Operator::nuclear) assert(nparam_sets == 1);

          auto p12 = 0;
          for(auto p1=0; p1!=nprim1; ++p1) {
            for(auto p2=0; p2!=nprim2; ++p2, ++p12) {
              compute_primdata(primdata_[p12],s1,s2,p1,p2,pset);
            }
          }
          primdata_[0].contrdepth = p12;

          if (compute_directly) {
            auto& result = cartesian_ints[0];
            switch (oper_) {
              case Operator::overlap:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12]._0_Overlap_0_x[0]
                          * primdata_[p12]._0_Overlap_0_y[0]
                          * primdata_[p12]._0_Overlap_0_z[0];
                  break;
              case Operator::nuclear:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12].LIBINT_T_S_ELECPOT_S(0)[0];
                  break;
              default:
                assert(false);
            }
            primdata_[0].targets[0] = cartesian_ints;
          }
          else {

            buildfnptrs_[s1.contr[0].l*hard_lmax_ + s2.contr[0].l](&primdata_[0]);
            cartesian_ints = primdata_[0].targets[0];

            if (accumulate_ints_in_scratch) {
              cartesian_ints = &scratch_[0];
              std::transform(primdata_[0].targets[0], primdata_[0].targets[0] + target_buf_size,
                             &scratch_[0],
                             &scratch_[0], std::plus<real_t>());

              // need to reconstruct derivatives of nuclear ints for each nucleus
              if (deriv_order_ > 0){
                const auto nints_per_center = target_buf_size/2;
                // first two blocks are derivatives with respect to Gaussian positions
                // rest are derivs with respect to nuclear coordinates
                auto dest = &scratch_[0] + (2+pset)*nints_per_center;
                auto src = primdata_[0].targets[0];
                for(auto i=0; i!=nints_per_center; ++i) {
                  dest[i] = -src[i];
                }
                src = primdata_[0].targets[0] + nints_per_center;
                for(auto i=0; i!=nints_per_center; ++i) {
                  dest[i] -= src[i];
                }
                num_shellsets_computed+=3; // we just added 3 shell sets
              } // reconstruct derivatives

            }
          } // ltot != 0

        } // pset (accumulation batches)

        auto result = cartesian_ints; // will be adjusted as we proceed
        if (tform_to_solids) {
          // where do spherical ints go?
          auto* spherical_ints = (cartesian_ints == &scratch_[0]) ? scratch2_ : &scratch_[0];
          result = spherical_ints;

          // transform to solid harmonics, one shell set at a time:
          // for each computed shell set ...
          for(auto s=0ul; s!=num_shellsets_computed; ++s, cartesian_ints+=ncart12) {
            // ... find its index in the target shell set:
            // 1. if regular ints do nothing
            // 2. for derivatives the target set includes derivatives w.r.t omitted centers,
            //    to be computed later (for all ints) or already computed (for nuclear);
            //    in the former case the "omitted" set of derivatives always comes at the end
            //    hence the index of the current shell set does not change (for 2-body ints
            //    the rules are different, but Libint will eliminate the need to reconstruct via
            //    translational invariance soon, so this logic will be unnecessary).
            auto s_target = s;
            // .. and compute the destination
            spherical_ints = result + n12 * s_target;
            if (s1.contr[0].pure && s2.contr[0].pure) {
              libint2::solidharmonics::tform(l1, l2, cartesian_ints, spherical_ints);
            }
            else {
              if (s1.contr[0].pure)
                libint2::solidharmonics::tform_rows(l1, n2, cartesian_ints, spherical_ints);
              else
                libint2::solidharmonics::tform_cols(n1, l2, cartesian_ints, spherical_ints);
            }
          } // loop cartesian shell set

        } // tform to solids

        // if computing derivatives of ints of geometry-independent operators
        // compute the omitted derivatives using translational invariance
        if (deriv_order_ > 0 && geometry_independent_operator) {
          assert(deriv_order_ == 1); // assuming 1st-order derivs here, arbitrary order later

          const auto nints_computed = n12*num_shellsets_computed; // target # of ints is twice this

          // make sure there is enough room left in libint stack
          // if not, copy into scratch2_
          if (not tform_to_solids) {
            const auto stack_size_remaining = stack_size_ - (result-primdata_[0].stack) - nints_computed;
            const auto copy_to_scratch2 = stack_size_remaining < nints_computed;
            if (copy_to_scratch2) {
              // this is tricky ... copy does not allow scratch2_ in [result, result + nints_computed)
              // but this would only happen in scratch2_ == result, but definition of scratch2_ ensures this
              std::copy(result, result + nints_computed, scratch2_);
              result = scratch2_;
            }
          }

          const auto src = result;
          const auto dest = result + nints_computed;
          for(auto f=0ul; f!=nints_computed; ++f) {
            dest[f] = -src[f];
          }
        } // rebuild omitted derivatives of Cartesian ints

        return result;
      }

      template <Operator oper, BraKet braket, size_t deriv_order>
      inline const real_t* compute2(const Shell& s1,
                                    const Shell& s2,
                                    const Shell& s3,
                                    const Shell& s4);

      typedef const real_t* (Engine::*compute2_ptr_type)(const Shell& s1,
          const Shell& s2,
          const Shell& s3,
          const Shell& s4);

      /** this specifies target precision for computing the integrals.
       *  target precision \f$ \epsilon \f$ is used in 3 ways:
       *  (1) to screen out primitive pairs in ShellPair object for which
       *      \f$ {\rm scr}_{12} = \max|c_1| \max|c_2| \exp(-\rho_{12} |AB|^2)/\gamma_{12} < \epsilon \f$ ;
       *  (2) to screen out primitive quartets outside compute_primdata() for which \f$ {\rm scr}_{12} {\rm scr}_{34} <  \epsilon \f$;
       *  (3) to screen out primitive quartets inside compute_primdata() for which the prefactor of \f$ F_m(\rho, T) \f$ is smaller
       *      than \f$ \epsilon \f$ .
       */
      void set_precision(real_t prec) {
        if (prec <= 0.) {
          precision_ = 0.;
          ln_precision_ = std::numeric_limits<real_t>::lowest();
        }
        else {
          precision_ = prec;
          ln_precision_ = std::log(precision_);
        }
      }
      /// @return the target precision for computing the integrals
      /// @sa set_precision(real_t)
      real_t precision() const {
        return precision_;
      }

      void print_timers() {
#ifdef LIBINT2_ENGINE_TIMERS
        std::cout << "timers: prereq = " << timers.read(0);
#  ifndef LIBINT2_PROFILE // if libint's profiling was on, engine's build timer will include its overhead
                          // do not report it, detailed profiling from libint will be printed below
        std::cout << " build = " << timers.read(1);
#  endif
        std::cout << " tform = " << timers.read(2) << std::endl;
#endif
#ifdef LIBINT2_PROFILE
        std::cout << "build timers: hrr = " << primdata_[0].timers->read(0)
                << " vrr = " << primdata_[0].timers->read(1) << std::endl;
#endif
#ifdef LIBINT2_ENGINE_TIMERS
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
        for(const auto& p: class_profiles) {
          printf("{\"%s\", %10.5lf, %10.5lf, %10.5lf, %10.5lf, %ld, %ld},\n",
                 p.first.to_string().c_str(),
                 p.second.prereqs,
                 p.second.build_vrr,
                 p.second.build_hrr,
                 p.second.tform,
                 p.second.nshellset,
                 p.second.nprimset);
        }
#  endif
#endif
      }

    private:
      Operator oper_;
      BraKet braket_;
      std::vector<Libint_t> primdata_;
      ShellPair spbra_, spket_;
      size_t stack_size_;  // amount allocated by libint2_init_xxx in primdata_[0].stack
      int lmax_;
      int hard_lmax_;      // max L supported by library for this operator type + 1
      size_t deriv_order_;
      real_t precision_;
      real_t ln_precision_;

      any core_eval_pack_;

      /// operator params
      any params_; // operator params
      /// for some operators need core ints params that are computed from operator params,
      /// e.g. integrals of \f$ f_{12}^2 \f$ are computed from parameters of \f$ f_{12} \f$
      any core_ints_params_;
      /// makes core ints params from the operator params
      void init_core_ints_params(const any& params);

      std::vector<real_t> scratch_; // for transposes and/or transforming to solid harmonics
      real_t* scratch2_;            // &scratch_[0] points to the first block large enough to hold all target ints
                                    // scratch2_ points to second such block. It could point into scratch_ or at primdata_[0].stack
      typedef void (*buildfnptr_t)(const Libint_t*);
      buildfnptr_t* buildfnptrs_;

      void reset_scratch() {
        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;
        const auto target_shellset_size = nshellsets() * std::pow(ncart_max, braket_rank());
        // need to be able to hold 2 sets of target shellsets: the worst case occurs when dealing with
        // 1-body Coulomb ints derivatives ... have 2+natom derivative sets that are stored in scratch
        // then need to transform to solids. To avoid copying back and forth make sure that there is enough
        // room to transform all ints and save them in correct order in single pass
        const auto need_extra_large_scratch = stack_size_ < target_shellset_size;
        scratch_.resize(need_extra_large_scratch ? 2*target_shellset_size : target_shellset_size);
        scratch2_ = need_extra_large_scratch ? &scratch_[target_shellset_size] : primdata_[0].stack;
      }

      inline void compute_primdata(Libint_t& primdata,
                                   const Shell& s1, const Shell& s2,
                                   size_t p1, size_t p2,
                                   size_t oset);

      /// 3-dim array of pointers to help dispatch efficiently based on oper_, braket_, and deriv_order_
      inline const std::vector<Engine::compute2_ptr_type>& compute2_ptrs() const;

      void initialize(size_t max_nprim = 0) {
        assert(libint2::initialized());
        assert(deriv_order_ <= LIBINT2_MAX_DERIV_ORDER);

        // initialize braket, if needed
        if (braket_ == BraKet::invalid) {
          switch(operator_rank()) {
            case 1: braket_ = BraKet::x_x; break;
            case 2: braket_ = BraKet::xx_xx; break;
            default: assert(false);
          }
        }

        if (max_nprim != 0)
          primdata_.resize(std::pow(max_nprim, braket_rank()));

#ifdef LIBINT2_ENGINE_TIMERS
        timers.set_now_overhead(25);
#endif
#ifdef LIBINT2_PROFILE
        primdata_[0].timers->set_now_overhead(25);
#endif

#define BOOST_PP_NBODYENGINE_MCR2(r,product)                                                                  \
         if (static_cast<int>(oper_) == BOOST_PP_TUPLE_ELEM(2,0,product) && deriv_order_ == BOOST_PP_TUPLE_ELEM(2,1,product) ) {\
           assert(lmax_ <= BOOST_PP_CAT(LIBINT2_MAX_AM_ ,                                                     \
                                        BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST,                        \
                                                         BOOST_PP_TUPLE_ELEM(2,0,product) )                   \
                                       ) );                                                                   \
           stack_size_ =                                                                                      \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_need_memory_ ,                                       \
                                      BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST,                          \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                 )(lmax_);                                                                    \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_init_ ,                                              \
                                      BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST,                          \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                )(&primdata_[0], lmax_, 0);                                                   \
           buildfnptrs_ = to_ptr1(                                                                            \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_build_ ,                                             \
                                      BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST,                          \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                )                                                                             \
                                 );                                                                           \
           hard_lmax_ =           BOOST_PP_CAT(                                                               \
                                    LIBINT2_MAX_AM_ ,                                                         \
                                    BOOST_PP_CAT(                                                             \
                                      BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST,                          \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) ),                    \
                                      BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),     \
                                                    BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()        \
                                                )                                                             \
                                    )                                                                         \
                                  ) + 1;                                                                      \
           reset_scratch();                                                                                   \
           return;                                                                                            \
         }

BOOST_PP_LIST_FOR_EACH_PRODUCT ( BOOST_PP_NBODYENGINE_MCR2, 2, (BOOST_PP_NBODY_OPERATOR_INDEX_LIST, BOOST_PP_NBODY_DERIV_ORDER_LIST) )

        assert(false); // either deriv_order_ or oper_ is wrong
      } // initialize()

      void finalize() {
        if (primdata_.size() != 0) {

#define BOOST_PP_NBODYENGINE_MCR3(r,product)                                                                  \
           LIBINT2_PREFIXED_NAME( BOOST_PP_CAT(                                                               \
                                    BOOST_PP_CAT(libint2_cleanup_ ,                                           \
                                      BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST,                          \
                                                       BOOST_PP_TUPLE_ELEM(2,0,product) )                     \
                                    ),                                                                        \
                                    BOOST_PP_IIF( BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(2,1,product),0),       \
                                                  BOOST_PP_TUPLE_ELEM(2,1,product), BOOST_PP_EMPTY()          \
                                                )                                                             \
                                  )                                                                           \
                                )(&primdata_[0]);                                                             \
           return;

BOOST_PP_LIST_FOR_EACH_PRODUCT ( BOOST_PP_NBODYENGINE_MCR3, 2, (BOOST_PP_NBODY_OPERATOR_INDEX_LIST, BOOST_PP_NBODY_DERIV_ORDER_LIST) )

        }
      } // finalize()

      //-------
      // utils
      //-------
      inline unsigned int nparams() const;
      inline unsigned int nopers() const;
      /// if Params == operator_traits<oper>::oper_params_type, will return any(params)
      /// else will set return any initialized with default value for operator_traits<type>::oper_params_type
      /// @param throw_if_wrong_type if true, and Params != operator_traits<type>::oper_params_type, will throw std::bad_cast
      template <typename Params>
      static any enforce_params_type(Operator oper,
                                     const Params& params,
                                     bool throw_if_wrong_type = not std::is_same<Params,empty_pod>::value);
      /// @return core eval pack corresponding to operator_traits<oper>::core_eval_type
      any make_core_eval_pack(Operator oper) const;

      //-------
      // profiling
      //-------
      static const bool skip_core_ints = false;

  }; // struct Engine

  namespace detail {
    inline std::vector<Engine::compute2_ptr_type> init_compute2_ptrs() {
      auto max_ncompute2_ptrs = nopers_2body * nbrakets_2body * nderivorders_2body;
      std::vector<Engine::compute2_ptr_type> result(max_ncompute2_ptrs, nullptr);

#define BOOST_PP_NBODYENGINE_MCR7(r,product)                                                                                  \
    if (BOOST_PP_TUPLE_ELEM(3,0,product) >= (int)Operator::first_2body_oper &&                                                \
        BOOST_PP_TUPLE_ELEM(3,0,product) <= (int)Operator::last_2body_oper &&                                                 \
        BOOST_PP_TUPLE_ELEM(3,1,product) >= (int)BraKet::first_2body_braket &&                                                \
        BOOST_PP_TUPLE_ELEM(3,1,product) <= (int)BraKet::last_2body_braket )                                                  \
      {  auto compute_ptr_idx = ((BOOST_PP_TUPLE_ELEM(3,0,product) - (int)Operator::first_2body_oper) * nbrakets_2body        \
                                 + (BOOST_PP_TUPLE_ELEM(3,1,product) - (int)BraKet::first_2body_braket)) * nderivorders_2body \
                                + BOOST_PP_TUPLE_ELEM(3,2,product) ;                                                          \
         result.at(compute_ptr_idx) = &Engine::compute2<static_cast<Operator>(BOOST_PP_TUPLE_ELEM(3,0,product)),              \
                                                     static_cast<BraKet>(BOOST_PP_TUPLE_ELEM(3,1,product)),                   \
                                                     BOOST_PP_TUPLE_ELEM(3,2,product)>;  }                                    \


BOOST_PP_LIST_FOR_EACH_PRODUCT ( BOOST_PP_NBODYENGINE_MCR7, 3, (BOOST_PP_NBODY_OPERATOR_INDEX_LIST, BOOST_PP_NBODY_BRAKET_INDEX_LIST, BOOST_PP_NBODY_DERIV_ORDER_LIST) )

      return result;
    }
  } // namespace detail

  inline const std::vector<Engine::compute2_ptr_type>& Engine::compute2_ptrs() const {
    static std::vector<compute2_ptr_type> compute2_ptrs_ = detail::init_compute2_ptrs();
    return compute2_ptrs_;
  }

  inline unsigned int Engine::nparams() const {
    switch (oper_) {
      case Operator::nuclear:
        return params_.as<operator_traits<Operator::nuclear>::oper_params_type>().size();
      default:
        return 1;
    }
    return 1;
  }
  inline unsigned int Engine::nopers() const {
    switch (static_cast<int>(oper_)) {
#define BOOST_PP_NBODYENGINE_MCR4(r,data,i,elem)  case i : return operator_traits< static_cast<Operator> ( i ) >::nopers;
BOOST_PP_LIST_FOR_EACH_I ( BOOST_PP_NBODYENGINE_MCR4, _, BOOST_PP_NBODY_OPERATOR_LIST)
      default: break;
    }
    assert(false); // unreachable
    return 0;
  }

  template <typename Params>
  inline any Engine::enforce_params_type(Operator oper,
                                  const Params& params,
                                  bool throw_if_wrong_type) {
    any result;
    switch(static_cast<int>(oper)) {

#define BOOST_PP_NBODYENGINE_MCR5(r,data,i,elem)                                                             \
      case i :                                                                                               \
      if (std::is_same<Params,operator_traits< static_cast<Operator> ( i ) >::oper_params_type>::value)      \
        result = params;                                                                                     \
      else {                                                                                                 \
        if (throw_if_wrong_type) throw std::bad_cast();                                                      \
        result = operator_traits<static_cast<Operator> ( i ) >::default_params();                            \
      }                                                                                                      \
      break;

BOOST_PP_LIST_FOR_EACH_I ( BOOST_PP_NBODYENGINE_MCR5, _, BOOST_PP_NBODY_OPERATOR_LIST)

      default:
        assert(false); // missed a case?
    }
    return result;
  }

  inline any Engine::make_core_eval_pack(Operator oper) const {
    any result;
    switch(static_cast<int>(oper)) {

#define BOOST_PP_NBODYENGINE_MCR6(r,data,i,elem)                                                             \
      case i :                                                                                               \
      result = libint2::detail::make_compressed_pair(                                                        \
          operator_traits<static_cast<Operator> ( i ) >::core_eval_type::instance(braket_rank()*lmax_ + deriv_order_,    \
                                                                                  std::numeric_limits<real_t>::epsilon()), \
          libint2::detail::CoreEvalScratch<operator_traits<static_cast<Operator> ( i ) >::core_eval_type>(braket_rank()*lmax_ + deriv_order_) \
      );                                                                                                     \
      break;

BOOST_PP_LIST_FOR_EACH_I ( BOOST_PP_NBODYENGINE_MCR6, _, BOOST_PP_NBODY_OPERATOR_LIST)

      default:
        assert(false); // missed a case?
    }
    return result;
  }

  inline void Engine::init_core_ints_params(const any& params) {
    if (oper_ == Operator::delcgtg2) {
      // [g12,[- \Del^2, g12] = 2 (\Del g12) \cdot (\Del g12)
      // (\Del exp(-a r_12^2) \cdot (\Del exp(-b r_12^2) = 4 a b (r_{12}^2 exp(- (a+b) r_{12}^2) )
      // i.e. need to scale each coefficient by 4 a b
      auto oparams = params.as<operator_traits<Operator::delcgtg2>::oper_params_type>();
      const auto ng = oparams.size();
      operator_traits<Operator::delcgtg2>::oper_params_type core_ints_params;
      core_ints_params.reserve(ng*(ng+1)/2);
      for(size_t b=0; b<ng; ++b)
        for(size_t k=0; k<=b; ++k) {
          const auto gexp = oparams[b].first + oparams[k].first;
          const auto gcoeff = oparams[b].second * oparams[k].second * (b == k ? 1 : 2); // if a != b include ab and ba
          const auto gcoeff_rescaled = 4 * oparams[b].first * oparams[k].first * gcoeff;
          core_ints_params.push_back(std::make_pair(gexp, gcoeff_rescaled));
        }
      core_ints_params_ = core_ints_params;
    }
    else
      core_ints_params_ = params;
  }

  inline void Engine::compute_primdata(Libint_t& primdata,
                                       const Shell& s1, const Shell& s2,
                                       size_t p1, size_t p2,
                                       size_t oset) {

    const auto& A = s1.O;
    const auto& B = s2.O;

    const auto alpha1 = s1.alpha[p1];
    const auto alpha2 = s2.alpha[p2];

    const auto c1 = s1.contr[0].coeff[p1];
    const auto c2 = s2.contr[0].coeff[p2];

    const auto gammap = alpha1 + alpha2;
    const auto oogammap = 1.0 / gammap;
    const auto rhop_over_alpha1 = alpha2 * oogammap;
    const auto rhop = alpha1 * rhop_over_alpha1;
    const auto Px = (alpha1 * A[0] + alpha2 * B[0]) * oogammap;
    const auto Py = (alpha1 * A[1] + alpha2 * B[1]) * oogammap;
    const auto Pz = (alpha1 * A[2] + alpha2 * B[2]) * oogammap;
    const auto AB_x = A[0] - B[0];
    const auto AB_y = A[1] - B[1];
    const auto AB_z = A[2] - B[2];
    const auto AB2_x = AB_x * AB_x;
    const auto AB2_y = AB_y * AB_y;
    const auto AB2_z = AB_z * AB_z;
    const auto AB2 = AB2_x + AB2_y + AB2_z;

    assert (LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD);

    // overlap and kinetic energy ints don't use HRR, hence VRR on both centers
    // Coulomb potential do HRR on center 1 only
#if LIBINT2_DEFINED(eri,PA_x)
      primdata.PA_x[0] = Px - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
      primdata.PA_y[0] = Py - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
      primdata.PA_z[0] = Pz - A[2];
#endif

    if (oper_ != Operator::nuclear) {

#if LIBINT2_DEFINED(eri,PB_x)
      primdata.PB_x[0] = Px - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
      primdata.PB_y[0] = Py - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
      primdata.PB_z[0] = Pz - B[2];
#endif
    }

    if (oper_ == Operator::emultipole1 ||
        oper_ == Operator::emultipole2 ||
        oper_ == Operator::emultipole3) {
      auto& O = params_.as<operator_traits<Operator::emultipole1>::oper_params_type>(); // same as emultipoleX
#if LIBINT2_DEFINED(eri,BO_x)
      primdata.BO_x[0] = B[0] - O[0];
#endif
#if LIBINT2_DEFINED(eri,BO_y)
      primdata.BO_y[0] = B[1] - O[1];
#endif
#if LIBINT2_DEFINED(eri,BO_z)
      primdata.BO_z[0] = B[2] - O[2];
#endif
    }

#if LIBINT2_DEFINED(eri,oo2z)
    primdata.oo2z[0] = 0.5*oogammap;
#endif

    if (oper_ == Operator::nuclear) { // additional factor for electrostatic potential
      auto& params = params_.as<operator_traits<Operator::nuclear>::oper_params_type>();
      const auto& C = params[oset].second;
#if LIBINT2_DEFINED(eri,PC_x)
      primdata.PC_x[0] = Px - C[0];
#endif
#if LIBINT2_DEFINED(eri,PC_y)
      primdata.PC_y[0] = Py - C[1];
#endif
#if LIBINT2_DEFINED(eri,PC_z)
      primdata.PC_z[0] = Pz - C[2];
#endif
      // elecpot uses HRR
#if LIBINT2_DEFINED(eri,AB_x)
      primdata.AB_x[0] = A[0] - B[0];
#endif
#if LIBINT2_DEFINED(eri,AB_y)
      primdata.AB_y[0] = A[1] - B[1];
#endif
#if LIBINT2_DEFINED(eri,AB_z)
      primdata.AB_z[0] = A[2] - B[2];
#endif

    }

    decltype(c1) sqrt_PI(1.77245385090551602729816748334);
    const auto xyz_pfac = sqrt_PI * sqrt(oogammap);
    const auto ovlp_ss_x = exp(- rhop * AB2_x) * xyz_pfac * c1 * c2;
    const auto ovlp_ss_y = exp(- rhop * AB2_y) * xyz_pfac;
    const auto ovlp_ss_z = exp(- rhop * AB2_z) * xyz_pfac;

    primdata._0_Overlap_0_x[0] = ovlp_ss_x;
    primdata._0_Overlap_0_y[0] = ovlp_ss_y;
    primdata._0_Overlap_0_z[0] = ovlp_ss_z;

    if (oper_ == Operator::kinetic || (deriv_order_ > 0)) {
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
      primdata.two_alpha0_bra[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
      primdata.two_alpha0_ket[0] = 2.0 * alpha2;
#endif
    }

    if (oper_ == Operator::nuclear) {
#if LIBINT2_DEFINED(eri,rho12_over_alpha1) || LIBINT2_DEFINED(eri,rho12_over_alpha2)
      if (deriv_order_ > 0) {
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
        primdata.rho12_over_alpha1[0] = rhop_over_alpha1;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
        primdata.rho12_over_alpha2[0] = alpha1 * oogammap;
#endif
      }
#endif
#if LIBINT2_DEFINED(eri,PC_x) && LIBINT2_DEFINED(eri,PC_y) && LIBINT2_DEFINED(eri,PC_z)
      const auto PC2 = primdata.PC_x[0] * primdata.PC_x[0] +
                       primdata.PC_y[0] * primdata.PC_y[0] +
                       primdata.PC_z[0] * primdata.PC_z[0];
      const auto U = gammap * PC2;
      const auto mmax = s1.contr[0].l + s2.contr[0].l + deriv_order_;
      auto* fm_ptr = &(primdata.LIBINT_T_S_ELECPOT_S(0)[0]);
      auto fm_engine_ptr = core_eval_pack_.as<detail::core_eval_pack_type<Operator::nuclear>>().first();
      fm_engine_ptr->eval(fm_ptr, U, mmax);

      decltype(U) two_o_sqrt_PI(1.12837916709551257389615890312);
      const auto q = params_.as<operator_traits<Operator::nuclear>::oper_params_type>()[oset].first;
      const auto pfac = - q * sqrt(gammap) * two_o_sqrt_PI * ovlp_ss_x * ovlp_ss_y * ovlp_ss_z;
      const auto m_fence = mmax + 1;
      for(auto m=0; m!=m_fence; ++m) {
        fm_ptr[m] *= pfac;
      }
#endif
    }

  } // Engine::compute_primdata()


  /// computes shell set of integrals of 2-body operator
  /// \note result is stored in the "chemists" form, i.e. (tbra1 tbra2 |tket1 tket2), in row-major order
  template <Operator oper, BraKet braket, size_t deriv_order>
  inline const real_t*
  Engine::compute2(const libint2::Shell& tbra1,
                   const libint2::Shell& tbra2,
                   const libint2::Shell& tket1,
                   const libint2::Shell& tket2) {

    //
    // i.e. bra and ket refer to chemists bra and ket
    //

    // can only handle 1 contraction at a time
    assert(tbra1.ncontr() == 1 && tbra2.ncontr() == 1 &&
           tket1.ncontr() == 1 && tket2.ncontr() == 1);

    // derivatives not supported for now
    assert(deriv_order == 0);

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // standard angular momentum ordering
    auto swap_bra = (tbra1.contr[0].l < tbra2.contr[0].l);
    auto swap_ket = (tket1.contr[0].l < tket2.contr[0].l);
    auto swap_braket = (tbra1.contr[0].l + tbra2.contr[0].l > tket1.contr[0].l + tket2.contr[0].l);
#else // orca angular momentum ordering
    auto swap_bra = (tbra1.contr[0].l > tbra2.contr[0].l);
    auto swap_ket = (tket1.contr[0].l > tket2.contr[0].l);
    auto swap_braket = (tbra1.contr[0].l + tbra2.contr[0].l < tket1.contr[0].l + tket2.contr[0].l);
#endif
    const auto& bra1 = swap_braket ? (swap_ket ? tket2 : tket1) : (swap_bra ? tbra2 : tbra1);
    const auto& bra2 = swap_braket ? (swap_ket ? tket1 : tket2) : (swap_bra ? tbra1 : tbra2);
    const auto& ket1 = swap_braket ? (swap_bra ? tbra2 : tbra1) : (swap_ket ? tket2 : tket1);
    const auto& ket2 = swap_braket ? (swap_bra ? tbra1 : tbra2) : (swap_ket ? tket1 : tket2);

    const bool tform = tbra1.contr[0].pure || tbra2.contr[0].pure || tket1.contr[0].pure || tket2.contr[0].pure;
    const bool use_scratch = (swap_braket || swap_bra || swap_ket || tform);

    // assert # of primitive pairs
    auto nprim_bra1 = bra1.nprim();
    auto nprim_bra2 = bra2.nprim();
    auto nprim_ket1 = ket1.nprim();
    auto nprim_ket2 = ket2.nprim();

    // adjust max angular momentum, if needed
    auto lmax = std::max(std::max(bra1.contr[0].l, bra2.contr[0].l), std::max(ket1.contr[0].l, ket2.contr[0].l));
    assert (lmax <= lmax_);
    if (lmax == 0) // (ss|ss) ints will be accumulated in the first element of stack
      primdata_[0].stack[0] = 0;

#ifdef LIBINT2_ENGINE_PROFILE_CLASS
    class_id id(bra1.contr[0].l, bra2.contr[0].l, ket1.contr[0].l, ket2.contr[0].l);
    if (class_profiles.find(id) == class_profiles.end()) {
      class_profile dummy;
      class_profiles[id] = dummy;
    }
#endif

    // compute primitive data
#ifdef LIBINT2_ENGINE_TIMERS
    timers.start(0);
#endif
    {
      auto p = 0;
      // this is far less aggressive than should be, but proper analysis
      // involves both bra and ket *bases* and thus cannot be done on shell-set basis
      // probably ln_precision_/2 - 10 is enough
      spbra_.init(bra1, bra2, ln_precision_);
      spket_.init(ket1, ket2, ln_precision_);
      const auto npbra = spbra_.primpairs.size();
      const auto npket = spket_.primpairs.size();
      for(auto pb=0; pb!=npbra; ++pb) {
        for(auto pk=0; pk!=npket; ++pk) {
          if (spbra_.primpairs[pb].scr + spket_.primpairs[pk].scr > ln_precision_) {
            Libint_t& primdata = primdata_[p];
            const auto& sbra1 = bra1;
            const auto& sbra2 = bra2;
            const auto& sket1 = ket1;
            const auto& sket2 = ket2;
            const auto& spbra = spbra_;
            const auto& spket = spket_;
            auto pbra = pb;
            auto pket = pk;

              const auto& A = sbra1.O;
              const auto& B = sbra2.O;
              const auto& C = sket1.O;
              const auto& D = sket2.O;
              const auto& AB = spbra.AB;
              const auto& CD = spket.AB;

              const auto& spbrapp = spbra.primpairs[pbra];
              const auto& spketpp = spket.primpairs[pket];
              const auto& pbra1 = spbrapp.p1;
              const auto& pbra2 = spbrapp.p2;
              const auto& pket1 = spketpp.p1;
              const auto& pket2 = spketpp.p2;

              const auto alpha0 = sbra1.alpha[pbra1];
              const auto alpha1 = sbra2.alpha[pbra2];
              const auto alpha2 = sket1.alpha[pket1];
              const auto alpha3 = sket2.alpha[pket2];

              const auto c0 = sbra1.contr[0].coeff[pbra1];
              const auto c1 = sbra2.contr[0].coeff[pbra2];
              const auto c2 = sket1.contr[0].coeff[pket1];
              const auto c3 = sket2.contr[0].coeff[pket2];

              const auto amtot = sbra1.contr[0].l + sket1.contr[0].l +
                                 sbra2.contr[0].l + sket2.contr[0].l;

              const auto gammap = alpha0 + alpha1;
              const auto oogammap = spbrapp.one_over_gamma;
              const auto rhop = alpha0 * alpha1 * oogammap;

              const auto gammaq = alpha2 + alpha3;
              const auto oogammaq = spketpp.one_over_gamma;
              const auto rhoq = alpha2 * alpha3 * oogammaq;

              const auto& P = spbrapp.P;
              const auto& Q = spketpp.P;
              const auto PQx = P[0] - Q[0];
              const auto PQy = P[1] - Q[1];
              const auto PQz = P[2] - Q[2];
              const auto PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

              const auto K12 = spbrapp.K * spketpp.K;
              decltype(K12) two_times_M_PI_to_25(34.986836655249725693); // (2 \pi)^{5/2}
              const auto gammapq = gammap + gammaq;
              const auto sqrt_gammapq = sqrt(gammapq);
              const auto oogammapq = 1.0 / (gammapq);
              auto pfac = two_times_M_PI_to_25 * K12 * sqrt_gammapq * oogammapq;
              pfac *= c0 * c1 * c2 * c3;

              if (std::abs(pfac) >= precision_) {

              const auto rho = gammap * gammaq * oogammapq;
              const auto T = PQ2*rho;
              auto* gm_ptr = &(primdata.LIBINT_T_SS_EREP_SS(0)[0]);
              const auto mmax = amtot + deriv_order;

              if (!skip_core_ints) {
                switch (oper) {
                  case Operator::coulomb: {
                    const auto& core_eval_ptr = core_eval_pack_.as<detail::core_eval_pack_type<Operator::coulomb>>().first();
                    core_eval_ptr->eval(gm_ptr, T, mmax);
                  } break;
                  case Operator::cgtg_x_coulomb: {
                    const auto& core_eval_ptr = core_eval_pack_.as<detail::core_eval_pack_type<Operator::cgtg_x_coulomb>>().first();
                    auto& core_eval_scratch = core_eval_pack_.as<detail::core_eval_pack_type<Operator::cgtg_x_coulomb>>().second();
                    const auto& core_ints_params = core_ints_params_.as<typename operator_traits<Operator::cgtg>::oper_params_type>();
                    core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params, &core_eval_scratch);
                  } break;
                  case Operator::cgtg: {
                    const auto&  core_eval_ptr = core_eval_pack_.as<detail::core_eval_pack_type<Operator::cgtg>>().first();
                    const auto& core_ints_params = core_ints_params_.as<typename operator_traits<Operator::cgtg>::oper_params_type>();
                    core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params);
                  } break;
                  case Operator::delcgtg2: {
                    const auto& core_eval_ptr = core_eval_pack_.as<detail::core_eval_pack_type<Operator::delcgtg2>>().first();
                    const auto& core_ints_params = core_ints_params_.as<typename operator_traits<Operator::cgtg>::oper_params_type>();
                    core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params);
                  } break;
                  case Operator::delta: {
                    const auto& core_eval_ptr = core_eval_pack_.as<detail::core_eval_pack_type<Operator::delta>>().first();
                    core_eval_ptr->eval(gm_ptr, rho, T, mmax);
                  } break;
                  default:
                    assert(false); // unreachable
                }
              }

              for(auto m=0; m!=mmax+1; ++m) {
                gm_ptr[m] *= pfac;
              }

              if (mmax != 0) {

          #if LIBINT2_DEFINED(eri,PA_x)
              primdata.PA_x[0] = P[0] - A[0];
          #endif
          #if LIBINT2_DEFINED(eri,PA_y)
              primdata.PA_y[0] = P[1] - A[1];
          #endif
          #if LIBINT2_DEFINED(eri,PA_z)
              primdata.PA_z[0] = P[2] - A[2];
          #endif
          #if LIBINT2_DEFINED(eri,PB_x)
              primdata.PB_x[0] = P[0] - B[0];
          #endif
          #if LIBINT2_DEFINED(eri,PB_y)
              primdata.PB_y[0] = P[1] - B[1];
          #endif
          #if LIBINT2_DEFINED(eri,PB_z)
              primdata.PB_z[0] = P[2] - B[2];
          #endif

          #if LIBINT2_DEFINED(eri,QC_x)
              primdata.QC_x[0] = Q[0] - C[0];
          #endif
          #if LIBINT2_DEFINED(eri,QC_y)
              primdata.QC_y[0] = Q[1] - C[1];
          #endif
          #if LIBINT2_DEFINED(eri,QC_z)
              primdata.QC_z[0] = Q[2] - C[2];
          #endif
          #if LIBINT2_DEFINED(eri,QD_x)
              primdata.QD_x[0] = Q[0] - D[0];
          #endif
          #if LIBINT2_DEFINED(eri,QD_y)
              primdata.QD_y[0] = Q[1] - D[1];
          #endif
          #if LIBINT2_DEFINED(eri,QD_z)
              primdata.QD_z[0] = Q[2] - D[2];
          #endif

          #if LIBINT2_DEFINED(eri,AB_x)
              primdata.AB_x[0] = AB[0];
          #endif
          #if LIBINT2_DEFINED(eri,AB_y)
              primdata.AB_y[0] = AB[1];
          #endif
          #if LIBINT2_DEFINED(eri,AB_z)
              primdata.AB_z[0] = AB[2];
          #endif
          #if LIBINT2_DEFINED(eri,BA_x)
              primdata.BA_x[0] = -AB[0];
          #endif
          #if LIBINT2_DEFINED(eri,BA_y)
              primdata.BA_y[0] = -AB[1];
          #endif
          #if LIBINT2_DEFINED(eri,BA_z)
              primdata.BA_z[0] = -AB[2];
          #endif

          #if LIBINT2_DEFINED(eri,CD_x)
              primdata.CD_x[0] = CD[0];
          #endif
          #if LIBINT2_DEFINED(eri,CD_y)
              primdata.CD_y[0] = CD[1];
          #endif
          #if LIBINT2_DEFINED(eri,CD_z)
              primdata.CD_z[0] = CD[2];
          #endif
          #if LIBINT2_DEFINED(eri,DC_x)
              primdata.DC_x[0] = -CD[0];
          #endif
          #if LIBINT2_DEFINED(eri,DC_y)
              primdata.DC_y[0] = -CD[1];
          #endif
          #if LIBINT2_DEFINED(eri,DC_z)
              primdata.DC_z[0] = -CD[2];
          #endif

              const auto gammap_o_gammapgammaq = oogammapq * gammap;
              const auto gammaq_o_gammapgammaq = oogammapq * gammaq;

              const auto Wx = (gammap_o_gammapgammaq * P[0] + gammaq_o_gammapgammaq * Q[0]);
              const auto Wy = (gammap_o_gammapgammaq * P[1] + gammaq_o_gammapgammaq * Q[1]);
              const auto Wz = (gammap_o_gammapgammaq * P[2] + gammaq_o_gammapgammaq * Q[2]);

          #if LIBINT2_DEFINED(eri,WP_x)
              primdata.WP_x[0] = Wx - P[0];
          #endif
          #if LIBINT2_DEFINED(eri,WP_y)
              primdata.WP_y[0] = Wy - P[1];
          #endif
          #if LIBINT2_DEFINED(eri,WP_z)
              primdata.WP_z[0] = Wz - P[2];
          #endif
          #if LIBINT2_DEFINED(eri,WQ_x)
              primdata.WQ_x[0] = Wx - Q[0];
          #endif
          #if LIBINT2_DEFINED(eri,WQ_y)
              primdata.WQ_y[0] = Wy - Q[1];
          #endif
          #if LIBINT2_DEFINED(eri,WQ_z)
              primdata.WQ_z[0] = Wz - Q[2];
          #endif
          #if LIBINT2_DEFINED(eri,oo2z)
              primdata.oo2z[0] = 0.5*oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,oo2e)
              primdata.oo2e[0] = 0.5*oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,oo2ze)
              primdata.oo2ze[0] = 0.5*oogammapq;
          #endif
          #if LIBINT2_DEFINED(eri,roz)
              primdata.roz[0] = rho*oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,roe)
              primdata.roe[0] = rho*oogammaq;
          #endif

              // using ITR?
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
              primdata.TwoPRepITR_pfac0_0_0_x[0] = - (alpha1*AB[0] + alpha3*CD[0]) * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
              primdata.TwoPRepITR_pfac0_0_0_y[0] = - (alpha1*AB[1] + alpha3*CD[1]) * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
              primdata.TwoPRepITR_pfac0_0_0_z[0] = - (alpha1*AB[2] + alpha3*CD[2]) * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
              primdata.TwoPRepITR_pfac0_1_0_x[0] = - (alpha1*AB[0] + alpha3*CD[0]) * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
              primdata.TwoPRepITR_pfac0_1_0_y[0] = - (alpha1*AB[1] + alpha3*CD[1]) * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
              primdata.TwoPRepITR_pfac0_1_0_z[0] = - (alpha1*AB[2] + alpha3*CD[2]) * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
              primdata.TwoPRepITR_pfac0_0_1_x[0] = (alpha0*AB[0] + alpha2*CD[0]) * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
              primdata.TwoPRepITR_pfac0_0_1_y[0] = (alpha0*AB[1] + alpha2*CD[1]) * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
              primdata.TwoPRepITR_pfac0_0_1_z[0] = (alpha0*AB[2] + alpha2*CD[2]) * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
              primdata.TwoPRepITR_pfac0_1_1_x[0] = (alpha0*AB[0] + alpha2*CD[0]) * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
              primdata.TwoPRepITR_pfac0_1_1_y[0] = (alpha0*AB[1] + alpha2*CD[1]) * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
              primdata.TwoPRepITR_pfac0_1_1_z[0] = (alpha0*AB[2] + alpha2*CD[2]) * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,eoz)
              primdata.eoz[0] = gammaq*oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,zoe)
              primdata.zoe[0] = gammap*oogammaq;
          #endif

              // prefactors for derivative ERI relations
              if (deriv_order > 0) {
          #if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
                primdata.alpha1_rho_over_zeta2[0] = alpha0 * (oogammap * gammaq_o_gammapgammaq);
          #endif
          #if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
                primdata.alpha2_rho_over_zeta2[0] = alpha1 * (oogammap * gammaq_o_gammapgammaq);
          #endif
          #if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
                primdata.alpha3_rho_over_eta2[0] = alpha2 * (oogammaq * gammap_o_gammapgammaq);
          #endif
          #if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
                primdata.alpha4_rho_over_eta2[0] = alpha3 * (oogammaq * gammap_o_gammapgammaq);
          #endif
          #if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
                primdata.alpha1_over_zetapluseta[0] = alpha0 * oogammapq;
          #endif
          #if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
                primdata.alpha2_over_zetapluseta[0] = alpha1 * oogammapq;
          #endif
          #if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
                primdata.alpha3_over_zetapluseta[0] = alpha2 * oogammapq;
          #endif
          #if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
                primdata.alpha4_over_zetapluseta[0] = alpha3 * oogammapq;
          #endif
          #if LIBINT2_DEFINED(eri,rho12_over_alpha1)
                primdata.rho12_over_alpha1[0] = alpha1 * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,rho12_over_alpha2)
                primdata.rho12_over_alpha2[0] = alpha0 * oogammap;
          #endif
          #if LIBINT2_DEFINED(eri,rho34_over_alpha3)
                primdata.rho34_over_alpha3[0] = alpha3 * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,rho34_over_alpha4)
                primdata.rho34_over_alpha4[0] = alpha2 * oogammaq;
          #endif
          #if LIBINT2_DEFINED(eri,two_alpha0_bra)
                primdata.two_alpha0_bra[0] = 2.0 * alpha0;
          #endif
          #if LIBINT2_DEFINED(eri,two_alpha0_ket)
                primdata.two_alpha0_ket[0] = 2.0 * alpha1;
          #endif
          #if LIBINT2_DEFINED(eri,two_alpha1_bra)
                primdata.two_alpha1_bra[0] = 2.0 * alpha2;
          #endif
          #if LIBINT2_DEFINED(eri,two_alpha1_ket)
                primdata.two_alpha1_ket[0] = 2.0 * alpha3;
          #endif
              }

              } // m != 0

              ++p;
              } // prefac-based prim quartet screen

          } // rough prim quartet screen based on pair values
        } // ket prim pair
      } // bra prim pair
      primdata_[0].contrdepth = p;
    }

#ifdef LIBINT2_ENGINE_TIMERS
    const auto t0 = timers.stop(0);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
    class_profiles[id].prereqs += t0.count();
    if (primdata_[0].contrdepth != 0) {
      class_profiles[id].nshellset += 1;
      class_profiles[id].nprimset += primdata_[0].contrdepth;
    }
#  endif
#endif

    // all primitive combinations screened out? return zeroes
    if (primdata_[0].contrdepth == 0) {
      const size_t n = bra1.size() * bra2.size() * ket1.size() * ket2.size();
      memset(primdata_[0].stack, 0, sizeof(real_t)*n);
      return primdata_[0].stack;
    }

    real_t* result = nullptr;

    if (lmax == 0) { // (ss|ss)
#ifdef LIBINT2_ENGINE_TIMERS
      timers.start(1);
#endif
      auto& stack = primdata_[0].stack[0];
      stack = 0;
      for(auto p=0; p != primdata_[0].contrdepth; ++p)
        stack += primdata_[p].LIBINT_T_SS_EREP_SS(0)[0];
      primdata_[0].targets[0] = primdata_[0].stack;
#ifdef LIBINT2_ENGINE_TIMERS
      const auto t1 = timers.stop(1);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
      class_profiles[id].build_vrr += t1.count();
#  endif
#endif

      result = primdata_[0].targets[0];
    }
    else { // not (ss|ss)
#ifdef LIBINT2_ENGINE_TIMERS
#    ifdef LIBINT2_PROFILE
      const auto t1_hrr_start = primdata_[0].timers->read(0);
      const auto t1_vrr_start = primdata_[0].timers->read(1);
#    endif
      timers.start(1);
#endif
      LIBINT2_PREFIXED_NAME(libint2_build_eri)[bra1.contr[0].l][bra2.contr[0].l][ket1.contr[0].l][ket2.contr[0].l](&primdata_[0]);
#ifdef LIBINT2_ENGINE_TIMERS
      const auto t1 = timers.stop(1);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
#    ifndef LIBINT2_PROFILE
      class_profiles[id].build_vrr += t1.count();
#    else
      class_profiles[id].build_hrr += primdata_[0].timers->read(0) - t1_hrr_start;
      class_profiles[id].build_vrr += primdata_[0].timers->read(1) - t1_vrr_start;
#    endif
#  endif
#endif

      result = primdata_[0].targets[0];

#ifdef LIBINT2_ENGINE_TIMERS
      timers.start(2);
#endif

      // if needed, permute and transform
      if (use_scratch) {

        constexpr auto using_scalar_real = std::is_same<double,real_t>::value || std::is_same<float,real_t>::value;
        static_assert(using_scalar_real, "Libint2 C++11 API only supports fundamental real types");
        typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;

        // a 2-d view of the 4-d source tensor
        const auto nr1_cart = bra1.cartesian_size();
        const auto nr2_cart = bra2.cartesian_size();
        const auto nc1_cart = ket1.cartesian_size();
        const auto nc2_cart = ket2.cartesian_size();
        const auto ncol_cart = nc1_cart * nc2_cart;
        const auto nr1 = bra1.size();
        const auto nr2 = bra2.size();
        const auto nc1 = ket1.size();
        const auto nc2 = ket2.size();
        const auto nrow = nr1 * nr2;
        const auto ncol = nc1 * nc2;

        // a 2-d view of the 4-d target tensor
        const auto nr1_tgt = tbra1.size();
        const auto nr2_tgt = tbra2.size();
        const auto nc1_tgt = tket1.size();
        const auto nc2_tgt = tket2.size();
        const auto ncol_tgt = nc1_tgt * nc2_tgt;

        // transform to solid harmonics first, then unpermute, if necessary
        auto mainbuf = result;
        auto scratchbuf = &scratch_[0];
        if (bra1.contr[0].pure) {
          libint2::solidharmonics::transform_first(bra1.contr[0].l, nr2_cart*ncol_cart,
                                                   mainbuf, scratchbuf);
          std::swap(mainbuf, scratchbuf);
        }
        if (bra2.contr[0].pure) {
          libint2::solidharmonics::transform_inner(bra1.size(), bra2.contr[0].l, ncol_cart,
                                                   mainbuf, scratchbuf);
          std::swap(mainbuf, scratchbuf);
        }
        if (ket1.contr[0].pure) {
          libint2::solidharmonics::transform_inner(nrow, ket1.contr[0].l, nc2_cart,
                                                   mainbuf, scratchbuf);
          std::swap(mainbuf, scratchbuf);
        }
        if (ket2.contr[0].pure) {
          libint2::solidharmonics::transform_last(bra1.size()*bra2.size()*ket1.size(), ket2.contr[0].l,
                                                  mainbuf, scratchbuf);
          std::swap(mainbuf, scratchbuf);
        }

        // loop over rows of the source matrix
        const auto* src_row_ptr = mainbuf;
        auto tgt_ptr = scratchbuf;
        for(auto r1=0; r1!=nr1; ++r1) {
          for(auto r2=0; r2!=nr2; ++r2, src_row_ptr+=ncol) {

            typedef Eigen::Map<const Matrix> ConstMap;
            typedef Eigen::Map<Matrix> Map;
            typedef Eigen::Map<Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > StridedMap;

            // represent this source row as a matrix
            ConstMap src_blk_mat(src_row_ptr, nc1, nc2);

            // and copy to the block of the target matrix
            if (swap_braket) {
              // if swapped bra and ket, a row of source becomes a column of target
              // source row {r1,r2} is mapped to target column {r1,r2} if !swap_ket, else to {r2,r1}
              const auto tgt_col_idx = !swap_ket ? r1 * nr2 + r2 : r2 * nr1 + r1;
              StridedMap tgt_blk_mat(tgt_ptr + tgt_col_idx, nr1_tgt, nr2_tgt,
                                     Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(nr2_tgt*ncol_tgt,ncol_tgt));
              if (swap_bra)
                tgt_blk_mat = src_blk_mat.transpose();
              else
                tgt_blk_mat = src_blk_mat;
            }
            else {
              // source row {r1,r2} is mapped to target row {r1,r2} if !swap_bra, else to {r2,r1}
              const auto tgt_row_idx = !swap_bra ? r1 * nr2 + r2 : r2 * nr1 + r1;
              Map tgt_blk_mat(tgt_ptr + tgt_row_idx*ncol, nc1_tgt, nc2_tgt);
              if (swap_ket)
                tgt_blk_mat = src_blk_mat.transpose();
              else
                tgt_blk_mat = src_blk_mat;
            }

          } // end of loop
        }   // over rows of source

        result = scratchbuf;

      } // if need_scratch => needed to transpose

#ifdef LIBINT2_ENGINE_TIMERS
      const auto t2 = timers.stop(2);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
      class_profiles[id].tform += t2.count();
#  endif
#endif

    } // not (ss|ss)

    return result;
  }

#undef BOOST_PP_NBODY_OPERATOR_LIST
#undef BOOST_PP_NBODY_OPERATOR_INDEX_TUPLE
#undef BOOST_PP_NBODY_OPERATOR_INDEX_LIST
#undef BOOST_PP_NBODY_BRAKET_INDEX_TUPLE
#undef BOOST_PP_NBODY_BRAKET_INDEX_LIST
#undef BOOST_PP_NBODY_DERIV_ORDER_TUPLE
#undef BOOST_PP_NBODY_DERIV_ORDER_LIST
#undef BOOST_PP_NBODYENGINE_MCR2
#undef BOOST_PP_NBODYENGINE_MCR3
#undef BOOST_PP_NBODYENGINE_MCR4
#undef BOOST_PP_NBODYENGINE_MCR5
#undef BOOST_PP_NBODYENGINE_MCR6
#undef BOOST_PP_NBODYENGINE_MCR7

#endif // new engine

} // namespace libint2

#endif /* _libint2_src_lib_libint_engine_h_ */
