/*
 *  Copyright (C) 2004-2020 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_engineimpl_h_
#define _libint2_src_lib_libint_engineimpl_h_

#include "./engine.h"
#include "./deriv_map.h"

#include <iterator>

#pragma GCC diagnostic push
#pragma GCC system_header
#include <Eigen/Core>
#pragma GCC diagnostic pop

#include <libint2/boys.h>
#if LIBINT_HAS_SYSTEM_BOOST_PREPROCESSOR_VARIADICS
# include <boost/preprocessor.hpp>
# include <boost/preprocessor/facilities/is_1.hpp>
#else  // use bundled boost
#  include <libint2/boost/preprocessor.hpp>
#  include <libint2/boost/preprocessor/facilities/is_1.hpp>
#endif

// extra PP macros

#define BOOST_PP_MAKE_TUPLE_INTERNAL(z, i, last) \
  i BOOST_PP_COMMA_IF(BOOST_PP_NOT_EQUAL(i, last))
/// BOOST_PP_MAKE_TUPLE(n) returns (0,1,....n-1)
#define BOOST_PP_MAKE_TUPLE(n) \
  (BOOST_PP_REPEAT(n, BOOST_PP_MAKE_TUPLE_INTERNAL, BOOST_PP_DEC(n)))

// the engine will be profiled by default if library was configured with
// --enable-profile
#ifdef LIBINT2_PROFILE
#define LIBINT2_ENGINE_TIMERS
// uncomment if want to profile each integral class
#define LIBINT2_ENGINE_PROFILE_CLASS
#endif
// uncomment if want to profile the engine even if library was configured
// without --enable-profile
//#  define LIBINT2_ENGINE_TIMERS

namespace libint2 {

template <typename T, unsigned N>
typename std::remove_all_extents<T>::type* to_ptr1(T (&a)[N]) {
  return reinterpret_cast<typename std::remove_all_extents<T>::type*>(&a);
}

/// list of libint task names for each Operator type.
/// These MUST appear in the same order as in Operator.
/// You must also update BOOST_PP_NBODY_OPERATOR_LAST_ONEBODY_INDEX when you add
/// one-body ints
#define BOOST_PP_NBODY_OPERATOR_LIST \
  (overlap,                          \
   (kinetic,                         \
    (elecpot,                        \
     (elecpot,                       \
      (elecpot,                      \
       (1emultipole,                 \
        (2emultipole,                \
         (3emultipole,               \
           (sphemultipole,           \
          (eri, (eri, (eri, (eri, (eri, (eri, (eri, (eri, (eri, (eri, BOOST_PP_NIL)))))))))))))))))))

#define BOOST_PP_NBODY_OPERATOR_INDEX_TUPLE \
  BOOST_PP_MAKE_TUPLE(BOOST_PP_LIST_SIZE(BOOST_PP_NBODY_OPERATOR_LIST))
#define BOOST_PP_NBODY_OPERATOR_INDEX_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_NBODY_OPERATOR_INDEX_TUPLE)
#define BOOST_PP_NBODY_OPERATOR_LAST_ONEBODY_INDEX \
  8  // sphemultipole, the 9th member of BOOST_PP_NBODY_OPERATOR_LIST, is the last
     // 1-body operator

// make list of braket indices for n-body ints
#define BOOST_PP_NBODY_BRAKET_INDEX_TUPLE \
  BOOST_PP_MAKE_TUPLE(BOOST_PP_INC(BOOST_PP_NBODY_BRAKET_MAX_INDEX))
#define BOOST_PP_NBODY_BRAKET_INDEX_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_NBODY_BRAKET_INDEX_TUPLE)
#define BOOST_PP_NBODY_BRAKET_RANK_TUPLE (2, 3, 4)
#define BOOST_PP_NBODY_BRAKET_RANK_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_NBODY_BRAKET_RANK_TUPLE)

// make list of derivative orders for n-body ints
#define BOOST_PP_NBODY_DERIV_ORDER_TUPLE \
  BOOST_PP_MAKE_TUPLE(BOOST_PP_INC(LIBINT2_MAX_DERIV_ORDER))
#define BOOST_PP_NBODY_DERIV_ORDER_LIST \
  BOOST_PP_TUPLE_TO_LIST(BOOST_PP_NBODY_DERIV_ORDER_TUPLE)


/// the runtime version of \c operator_traits<oper>::default_params()
__libint2_engine_inline libint2::any
default_params(const Operator& oper) {
  switch (static_cast<int>(oper)) {
#define BOOST_PP_NBODYENGINE_MCR1(r, data, i, elem) \
  case i:                                           \
    return operator_traits<static_cast<Operator>(i)>::default_params();
    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_NBODYENGINE_MCR1, _,
                             BOOST_PP_NBODY_OPERATOR_LIST)
    default:
      break;
  }
  assert(false && "missing case in switch");  // unreachable
  return libint2::any();
}

/// Computes target shell sets of integrals.

/// @return vector of pointers to target shell sets, the number of sets =
/// Engine::nshellsets();
///         if the first pointer equals \c nullptr then all elements were
///         screened out.
/// \note resulting shell sets are stored in row-major order.
/// \note Call Engine::compute1() or Engine::compute2() directly to avoid extra
/// copies.
template <typename... ShellPack>
__libint2_engine_inline const Engine::target_ptr_vec& Engine::compute(
    const libint2::Shell& first_shell, const ShellPack&... rest_of_shells) {
  constexpr auto nargs = 1 + sizeof...(rest_of_shells);
  assert(nargs == braket_rank() && "# of arguments to compute() does not match the braket type");

  std::array<std::reference_wrapper<const Shell>, nargs> shells{{
      first_shell, rest_of_shells...}};

  if (operator_rank() == 1) {
    if (nargs == 2) return compute1(shells[0], shells[1]);
  } else if (operator_rank() == 2) {
    auto compute_ptr_idx = ((static_cast<int>(oper_) -
                             static_cast<int>(Operator::first_2body_oper)) *
                                nbrakets_2body +
                            (static_cast<int>(braket_) -
                             static_cast<int>(BraKet::first_2body_braket))) *
                               nderivorders_2body +
                           deriv_order_;
    auto compute_ptr = compute2_ptrs().at(compute_ptr_idx);
    assert(compute_ptr != nullptr && "2-body compute function not found");
    if (nargs == 2)
      return (this->*compute_ptr)(shells[0], Shell::unit(), shells[1],
                                  Shell::unit(), nullptr, nullptr);
    if (nargs == 3)
      return (this->*compute_ptr)(shells[0], Shell::unit(), shells[1],
                                  shells[2], nullptr, nullptr);
    if (nargs == 4)
      return (this->*compute_ptr)(shells[0], shells[1], shells[2], shells[3], nullptr, nullptr);
  }

  assert(false && "missing feature");  // only reached if missing a feature
  return targets_;
}

/// Computes target shell sets of 1-body integrals.
/// @return vector of pointers to target shell sets, the number of sets =
/// Engine::nshellsets()
/// \note resulting shell sets are stored in row-major order
__libint2_engine_inline const Engine::target_ptr_vec& Engine::compute1(
    const libint2::Shell& s1, const libint2::Shell& s2) {
  // can only handle 1 contraction at a time
  assert((s1.ncontr() == 1 && s2.ncontr() == 1) &&
         "generally-contracted shells not yet supported");

  const auto oper_is_nuclear =
      (oper_ == Operator::nuclear || oper_ == Operator::erf_nuclear ||
       oper_ == Operator::erfc_nuclear);

  const auto l1 = s1.contr[0].l;
  const auto l2 = s2.contr[0].l;
  assert(l1 <= lmax_ && "the angular momentum limit is exceeded");
  assert(l2 <= lmax_ && "the angular momentum limit is exceeded");

  // if want nuclear, make sure there is at least one nucleus .. otherwise the
  // user likely forgot to call set_params
  if (oper_is_nuclear && nparams() == 0)
    throw std::logic_error(
        "Engine<*nuclear>, but no charges found; forgot to call "
        "set_params()?");

  const auto n1 = s1.size();
  const auto n2 = s2.size();
  const auto n12 = n1 * n2;
  const auto ncart1 = s1.cartesian_size();
  const auto ncart2 = s2.cartesian_size();
  const auto ncart12 = ncart1 * ncart2;

  // assert # of primitive pairs
  const auto nprim1 = s1.nprim();
  const auto nprim2 = s2.nprim();
  const auto nprimpairs = nprim1 * nprim2;
  assert(nprimpairs <= primdata_.size() && "the max number of primitive pairs exceeded");

  auto nparam_sets = nparams();

  // keep track if need to set targets_ explicitly
  bool set_targets = set_targets_;

  // # of targets computed by libint
  const auto ntargets = nopers() * num_geometrical_derivatives(2, deriv_order_);

  // Libint computes derivatives with respect to basis functions only, must
  // must use translational invariance to recover derivatives w.r.t. operator
  // degrees of freedom
  // will compute derivs w.r.t. 2 Gaussian centers + (if nuclear) nparam_sets
  // operator centers
  const auto nderivcenters_shset = 2 + (oper_is_nuclear ? nparam_sets : 0);
  const auto nderivcoord = 3 * nderivcenters_shset;
  const auto num_shellsets_computed =
      nopers() * num_geometrical_derivatives(nderivcenters_shset, deriv_order_);

  // will use scratch_ if:
  // - Coulomb ints are computed 1 charge at a time, contributions are
  // accumulated in scratch_ (unless la==lb==0)
  // - derivatives on the missing center need to be reconstructed (no need to
  // accumulate into scratch though)
  // NB ints in scratch are packed in order
  const auto accumulate_ints_in_scratch = oper_is_nuclear;

  // adjust max angular momentum, if needed
  const auto lmax = std::max(l1, l2);
  assert(lmax <= lmax_ && "the angular momentum limit is exceeded");

  // N.B. for l=0 no need to transform to solid harmonics
  // this is a workaround for the corner case of oper_ == Operator::*nuclear,
  // and solid harmonics (s|s) integral ... beware the integral storage state
  // machine
  const auto tform_to_solids =
      (s1.contr[0].pure || s2.contr[0].pure) && lmax != 0;

  // simple (s|s) ints will be computed directly and accumulated in the first
  // element of stack
  const auto compute_directly =
      lmax == 0 && deriv_order_ == 0 &&
      (oper_ == Operator::overlap || oper_is_nuclear);
  if (compute_directly) {
    primdata_[0].stack[0] = 0;
    targets_[0] = primdata_[0].stack;
  }

  if (accumulate_ints_in_scratch)
    std::fill(std::begin(scratch_),
              std::begin(scratch_) + num_shellsets_computed * ncart12, 0.0);

  // loop over accumulation batches
  for (auto pset = 0u; pset != nparam_sets; ++pset) {
    if (!oper_is_nuclear)
      assert(nparam_sets == 1 && "unexpected number of operator parameters");

    auto p12 = 0;
    for (auto p1 = 0; p1 != nprim1; ++p1) {
      for (auto p2 = 0; p2 != nprim2; ++p2, ++p12) {
        compute_primdata(primdata_[p12], s1, s2, p1, p2, pset);
      }
    }
    primdata_[0].contrdepth = p12;

    if (compute_directly) {
      auto& result = primdata_[0].stack[0];
      switch (oper_) {
        case Operator::overlap:
          for (auto p12 = 0; p12 != primdata_[0].contrdepth; ++p12)
            result += primdata_[p12]._0_Overlap_0_x[0] *
                      primdata_[p12]._0_Overlap_0_y[0] *
                      primdata_[p12]._0_Overlap_0_z[0];
          break;
        case Operator::nuclear:
        case Operator::erf_nuclear:
        case Operator::erfc_nuclear:
          for (auto p12 = 0; p12 != primdata_[0].contrdepth; ++p12)
            result += primdata_[p12].LIBINT_T_S_ELECPOT_S(0)[0];
          break;
        default:
          assert(false && "missing case in switch");
      }
      primdata_[0].targets[0] = &result;
    } else {
      const auto buildfnidx = s1.contr[0].l * hard_lmax_ + s2.contr[0].l;
      assert(buildfnptrs_[buildfnidx] && "null build function ptr");
      buildfnptrs_[buildfnidx](&primdata_[0]);

      if (accumulate_ints_in_scratch) {
        set_targets = true;
        // - for non-derivative ints and first derivative ints the target
        //   ints computed by libint will appear at the front of targets_
        // - for second and higher derivs need to re-index targets, hence
        //   will accumulate later, when computing operator derivatives via
        //   transinv
        if (deriv_order_ <= 1) {
          // accumulate targets computed by libint for this pset into the
          // accumulated targets in scratch
          auto s_target = &scratch_[0];
          for (auto s = 0; s != ntargets; ++s, s_target += ncart12)
            if (pset != 0)
              std::transform(primdata_[0].targets[s],
                             primdata_[0].targets[s] + ncart12, s_target,
                             s_target, std::plus<value_type>());
            else
              std::copy(primdata_[0].targets[s],
                        primdata_[0].targets[s] + ncart12, s_target);
        }

        // 2. reconstruct derivatives of nuclear ints for each nucleus
        //    using translational invariance
        // NB this is done in cartesian basis, otherwise would have to tform
        // to solids contributions from every atom, rather than the running
        // total at the end
        if (deriv_order_ > 0) {
          switch (deriv_order_) {
            case 1: {
              // first 6 shellsets are derivatives with respect to Gaussian
              // positions
              // following them are derivs with respect to nuclear coordinates
              // (3 per nucleus)
              assert(ntargets == 6 && "unexpected # of targets");
              auto dest = &scratch_[0] + (6 + pset * 3) * ncart12;
              for (auto s = 0; s != 3; ++s, dest += ncart12) {
                auto src = primdata_[0].targets[s];
                for (auto i = 0; i != ncart12; ++i) {
                  dest[i] = -src[i];
                }
              }
              dest -= 3 * ncart12;
              for (auto s = 3; s != 6; ++s, dest += ncart12) {
                auto src = primdata_[0].targets[s];
                for (auto i = 0; i != ncart12; ++i) {
                  dest[i] -= src[i];
                }
              }
            } break;

            case 2: {
// computes upper triangle index
// n2 = matrix size times 2
// i,j = indices, i<j
#define upper_triangle_index_ord(n2, i, j) ((i) * ((n2) - (i)-1) / 2 + (j))
// same as above, but orders i and j
#define upper_triangle_index(n2, i, j) \
  upper_triangle_index_ord(n2, std::min((i), (j)), std::max((i), (j)))

              // accumulate ints for this pset to scratch in locations
              // remapped to overall deriv index
              const auto ncoords_times_two = nderivcoord * 2;
              for (auto d0 = 0, d01 = 0; d0 != 6; ++d0) {
                for (auto d1 = d0; d1 != 6; ++d1, ++d01) {
                  const auto d01_full =
                      upper_triangle_index_ord(ncoords_times_two, d0, d1);
                  auto tgt = &scratch_[d01_full * ncart12];
                  if (pset != 0)
                    std::transform(primdata_[0].targets[d01],
                                   primdata_[0].targets[d01] + ncart12, tgt,
                                   tgt, std::plus<value_type>());
                  else
                    std::copy(primdata_[0].targets[d01],
                              primdata_[0].targets[d01] + ncart12, tgt);
                }
              }

              // use translational invariance to build derivatives w.r.t.
              // operator centers
              {
                // mixed derivatives: first deriv w.r.t. Gaussian, second
                // w.r.t. operator coord pset
                const auto c1 = 2 + pset;
                for (auto c0 = 0; c0 != 2; ++c0) {
                  for (auto xyz0 = 0; xyz0 != 3; ++xyz0) {
                    const auto coord0 = c0 * 3 + xyz0;
                    for (auto xyz1 = 0; xyz1 != 3; ++xyz1) {
                      const auto coord1 = c1 * 3 + xyz1;  // coord1 > coord0

                      const auto coord01_abs = upper_triangle_index_ord(
                          ncoords_times_two, coord0, coord1);
                      auto tgt = &scratch_[coord01_abs * ncart12];

                      // d2 / dAi dOj = - d2 / dAi dAj
                      {
                        auto coord1_A = xyz1;
                        const auto coord01_A =
                            upper_triangle_index(12, coord0, coord1_A);
                        const auto src = primdata_[0].targets[coord01_A];
                        for (auto i = 0; i != ncart12; ++i) tgt[i] = -src[i];
                      }

                      // d2 / dAi dOj -= d2 / dAi dBj
                      {
                        auto coord1_B = 3 + xyz1;
                        const auto coord01_B =
                            upper_triangle_index(12, coord0, coord1_B);
                        const auto src = primdata_[0].targets[coord01_B];
                        for (auto i = 0; i != ncart12; ++i) tgt[i] -= src[i];
                      }
                    }
                  }
                }
              }  // mixed derivs
              {
                // operator derivs
                const auto c0 = 2 + pset;
                const auto c1 = c0;
                for (auto xyz0 = 0; xyz0 != 3; ++xyz0) {
                  const auto coord0 = c0 * 3 + xyz0;
                  for (auto xyz1 = xyz0; xyz1 != 3; ++xyz1) {
                    const auto coord1 = c1 * 3 + xyz1;  // coord1 > coord0

                    const auto coord01_abs = upper_triangle_index_ord(
                        ncoords_times_two, coord0, coord1);
                    auto tgt = &scratch_[coord01_abs * ncart12];

                    // d2 / dOi dOj = d2 / dAi dAj
                    {
                      auto coord0_A = xyz0;
                      auto coord1_A = xyz1;
                      const auto coord01_AA =
                          upper_triangle_index_ord(12, coord0_A, coord1_A);
                      const auto src = primdata_[0].targets[coord01_AA];
                      for (auto i = 0; i != ncart12; ++i) tgt[i] = src[i];
                    }

                    // d2 / dOi dOj += d2 / dAi dBj
                    {
                      auto coord0_A = xyz0;
                      auto coord1_B = 3 + xyz1;
                      const auto coord01_AB =
                          upper_triangle_index_ord(12, coord0_A, coord1_B);
                      const auto src = primdata_[0].targets[coord01_AB];
                      for (auto i = 0; i != ncart12; ++i) tgt[i] += src[i];
                    }

                    // d2 / dOi dOj += d2 / dBi dAj
                    {
                      auto coord0_B = 3 + xyz0;
                      auto coord1_A = xyz1;
                      const auto coord01_BA =
                          upper_triangle_index_ord(12, coord1_A, coord0_B);
                      const auto src = primdata_[0].targets[coord01_BA];
                      for (auto i = 0; i != ncart12; ++i) tgt[i] += src[i];
                    }

                    // d2 / dOi dOj += d2 / dBi dBj
                    {
                      auto coord0_B = 3 + xyz0;
                      auto coord1_B = 3 + xyz1;
                      const auto coord01_BB =
                          upper_triangle_index_ord(12, coord0_B, coord1_B);
                      const auto src = primdata_[0].targets[coord01_BB];
                      for (auto i = 0; i != ncart12; ++i) tgt[i] += src[i];
                    }
                  }
                }
              }  // operator derivs

#undef upper_triangle_index
            } break;

            default: {
              assert(deriv_order_ <= 2 && "feature not implemented");

              // 1. since # of derivatives changes, remap derivatives computed
              //    by libint; targets_ will hold the "remapped" pointers to
              //    the data
              using ShellSetDerivIterator =
                  libint2::FixedOrderedIntegerPartitionIterator<
                      std::vector<unsigned int>>;
              ShellSetDerivIterator shellset_gaussian_diter(deriv_order_, 2);
              ShellSetDerivIterator shellset_full_diter(deriv_order_,
                                                        nderivcenters_shset);
              std::vector<unsigned int> full_deriv(3 * nderivcenters_shset, 0);
              std::size_t s = 0;
              while (shellset_gaussian_diter) {  // loop over derivs computed
                                                 // by libint
                const auto& s1s2_deriv = *shellset_gaussian_diter;
                std::copy(std::begin(s1s2_deriv), std::end(s1s2_deriv),
                          std::begin(full_deriv));
                const auto full_rank = ShellSetDerivIterator::rank(full_deriv);
                targets_[full_rank] = primdata_[0].targets[s];
              }
              // use translational invariance to build derivatives w.r.t.
              // operator centers
            }

          }  // deriv_order_ switch
        }    // reconstruct derivatives
      }
    }  // ltot != 0
  }    // pset (accumulation batches)

  if (tform_to_solids) {
    set_targets = false;
    // where do spherical ints go?
    auto* spherical_ints =
        (accumulate_ints_in_scratch) ? scratch2_ : &scratch_[0];

    // transform to solid harmonics, one shell set at a time:
    // for each computed shell set ...
    for (auto s = 0ul; s != num_shellsets_computed;
         ++s, spherical_ints += n12) {
      auto cartesian_ints = accumulate_ints_in_scratch
                                ? &scratch_[s * ncart12]
                                : primdata_[0].targets[s];
      // transform
      if (s1.contr[0].pure && s2.contr[0].pure) {
        libint2::solidharmonics::tform(l1, l2, cartesian_ints, spherical_ints);
      } else {
        if (s1.contr[0].pure)
          libint2::solidharmonics::tform_rows(l1, n2, cartesian_ints,
                                              spherical_ints);
        else
          libint2::solidharmonics::tform_cols(n1, l2, cartesian_ints,
                                              spherical_ints);
      }
      // .. and compute the destination
      targets_[s] = spherical_ints;
    }  // loop cartesian shell set
  }    // tform to solids

  if (set_targets) {
    for (auto s = 0ul; s != num_shellsets_computed; ++s) {
      auto cartesian_ints = accumulate_ints_in_scratch
                                ? &scratch_[s * ncart12]
                                : primdata_[0].targets[s];
      targets_[s] = cartesian_ints;
    }
  }

  if (cartesian_shell_normalization() == CartesianShellNormalization::uniform) {
    std::array<std::reference_wrapper<const Shell>, 2> shells{s1,s2};
    for (auto s = 0ul; s != num_shellsets_computed; ++s) {
      uniform_normalize_cartesian_shells(const_cast<value_type*>(targets_[s]), shells);
    }
  }

  return targets_;
}

// generic _initializer
__libint2_engine_inline void Engine::_initialize() {
#define BOOST_PP_NBODYENGINE_MCR3_ncenter(product) \
  BOOST_PP_TUPLE_ELEM(3, 1, product)

#define BOOST_PP_NBODYENGINE_MCR3_default_ncenter(product)                   \
  BOOST_PP_IIF(BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(3, 0, product),          \
                                BOOST_PP_NBODY_OPERATOR_LAST_ONEBODY_INDEX), \
               4, 2)

#define BOOST_PP_NBODYENGINE_MCR3_NCENTER(product)                            \
  BOOST_PP_IIF(                                                               \
      BOOST_PP_NOT_EQUAL(BOOST_PP_NBODYENGINE_MCR3_ncenter(product),          \
                         BOOST_PP_NBODYENGINE_MCR3_default_ncenter(product)), \
      BOOST_PP_NBODYENGINE_MCR3_ncenter(product), BOOST_PP_EMPTY())

#define BOOST_PP_NBODYENGINE_MCR3_OPER(product)  \
  BOOST_PP_LIST_AT(BOOST_PP_NBODY_OPERATOR_LIST, \
                   BOOST_PP_TUPLE_ELEM(3, 0, product))

#define BOOST_PP_NBODYENGINE_MCR3_DERIV(product)                        \
  BOOST_PP_IIF(BOOST_PP_GREATER(BOOST_PP_TUPLE_ELEM(3, 2, product), 0), \
               BOOST_PP_TUPLE_ELEM(3, 2, product), BOOST_PP_EMPTY())

#define BOOST_PP_NBODYENGINE_MCR3_task(product)                         \
  BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_NBODYENGINE_MCR3_ncenter(product), \
                            BOOST_PP_NBODYENGINE_MCR3_OPER(product)),   \
               BOOST_PP_NBODYENGINE_MCR3_DERIV(product))

#define BOOST_PP_NBODYENGINE_MCR3_TASK(product)                             \
  BOOST_PP_IIF(                                                             \
      BOOST_PP_CAT(LIBINT2_TASK_EXISTS_,                                    \
                   BOOST_PP_NBODYENGINE_MCR3_task(product)),                \
      BOOST_PP_CAT(BOOST_PP_CAT(BOOST_PP_NBODYENGINE_MCR3_NCENTER(product), \
                                BOOST_PP_NBODYENGINE_MCR3_OPER(product)),   \
                   BOOST_PP_NBODYENGINE_MCR3_DERIV(product)),               \
      default)

#define BOOST_PP_NBODYENGINE_MCR3(r, product)                                  \
  if (static_cast<int>(oper_) == BOOST_PP_TUPLE_ELEM(3, 0, product) &&         \
      static_cast<int>(rank(braket_)) == BOOST_PP_TUPLE_ELEM(3, 1, product) && \
      deriv_order_ == BOOST_PP_TUPLE_ELEM(3, 2, product)) {                    \
    hard_lmax_ = BOOST_PP_CAT(LIBINT2_MAX_AM_,                                 \
                              BOOST_PP_NBODYENGINE_MCR3_TASK(product)) +       \
                 1;                                                            \
    hard_default_lmax_ =                                                       \
    BOOST_PP_IF(BOOST_PP_IS_1(BOOST_PP_CAT(LIBINT2_CENTER_DEPENDENT_MAX_AM_,   \
                              BOOST_PP_NBODYENGINE_MCR3_task(product))),       \
                              BOOST_PP_CAT(LIBINT2_MAX_AM_,                    \
                                           BOOST_PP_CAT(default,               \
                                                        BOOST_PP_NBODYENGINE_MCR3_DERIV(product) \
                                                       )                       \
                                          ) + 1, std::numeric_limits<int>::max()); \
    const auto lmax =                                                          \
    BOOST_PP_IF(BOOST_PP_IS_1(BOOST_PP_CAT(LIBINT2_CENTER_DEPENDENT_MAX_AM_,   \
                              BOOST_PP_NBODYENGINE_MCR3_task(product))),       \
      std::max(hard_lmax_,hard_default_lmax_), hard_lmax_);                    \
    if (lmax_ >= lmax) {                                                       \
      throw Engine::lmax_exceeded(                                             \
          BOOST_PP_STRINGIZE(BOOST_PP_NBODYENGINE_MCR3_TASK(product)),         \
          lmax, lmax_);                                                        \
    }                                                                          \
    if (stack_size_ > 0)                                                       \
      libint2_cleanup_default(&primdata_[0]);                                  \
    stack_size_ = LIBINT2_PREFIXED_NAME(BOOST_PP_CAT(                          \
        libint2_need_memory_, BOOST_PP_NBODYENGINE_MCR3_TASK(product)))(       \
        lmax_);                                                                \
    LIBINT2_PREFIXED_NAME(                                                     \
        BOOST_PP_CAT(libint2_init_, BOOST_PP_NBODYENGINE_MCR3_TASK(product)))  \
    (&primdata_[0], lmax_, 0);                                                 \
    BOOST_PP_IF(BOOST_PP_IS_1(LIBINT2_FLOP_COUNT),                             \
      LIBINT2_PREFIXED_NAME(libint2_init_flopcounter)                          \
    (&primdata_[0], primdata_.size()), BOOST_PP_EMPTY());                      \
    buildfnptrs_ = to_ptr1(LIBINT2_PREFIXED_NAME(BOOST_PP_CAT(                 \
        libint2_build_, BOOST_PP_NBODYENGINE_MCR3_TASK(product))));            \
    reset_scratch();                                                           \
    return;                                                                    \
  }

  BOOST_PP_LIST_FOR_EACH_PRODUCT(
      BOOST_PP_NBODYENGINE_MCR3, 3,
      (BOOST_PP_NBODY_OPERATOR_INDEX_LIST, BOOST_PP_NBODY_BRAKET_RANK_LIST,
       BOOST_PP_NBODY_DERIV_ORDER_LIST))

  assert(
      false &&
      "missing case in switch");  // either deriv_order_ or oper_ is wrong
}  // _initialize<R>()

__libint2_engine_inline void Engine::initialize(size_t max_nprim) {
  assert(libint2::initialized() && "libint is not initialized");
  assert(deriv_order_ <= LIBINT2_MAX_DERIV_ORDER &&
         "exceeded the max derivative order of the library");

  // validate braket
#ifndef INCLUDE_ONEBODY
  assert(braket_ != BraKet::x_x &&
         "this braket type not supported by the library; give --enable-1body to configure");
#endif
#ifndef INCLUDE_ERI
  assert(braket_ != BraKet::xx_xx &&
         "this braket type not supported by the library; give --enable-eri to configure");
#endif
#ifndef INCLUDE_ERI3
  assert((braket_ != BraKet::xs_xx && braket_ != BraKet::xx_xs) &&
         "this braket type not supported by the library; give --enable-eri3 to configure");
#endif
#ifndef INCLUDE_ERI2
  assert(braket_ != BraKet::xs_xs &&
         "this braket type not supported by the library; give --enable-eri2 to configure");
#endif

  // make sure it's no default initialized
  if (lmax_ < 0)
    throw using_default_initialized();

  // initialize braket, if needed
  if (braket_ == BraKet::invalid) braket_ = default_braket(oper_);

  if (max_nprim != 0) primdata_.resize(std::pow(max_nprim, braket_rank()));

  // Grab pointer to derivative index permutation map
  if (deriv_order_ > 0) {
    int ncenters;
    if (braket_ == BraKet::xx_xx) ncenters = 4; 
    if (braket_ == BraKet::xs_xx) ncenters = 3; 
    mapDerivIndex = libint2::derivmap::DerivMapGenerator::instance(deriv_order_, ncenters);
  }

  // initialize targets
  {
    decltype(targets_)::allocator_type alloc(primdata_[0].targets);
    targets_ = decltype(targets_)(alloc);
    // in some cases extra memory use can be avoided if targets_ manages its own
    // memory
    // the only instance is where we permute derivative integrals, this calls
    // for permuting
    // target indices.
    const auto permutable_targets =
        deriv_order_ > 0 &&
        (braket_ == BraKet::xx_xx || braket_ == BraKet::xs_xx ||
         braket_ == BraKet::xx_xs);
    if (permutable_targets)
      targets_.reserve(max_ntargets + 1);
    else
      targets_.reserve(max_ntargets);
    // will be resized to appropriate size in reset_scratch via _initialize
  }

#ifdef LIBINT2_ENGINE_TIMERS
  timers.set_now_overhead(25);
#endif
#ifdef LIBINT2_PROFILE
  primdata_[0].timers->set_now_overhead(25);
#endif

  _initialize();
}

namespace detail {
__libint2_engine_inline std::vector<Engine::compute2_ptr_type>
init_compute2_ptrs() {
  auto max_ncompute2_ptrs = nopers_2body * nbrakets_2body * nderivorders_2body;
  std::vector<Engine::compute2_ptr_type> result(max_ncompute2_ptrs, nullptr);

#define BOOST_PP_NBODYENGINE_MCR7(r, product)                                 \
  if (BOOST_PP_TUPLE_ELEM(3, 0, product) >=                                   \
          static_cast<int>(Operator::first_2body_oper) &&                     \
      BOOST_PP_TUPLE_ELEM(3, 0, product) <=                                   \
          static_cast<int>(Operator::last_2body_oper) &&                      \
      BOOST_PP_TUPLE_ELEM(3, 1, product) >=                                   \
          static_cast<int>(BraKet::first_2body_braket) &&                     \
      BOOST_PP_TUPLE_ELEM(3, 1, product) <=                                   \
          static_cast<int>(BraKet::last_2body_braket)) {                      \
    auto compute_ptr_idx = ((BOOST_PP_TUPLE_ELEM(3, 0, product) -             \
                             static_cast<int>(Operator::first_2body_oper)) *  \
                                nbrakets_2body +                              \
                            (BOOST_PP_TUPLE_ELEM(3, 1, product) -             \
                             static_cast<int>(BraKet::first_2body_braket))) * \
                               nderivorders_2body +                           \
                           BOOST_PP_TUPLE_ELEM(3, 2, product);                \
    result.at(compute_ptr_idx) = &Engine::compute2<                           \
        static_cast<Operator>(BOOST_PP_TUPLE_ELEM(3, 0, product)),            \
        static_cast<BraKet>(BOOST_PP_TUPLE_ELEM(3, 1, product)),              \
        BOOST_PP_TUPLE_ELEM(3, 2, product)>;                                  \
  }

  BOOST_PP_LIST_FOR_EACH_PRODUCT(
      BOOST_PP_NBODYENGINE_MCR7, 3,
      (BOOST_PP_NBODY_OPERATOR_INDEX_LIST, BOOST_PP_NBODY_BRAKET_INDEX_LIST,
       BOOST_PP_NBODY_DERIV_ORDER_LIST))

  return result;
}
}  // namespace detail

__libint2_engine_inline const std::vector<Engine::compute2_ptr_type>&
Engine::compute2_ptrs() const {
  static std::vector<compute2_ptr_type> compute2_ptrs_ =
      detail::init_compute2_ptrs();
  return compute2_ptrs_;
}

__libint2_engine_inline unsigned int Engine::nparams() const {
  switch (oper_) {
    case Operator::nuclear:
      return any_cast<const operator_traits<Operator::nuclear>::oper_params_type&>(params_)
          .size();
    case Operator::erf_nuclear:
    case Operator::erfc_nuclear:
      return std::get<1>(any_cast<const operator_traits<Operator::erfc_nuclear>::oper_params_type&>(params_))
          .size();
    default:
      return 1;
  }
  return 1;
}
__libint2_engine_inline unsigned int Engine::nopers() const {
  switch (static_cast<int>(oper_)) {
#define BOOST_PP_NBODYENGINE_MCR4(r, data, i, elem) \
  case i:                                           \
    return operator_traits<static_cast<Operator>(i)>::nopers;
    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_NBODYENGINE_MCR4, _,
                             BOOST_PP_NBODY_OPERATOR_LIST)
    default:
      break;
  }
  assert(false && "missing case in switch");  // unreachable
  return 0;
}

template <>
__libint2_engine_inline any Engine::enforce_params_type<any>(
    Operator oper, const any& params, bool throw_if_wrong_type) {
  any result;
  switch (static_cast<int>(oper)) {
#define BOOST_PP_NBODYENGINE_MCR5A(r, data, i, elem)                           \
  case i:                                                                      \
    if (any_cast<operator_traits<static_cast<Operator>(i)>::oper_params_type>( \
            &params) != nullptr) {                                             \
      result = params;                                                         \
    } else {                                                                   \
      if (throw_if_wrong_type) throw bad_any_cast();                           \
      result = operator_traits<static_cast<Operator>(i)>::default_params();    \
    }                                                                          \
    break;

    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_NBODYENGINE_MCR5A, _,
                             BOOST_PP_NBODY_OPERATOR_LIST)

    default:
      assert(false && "missing case in switch");  // missed a case?
  }
  return result;
}

template <typename Params>
__libint2_engine_inline any Engine::enforce_params_type(
    Operator oper, const Params& params, bool throw_if_wrong_type) {
  any result;
  switch (static_cast<int>(oper)) {
#define BOOST_PP_NBODYENGINE_MCR5B(r, data, i, elem)                         \
  case i:                                                                   \
    if (std::is_same<Params, operator_traits<static_cast<Operator>(         \
                                 i)>::oper_params_type>::value) {           \
      result = params;                                                      \
    } else {                                                                \
      if (throw_if_wrong_type) throw std::bad_cast();                       \
      result = operator_traits<static_cast<Operator>(i)>::default_params(); \
    }                                                                       \
    break;

    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_NBODYENGINE_MCR5B, _,
                             BOOST_PP_NBODY_OPERATOR_LIST)

    default:
      assert(false && "missing case in switch");  // missed a case?
  }
  return result;
}

__libint2_engine_inline any Engine::make_core_eval_pack(Operator oper) const {
  any result;
  switch (static_cast<int>(oper)) {
#define BOOST_PP_NBODYENGINE_MCR6(r, data, i, elem)                          \
  case i:                                                                    \
    result = libint2::detail::make_compressed_pair(                          \
        operator_traits<static_cast<Operator>(i)>::core_eval_type::instance( \
            braket_rank() * lmax_ + deriv_order_,                            \
            std::numeric_limits<scalar_type>::epsilon()),                    \
        libint2::detail::CoreEvalScratch<                                    \
            operator_traits<static_cast<Operator>(i)>::core_eval_type>(      \
            braket_rank() * lmax_ + deriv_order_));                          \
    assert(any_cast<detail::core_eval_pack_type<static_cast<Operator>(i)>>(  \
               &result) != nullptr);                                         \
    break;

    BOOST_PP_LIST_FOR_EACH_I(BOOST_PP_NBODYENGINE_MCR6, _,
                             BOOST_PP_NBODY_OPERATOR_LIST)

    default:
      assert(false && "missing case in switch");  // missed a case?
  }
  return result;
}

__libint2_engine_inline void Engine::init_core_ints_params(const any& params) {
  if (oper_ == Operator::delcgtg2) {
    // [g12,[- \Del^2, g12] = 2 (\Del g12) \cdot (\Del g12)
    // (\Del exp(-a r_12^2) \cdot (\Del exp(-b r_12^2) = 4 a b (r_{12}^2 exp(-
    // (a+b) r_{12}^2) )
    // i.e. need to scale each coefficient by 4 a b
    auto oparams =
        any_cast<const operator_traits<Operator::delcgtg2>::oper_params_type&>(params);
    const auto ng = oparams.size();
    operator_traits<Operator::delcgtg2>::oper_params_type core_ints_params;
    core_ints_params.reserve(ng * (ng + 1) / 2);
    for (size_t b = 0; b < ng; ++b)
      for (size_t k = 0; k <= b; ++k) {
        const auto gexp = oparams[b].first + oparams[k].first;
        const auto gcoeff = oparams[b].second * oparams[k].second *
                            (b == k ? 1 : 2);  // if a != b include ab and ba
        const auto gcoeff_rescaled =
            4 * oparams[b].first * oparams[k].first * gcoeff;
        core_ints_params.push_back(std::make_pair(gexp, gcoeff_rescaled));
      }
    core_ints_params_ = core_ints_params;
  } else {
    core_ints_params_ = params;
  }
}

__libint2_engine_inline void Engine::compute_primdata(Libint_t& primdata, const Shell& s1,
                                     const Shell& s2, size_t p1, size_t p2,
                                     size_t oset) {
  const auto& A = s1.O;
  const auto& B = s2.O;

  const auto alpha1 = s1.alpha[p1];
  const auto alpha2 = s2.alpha[p2];

  const auto c1 = s1.contr[0].coeff[p1];
  const auto c2 = s2.contr[0].coeff[p2];

  const auto gammap = alpha1 + alpha2;
  const auto oogammap = 1 / gammap;
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

  assert(LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD && "non-standard shell ordering");

  const auto oper_is_nuclear =
      (oper_ == Operator::nuclear || oper_ == Operator::erf_nuclear ||
       oper_ == Operator::erfc_nuclear);

  // need to use HRR? see strategy.cc
  const auto l1 = s1.contr[0].l;
  const auto l2 = s2.contr[0].l;
  const bool use_hrr = (oper_is_nuclear || oper_ == Operator::sphemultipole) && l1 > 0 && l2 > 0;
  // unlike the 2-body ints, can go both ways, determine which way to go (the logic must match TwoCenter_OS_Tactic)
  const bool hrr_ket_to_bra = l1 >= l2;
  if (use_hrr) {
    if (hrr_ket_to_bra) {
#if LIBINT2_DEFINED(eri, AB_x)
    primdata.AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri, AB_y)
    primdata.AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri, AB_z)
    primdata.AB_z[0] = AB_z;
#endif
    }
    else {
#if LIBINT2_DEFINED(eri, BA_x)
    primdata.BA_x[0] = - AB_x;
#endif
#if LIBINT2_DEFINED(eri, BA_y)
    primdata.BA_y[0] = - AB_y;
#endif
#if LIBINT2_DEFINED(eri, BA_z)
    primdata.BA_z[0] = - AB_z;
#endif
    }
  }

  // figure out whether will do VRR on center A and/or B
//  if ((!use_hrr && l1 > 0) || hrr_ket_to_bra) {
#if LIBINT2_DEFINED(eri, PA_x)
  primdata.PA_x[0] = Px - A[0];
#endif
#if LIBINT2_DEFINED(eri, PA_y)
  primdata.PA_y[0] = Py - A[1];
#endif
#if LIBINT2_DEFINED(eri, PA_z)
  primdata.PA_z[0] = Pz - A[2];
#endif
//  }
//
//  if ((!use_hrr && l2 > 0) || !hrr_ket_to_bra) {
#if LIBINT2_DEFINED(eri, PB_x)
    primdata.PB_x[0] = Px - B[0];
#endif
#if LIBINT2_DEFINED(eri, PB_y)
    primdata.PB_y[0] = Py - B[1];
#endif
#if LIBINT2_DEFINED(eri, PB_z)
    primdata.PB_z[0] = Pz - B[2];
#endif
//  }

  if (oper_ == Operator::emultipole1 || oper_ == Operator::emultipole2 ||
      oper_ == Operator::emultipole3) {
    const auto& O = any_cast<const operator_traits<
        Operator::emultipole1>::oper_params_type&>(params_);  // same as emultipoleX
#if LIBINT2_DEFINED(eri, BO_x)
    primdata.BO_x[0] = B[0] - O[0];
#endif
#if LIBINT2_DEFINED(eri, BO_y)
    primdata.BO_y[0] = B[1] - O[1];
#endif
#if LIBINT2_DEFINED(eri, BO_z)
    primdata.BO_z[0] = B[2] - O[2];
#endif
  }
  if (oper_ == Operator::sphemultipole) {
    const auto& O = any_cast<const operator_traits<
        Operator::emultipole1>::oper_params_type&>(params_);
#if LIBINT2_DEFINED(eri, PO_x)
    primdata.PO_x[0] = Px - O[0];
#endif
#if LIBINT2_DEFINED(eri, PO_y)
    primdata.PO_y[0] = Py - O[1];
#endif
#if LIBINT2_DEFINED(eri, PO_z)
    primdata.PO_z[0] = Pz - O[2];
#endif
#if LIBINT2_DEFINED(eri, PO2)
    primdata.PO2[0] = (Px - O[0])*(Px - O[0]) + (Py - O[1])*(Py - O[1]) + (Pz - O[2])*(Pz - O[2]);
#endif
  }

#if LIBINT2_DEFINED(eri, oo2z)
  primdata.oo2z[0] = 0.5 * oogammap;
#endif

  decltype(c1) sqrt_PI(1.77245385090551602729816748334);
  const auto xyz_pfac = sqrt_PI * sqrt(oogammap);
  const auto ovlp_ss_x = exp(-rhop * AB2_x) * xyz_pfac * c1 * c2 * scale_;
  const auto ovlp_ss_y = exp(-rhop * AB2_y) * xyz_pfac;
  const auto ovlp_ss_z = exp(-rhop * AB2_z) * xyz_pfac;

  primdata._0_Overlap_0_x[0] = ovlp_ss_x;
  primdata._0_Overlap_0_y[0] = ovlp_ss_y;
  primdata._0_Overlap_0_z[0] = ovlp_ss_z;

  if (oper_ == Operator::kinetic || (deriv_order_ > 0)) {
#if LIBINT2_DEFINED(eri, two_alpha0_bra)
    primdata.two_alpha0_bra[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_ket)
    primdata.two_alpha0_ket[0] = 2.0 * alpha2;
#endif
  }

  if (oper_is_nuclear) {

    const auto& params = (oper_ == Operator::nuclear) ?
        any_cast<const operator_traits<Operator::nuclear>::oper_params_type&>(params_) :
        std::get<1>(any_cast<const operator_traits<Operator::erfc_nuclear>::oper_params_type&>(params_));

    const auto& C = params[oset].second;
    const auto& q = params[oset].first;
#if LIBINT2_DEFINED(eri, PC_x)
    primdata.PC_x[0] = Px - C[0];
#endif
#if LIBINT2_DEFINED(eri, PC_y)
    primdata.PC_y[0] = Py - C[1];
#endif
#if LIBINT2_DEFINED(eri, PC_z)
    primdata.PC_z[0] = Pz - C[2];
#endif

#if LIBINT2_DEFINED(eri, rho12_over_alpha1) || \
    LIBINT2_DEFINED(eri, rho12_over_alpha2)
    if (deriv_order_ > 0) {
#if LIBINT2_DEFINED(eri, rho12_over_alpha1)
      primdata.rho12_over_alpha1[0] = rhop_over_alpha1;
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
      primdata.rho12_over_alpha2[0] = alpha1 * oogammap;
#endif
    }
#endif
#if LIBINT2_DEFINED(eri, PC_x) && LIBINT2_DEFINED(eri, PC_y) && \
    LIBINT2_DEFINED(eri, PC_z)
    const auto PC2 = primdata.PC_x[0] * primdata.PC_x[0] +
                     primdata.PC_y[0] * primdata.PC_y[0] +
                     primdata.PC_z[0] * primdata.PC_z[0];
    const scalar_type U = gammap * PC2;
    const scalar_type rho = rhop;
    const auto mmax = s1.contr[0].l + s2.contr[0].l + deriv_order_;
    auto* fm_ptr = &(primdata.LIBINT_T_S_ELECPOT_S(0)[0]);
    if (oper_ == Operator::nuclear) {
      auto fm_engine_ptr =
          any_cast<const detail::core_eval_pack_type<Operator::nuclear>&>(core_eval_pack_)
          .first();
      fm_engine_ptr->eval(fm_ptr, U, mmax);
    } else if (oper_ == Operator::erf_nuclear) {
      const auto& core_eval_ptr =
          any_cast<const detail::core_eval_pack_type<Operator::erf_nuclear>&>(core_eval_pack_)
            .first();
      auto core_ints_params =
          std::get<0>(any_cast<const typename operator_traits<
            Operator::erf_nuclear>::oper_params_type&>(core_ints_params_));
      core_eval_ptr->eval(fm_ptr, rho, U, mmax, core_ints_params);
    } else if (oper_ == Operator::erfc_nuclear) {
      const auto& core_eval_ptr =
          any_cast<const detail::core_eval_pack_type<Operator::erfc_nuclear>&>(core_eval_pack_)
            .first();
      auto core_ints_params =
          std::get<0>(any_cast<const typename operator_traits<
            Operator::erfc_nuclear>::oper_params_type&>(core_ints_params_));
      core_eval_ptr->eval(fm_ptr, rho, U, mmax, core_ints_params);
    }

    decltype(U) two_o_sqrt_PI(1.12837916709551257389615890312);
    const auto pfac =
        -q * sqrt(gammap) * two_o_sqrt_PI * ovlp_ss_x * ovlp_ss_y * ovlp_ss_z;
    const auto m_fence = mmax + 1;
    for (auto m = 0; m != m_fence; ++m) {
      fm_ptr[m] *= pfac;
    }
#endif
  }
}  // Engine::compute_primdata()

/// computes shell set of integrals of 2-body operator
/// \note result is stored in the "chemists"/Mulliken form, (tbra1 tbra2 |tket1
/// tket2), i.e. bra and ket are in chemists meaning; result is packed in
/// row-major order.
template <Operator op, BraKet bk, size_t deriv_order>
__libint2_engine_inline const Engine::target_ptr_vec& Engine::compute2(
    const libint2::Shell& tbra1, const libint2::Shell& tbra2,
    const libint2::Shell& tket1, const libint2::Shell& tket2,
    const ShellPair* tspbra, const ShellPair* tspket) {
  assert(op == oper_ && "Engine::compute2 -- operator mismatch");
  assert(bk == braket_ && "Engine::compute2 -- braket mismatch");
  assert(deriv_order == deriv_order_ &&
         "Engine::compute2 -- deriv_order mismatch");
  assert(((tspbra == nullptr && tspket == nullptr) || (tspbra != nullptr && tspket != nullptr)) &&
         "Engine::compute2 -- expects zero or two ShellPair objects");

  //
  // i.e. bra and ket refer to chemists bra and ket
  //

  // can only handle 1 contraction at a time
  assert((tbra1.ncontr() == 1 && tbra2.ncontr() == 1 && tket1.ncontr() == 1 &&
          tket2.ncontr() == 1) && "generally-contracted shells are not yet supported");

  // angular momentum limit obeyed?
  assert(tbra1.contr[0].l <= lmax_ && "the angular momentum limit is exceeded");
  assert(tbra2.contr[0].l <= lmax_ && "the angular momentum limit is exceeded");
  assert(tket1.contr[0].l <= lmax_ && "the angular momentum limit is exceeded");
  assert(tket2.contr[0].l <= lmax_ && "the angular momentum limit is exceeded");

#if LIBINT2_SHELLQUARTET_SET == \
    LIBINT2_SHELLQUARTET_SET_STANDARD  // standard angular momentum ordering
  const auto swap_tbra = (tbra1.contr[0].l < tbra2.contr[0].l);
  const auto swap_tket = (tket1.contr[0].l < tket2.contr[0].l);
  const auto swap_braket =
      ((braket_ == BraKet::xx_xx) && (tbra1.contr[0].l + tbra2.contr[0].l >
                                     tket1.contr[0].l + tket2.contr[0].l)) ||
        braket_ == BraKet::xx_xs;
#else  // orca angular momentum ordering
  const auto swap_tbra = (tbra1.contr[0].l > tbra2.contr[0].l);
  const auto swap_tket = (tket1.contr[0].l > tket2.contr[0].l);
  const auto swap_braket =
      ((braket_ == BraKet::xx_xx) && (tbra1.contr[0].l + tbra2.contr[0].l <
                                     tket1.contr[0].l + tket2.contr[0].l)) ||
        braket_ == BraKet::xx_xs;
  assert(false && "feature not implemented");
#endif
  const auto& bra1 =
      swap_braket ? (swap_tket ? tket2 : tket1) : (swap_tbra ? tbra2 : tbra1);
  const auto& bra2 =
      swap_braket ? (swap_tket ? tket1 : tket2) : (swap_tbra ? tbra1 : tbra2);
  const auto& ket1 =
      swap_braket ? (swap_tbra ? tbra2 : tbra1) : (swap_tket ? tket2 : tket1);
  const auto& ket2 =
      swap_braket ? (swap_tbra ? tbra1 : tbra2) : (swap_tket ? tket1 : tket2);
  const auto swap_bra = swap_braket ? swap_tket : swap_tbra;
  const auto swap_ket = swap_braket ? swap_tbra : swap_tket;
  // "permute" also the user-provided shell pair data
  const auto* spbra_precomputed = swap_braket ? tspket : tspbra;
  const auto* spket_precomputed = swap_braket ? tspbra : tspket;

  const auto tform = bra1.contr[0].pure || bra2.contr[0].pure ||
                     ket1.contr[0].pure || ket2.contr[0].pure;
  const auto permute = swap_braket || swap_tbra || swap_tket;
  const auto use_scratch = permute || tform;

  // assert # of primitive pairs
  auto nprim_bra1 = bra1.nprim();
  auto nprim_bra2 = bra2.nprim();
  auto nprim_ket1 = ket1.nprim();
  auto nprim_ket2 = ket2.nprim();

  // adjust max angular momentum, if needed
  auto lmax = std::max(std::max(bra1.contr[0].l, bra2.contr[0].l),
                       std::max(ket1.contr[0].l, ket2.contr[0].l));
  assert(lmax <= lmax_ && "the angular momentum limit is exceeded");
  const auto lmax_bra = std::max(bra1.contr[0].l, bra2.contr[0].l);
  const auto lmax_ket = std::max(ket1.contr[0].l, ket2.contr[0].l);

#ifdef LIBINT2_ENGINE_PROFILE_CLASS
  class_id id(bra1.contr[0].l, bra2.contr[0].l, ket1.contr[0].l,
              ket2.contr[0].l);
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
    // initialize shell pairs, if not given ...
    // using ln_precision_ is far less aggressive than should be, but proper analysis
    // involves both bra and ket *bases* and thus cannot be done on shell-set
    // basis ... probably ln_precision_/2 - 10 is enough
    const ShellPair& spbra = spbra_precomputed ? *spbra_precomputed : (spbra_.init(bra1, bra2, ln_precision_), spbra_) ;
    const ShellPair& spket = spket_precomputed ? *spket_precomputed : (spket_.init(ket1, ket2, ln_precision_), spket_);
    // determine whether shell pair data refers to the actual ({bra1,bra2}) or swapped ({bra2,bra1}) pairs
    // if computed the shell pair data here then it's always in actual order, otherwise check swap_bra/swap_ket
    const auto spbra_is_swapped = spbra_precomputed ? swap_bra : false;
    const auto spket_is_swapped = spket_precomputed ? swap_ket : false;

    using real_t = Shell::real_t;
    // swapping bra turns AB into BA = -AB
    real_t BA[3];
    if (spbra_is_swapped) {
      for(auto xyz=0; xyz!=3; ++xyz)
        BA[xyz] = - spbra_precomputed->AB[xyz];
    }
    const auto& AB = spbra_is_swapped ? BA : spbra.AB;
    // swapping ket turns CD into DC = -CD
    real_t DC[3];
    if (spket_is_swapped) {
      for(auto xyz=0; xyz!=3; ++xyz)
        DC[xyz] = - spket_precomputed->AB[xyz];
    }
    const auto& CD = spket_is_swapped ? DC : spket.AB;

    const auto& A = bra1.O;
    const auto& B = bra2.O;
    const auto& C = ket1.O;
    const auto& D = ket2.O;

    // compute all primitive quartet data
    const auto npbra = spbra.primpairs.size();
    const auto npket = spket.primpairs.size();
    for (auto pb = 0; pb != npbra; ++pb) {
      for (auto pk = 0; pk != npket; ++pk) {
        // primitive quartet screening
        if (spbra.primpairs[pb].scr + spket.primpairs[pk].scr >
            ln_precision_) {
          Libint_t& primdata = primdata_[p];
          const auto& sbra1 = bra1;
          const auto& sbra2 = bra2;
          const auto& sket1 = ket1;
          const auto& sket2 = ket2;
          auto pbra = pb;
          auto pket = pk;

          const auto& spbrapp = spbra.primpairs[pbra];
          const auto& spketpp = spket.primpairs[pket];
          // if shell-pair data given by user
          const auto& pbra1 = spbra_is_swapped ? spbrapp.p2 : spbrapp.p1;
          const auto& pbra2 = spbra_is_swapped ? spbrapp.p1 : spbrapp.p2;
          const auto& pket1 = spket_is_swapped ? spketpp.p2 : spketpp.p1;
          const auto& pket2 = spket_is_swapped ? spketpp.p1 : spketpp.p2;

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
          decltype(K12) two_times_M_PI_to_25(
              34.986836655249725693);  // (2 \pi)^{5/2}
          const auto gammapq = gammap + gammaq;
          const auto sqrt_gammapq = sqrt(gammapq);
          const auto oogammapq = 1.0 / (gammapq);
          auto pfac = two_times_M_PI_to_25 * K12 * sqrt_gammapq * oogammapq;
          pfac *= c0 * c1 * c2 * c3 * scale_;

          if (std::abs(pfac) >= precision_) {
            const scalar_type rho = gammap * gammaq * oogammapq;
            const scalar_type T = PQ2 * rho;
            auto* gm_ptr = &(primdata.LIBINT_T_SS_EREP_SS(0)[0]);
            const auto mmax = amtot + deriv_order;

            if (!skip_core_ints) {
              switch (oper_) {
                case Operator::coulomb: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::coulomb>&>(core_eval_pack_)
                          .first();
                  core_eval_ptr->eval(gm_ptr, T, mmax);
                } break;
                case Operator::cgtg_x_coulomb: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<
                              Operator::cgtg_x_coulomb>&>(core_eval_pack_)
                          .first();
                  auto& core_eval_scratch = any_cast<detail::core_eval_pack_type<
                                                    Operator::cgtg_x_coulomb>&>(core_eval_pack_)
                                                .second();
                  const auto& core_ints_params =
                      any_cast<const typename operator_traits<
                      Operator::cgtg>::oper_params_type&>(core_ints_params_);
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params,
                                      &core_eval_scratch);
                } break;
                case Operator::cgtg: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::cgtg>&>(core_eval_pack_)
                          .first();
                  const auto& core_ints_params =
                      any_cast<const typename operator_traits<
                          Operator::cgtg>::oper_params_type&>(core_ints_params_);
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params);
                } break;
                case Operator::delcgtg2: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::delcgtg2>&>(core_eval_pack_)
                          .first();
                  const auto& core_ints_params =
                      any_cast<const typename operator_traits<
                          Operator::cgtg>::oper_params_type&>(core_ints_params_);
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params);
                } break;
                case Operator::delta: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::delta>&>(core_eval_pack_)
                          .first();
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax);
                } break;
                case Operator::r12: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::r12>&>(core_eval_pack_)
                          .first();
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax);
                } break;
                case Operator::erf_coulomb: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::erf_coulomb>&>(core_eval_pack_)
                          .first();
                  auto core_ints_params =
                      any_cast<const typename operator_traits<
                          Operator::erf_coulomb>::oper_params_type&>(core_ints_params_);
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params);
                } break;
                case Operator::erfc_coulomb: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::erfc_coulomb>&>(core_eval_pack_)
                          .first();
                  auto core_ints_params =
                      any_cast<const typename operator_traits<
                          Operator::erfc_coulomb>::oper_params_type&>(core_ints_params_);
                  core_eval_ptr->eval(gm_ptr, rho, T, mmax, core_ints_params);
                } break;
                case Operator::stg: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::stg>&>(core_eval_pack_)
                          .first();
                  auto zeta =
                      any_cast<const typename operator_traits<
                          Operator::stg>::oper_params_type&>(core_ints_params_);
                  const auto one_over_rho = gammapq * oogammap * oogammaq;
                  core_eval_ptr->eval_slater(gm_ptr, one_over_rho, T, mmax, zeta);
                } break;
                case Operator::stg_x_coulomb: {
                  const auto& core_eval_ptr =
                      any_cast<const detail::core_eval_pack_type<Operator::stg>&>(core_eval_pack_)
                          .first();
                  auto zeta =
                      any_cast<const typename operator_traits<
                          Operator::stg>::oper_params_type&>(core_ints_params_);
                  const auto one_over_rho = gammapq * oogammap * oogammaq;
                  core_eval_ptr->eval_yukawa(gm_ptr, one_over_rho, T, mmax, zeta);
                } break;
                default:
                  assert(false && "missing case in a switch");  // unreachable
              }
            }

            for (auto m = 0; m != mmax + 1; ++m) {
              gm_ptr[m] *= pfac;
            }

            if (mmax != 0) {
              if (braket_ == BraKet::xx_xx) {
#if LIBINT2_DEFINED(eri, PA_x)
                primdata.PA_x[0] = P[0] - A[0];
#endif
#if LIBINT2_DEFINED(eri, PA_y)
                primdata.PA_y[0] = P[1] - A[1];
#endif
#if LIBINT2_DEFINED(eri, PA_z)
                primdata.PA_z[0] = P[2] - A[2];
#endif
#if LIBINT2_DEFINED(eri, PB_x)
                primdata.PB_x[0] = P[0] - B[0];
#endif
#if LIBINT2_DEFINED(eri, PB_y)
                primdata.PB_y[0] = P[1] - B[1];
#endif
#if LIBINT2_DEFINED(eri, PB_z)
                primdata.PB_z[0] = P[2] - B[2];
#endif
              }

              if (braket_ != BraKet::xs_xs) {
#if LIBINT2_DEFINED(eri, QC_x)
                primdata.QC_x[0] = Q[0] - C[0];
#endif
#if LIBINT2_DEFINED(eri, QC_y)
                primdata.QC_y[0] = Q[1] - C[1];
#endif
#if LIBINT2_DEFINED(eri, QC_z)
                primdata.QC_z[0] = Q[2] - C[2];
#endif
#if LIBINT2_DEFINED(eri, QD_x)
                primdata.QD_x[0] = Q[0] - D[0];
#endif
#if LIBINT2_DEFINED(eri, QD_y)
                primdata.QD_y[0] = Q[1] - D[1];
#endif
#if LIBINT2_DEFINED(eri, QD_z)
                primdata.QD_z[0] = Q[2] - D[2];
#endif
              }

              if (braket_ == BraKet::xx_xx) {
#if LIBINT2_DEFINED(eri, AB_x)
                primdata.AB_x[0] = AB[0];
#endif
#if LIBINT2_DEFINED(eri, AB_y)
                primdata.AB_y[0] = AB[1];
#endif
#if LIBINT2_DEFINED(eri, AB_z)
                primdata.AB_z[0] = AB[2];
#endif
#if LIBINT2_DEFINED(eri, BA_x)
                primdata.BA_x[0] = -AB[0];
#endif
#if LIBINT2_DEFINED(eri, BA_y)
                primdata.BA_y[0] = -AB[1];
#endif
#if LIBINT2_DEFINED(eri, BA_z)
                primdata.BA_z[0] = -AB[2];
#endif
              }

              if (braket_ != BraKet::xs_xs) {
#if LIBINT2_DEFINED(eri, CD_x)
                primdata.CD_x[0] = CD[0];
#endif
#if LIBINT2_DEFINED(eri, CD_y)
                primdata.CD_y[0] = CD[1];
#endif
#if LIBINT2_DEFINED(eri, CD_z)
                primdata.CD_z[0] = CD[2];
#endif
#if LIBINT2_DEFINED(eri, DC_x)
                primdata.DC_x[0] = -CD[0];
#endif
#if LIBINT2_DEFINED(eri, DC_y)
                primdata.DC_y[0] = -CD[1];
#endif
#if LIBINT2_DEFINED(eri, DC_z)
                primdata.DC_z[0] = -CD[2];
#endif
              }

              const auto gammap_o_gammapgammaq = oogammapq * gammap;
              const auto gammaq_o_gammapgammaq = oogammapq * gammaq;

              const auto Wx =
                  (gammap_o_gammapgammaq * P[0] + gammaq_o_gammapgammaq * Q[0]);
              const auto Wy =
                  (gammap_o_gammapgammaq * P[1] + gammaq_o_gammapgammaq * Q[1]);
              const auto Wz =
                  (gammap_o_gammapgammaq * P[2] + gammaq_o_gammapgammaq * Q[2]);

              if (deriv_order > 0 || lmax_bra > 0) {
#if LIBINT2_DEFINED(eri, WP_x)
                primdata.WP_x[0] = Wx - P[0];
#endif
#if LIBINT2_DEFINED(eri, WP_y)
                primdata.WP_y[0] = Wy - P[1];
#endif
#if LIBINT2_DEFINED(eri, WP_z)
                primdata.WP_z[0] = Wz - P[2];
#endif
              }
              if (deriv_order > 0 || lmax_ket > 0) {
#if LIBINT2_DEFINED(eri, WQ_x)
                primdata.WQ_x[0] = Wx - Q[0];
#endif
#if LIBINT2_DEFINED(eri, WQ_y)
                primdata.WQ_y[0] = Wy - Q[1];
#endif
#if LIBINT2_DEFINED(eri, WQ_z)
                primdata.WQ_z[0] = Wz - Q[2];
#endif
              }
#if LIBINT2_DEFINED(eri, oo2z)
              primdata.oo2z[0] = 0.5 * oogammap;
#endif
#if LIBINT2_DEFINED(eri, oo2e)
              primdata.oo2e[0] = 0.5 * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, oo2ze)
              primdata.oo2ze[0] = 0.5 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri, roz)
              primdata.roz[0] = rho * oogammap;
#endif
#if LIBINT2_DEFINED(eri, roe)
              primdata.roe[0] = rho * oogammaq;
#endif

// using ITR?
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_x)
              primdata.TwoPRepITR_pfac0_0_0_x[0] =
                  -(alpha1 * AB[0] + alpha3 * CD[0]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_y)
              primdata.TwoPRepITR_pfac0_0_0_y[0] =
                  -(alpha1 * AB[1] + alpha3 * CD[1]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_0_z)
              primdata.TwoPRepITR_pfac0_0_0_z[0] =
                  -(alpha1 * AB[2] + alpha3 * CD[2]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_x)
              primdata.TwoPRepITR_pfac0_1_0_x[0] =
                  -(alpha1 * AB[0] + alpha3 * CD[0]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_y)
              primdata.TwoPRepITR_pfac0_1_0_y[0] =
                  -(alpha1 * AB[1] + alpha3 * CD[1]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_0_z)
              primdata.TwoPRepITR_pfac0_1_0_z[0] =
                  -(alpha1 * AB[2] + alpha3 * CD[2]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_x)
              primdata.TwoPRepITR_pfac0_0_1_x[0] =
                  (alpha0 * AB[0] + alpha2 * CD[0]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_y)
              primdata.TwoPRepITR_pfac0_0_1_y[0] =
                  (alpha0 * AB[1] + alpha2 * CD[1]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_0_1_z)
              primdata.TwoPRepITR_pfac0_0_1_z[0] =
                  (alpha0 * AB[2] + alpha2 * CD[2]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_x)
              primdata.TwoPRepITR_pfac0_1_1_x[0] =
                  (alpha0 * AB[0] + alpha2 * CD[0]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_y)
              primdata.TwoPRepITR_pfac0_1_1_y[0] =
                  (alpha0 * AB[1] + alpha2 * CD[1]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, TwoPRepITR_pfac0_1_1_z)
              primdata.TwoPRepITR_pfac0_1_1_z[0] =
                  (alpha0 * AB[2] + alpha2 * CD[2]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, eoz)
              primdata.eoz[0] = gammaq * oogammap;
#endif
#if LIBINT2_DEFINED(eri, zoe)
              primdata.zoe[0] = gammap * oogammaq;
#endif

              // prefactors for derivative ERI relations
              if (deriv_order > 0) {
#if LIBINT2_DEFINED(eri, alpha1_rho_over_zeta2)
                primdata.alpha1_rho_over_zeta2[0] =
                    alpha0 * (oogammap * gammaq_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha2_rho_over_zeta2)
                primdata.alpha2_rho_over_zeta2[0] =
                    alpha1 * (oogammap * gammaq_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha3_rho_over_eta2)
                primdata.alpha3_rho_over_eta2[0] =
                    alpha2 * (oogammaq * gammap_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha4_rho_over_eta2)
                primdata.alpha4_rho_over_eta2[0] =
                    alpha3 * (oogammaq * gammap_o_gammapgammaq);
#endif
#if LIBINT2_DEFINED(eri, alpha1_over_zetapluseta)
                primdata.alpha1_over_zetapluseta[0] = alpha0 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri, alpha2_over_zetapluseta)
                primdata.alpha2_over_zetapluseta[0] = alpha1 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri, alpha3_over_zetapluseta)
                primdata.alpha3_over_zetapluseta[0] = alpha2 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri, alpha4_over_zetapluseta)
                primdata.alpha4_over_zetapluseta[0] = alpha3 * oogammapq;
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha1)
                primdata.rho12_over_alpha1[0] = alpha1 * oogammap;
#endif
#if LIBINT2_DEFINED(eri, rho12_over_alpha2)
                primdata.rho12_over_alpha2[0] = alpha0 * oogammap;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha3)
                primdata.rho34_over_alpha3[0] = alpha3 * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, rho34_over_alpha4)
                primdata.rho34_over_alpha4[0] = alpha2 * oogammaq;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_bra)
                primdata.two_alpha0_bra[0] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri, two_alpha0_ket)
                primdata.two_alpha0_ket[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_bra)
                primdata.two_alpha1_bra[0] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri, two_alpha1_ket)
                primdata.two_alpha1_ket[0] = 2.0 * alpha3;
#endif
              }
            }  // m != 0

            ++p;
          }  // prefac-based prim quartet screen

        }  // rough prim quartet screen based on pair values
      }    // ket prim pair
    }      // bra prim pair
    primdata_[0].contrdepth = p;
  }

#ifdef LIBINT2_ENGINE_TIMERS
  const auto t0 = timers.stop(0);
#ifdef LIBINT2_ENGINE_PROFILE_CLASS
  class_profiles[id].prereqs += t0.count();
  if (primdata_[0].contrdepth != 0) {
    class_profiles[id].nshellset += 1;
    class_profiles[id].nprimset += primdata_[0].contrdepth;
  }
#endif
#endif

  // all primitive combinations screened out? set 1st target ptr to nullptr
  if (primdata_[0].contrdepth == 0) {
    targets_[0] = nullptr;
    return targets_;
  }

  // compute directly (ss|ss)
  const auto compute_directly = lmax == 0 && deriv_order == 0;

  if (compute_directly) {
#ifdef LIBINT2_ENGINE_TIMERS
    timers.start(1);
#endif
    auto& stack = primdata_[0].stack[0];
    stack = 0;
    for (auto p = 0; p != primdata_[0].contrdepth; ++p)
      stack += primdata_[p].LIBINT_T_SS_EREP_SS(0)[0];
    primdata_[0].targets[0] = primdata_[0].stack;
#ifdef LIBINT2_ENGINE_TIMERS
    const auto t1 = timers.stop(1);
#ifdef LIBINT2_ENGINE_PROFILE_CLASS
    class_profiles[id].build_vrr += t1.count();
#endif
#endif
  }       // compute directly
  else {  // call libint
#ifdef LIBINT2_ENGINE_TIMERS
#ifdef LIBINT2_PROFILE
    const auto t1_hrr_start = primdata_[0].timers->read(0);
    const auto t1_vrr_start = primdata_[0].timers->read(1);
#endif
    timers.start(1);
#endif

    size_t buildfnidx;
    switch (braket_) {
      case BraKet::xx_xx:
        buildfnidx =
            ((bra1.contr[0].l * hard_lmax_ + bra2.contr[0].l) * hard_lmax_ +
             ket1.contr[0].l) *
                hard_lmax_ +
            ket2.contr[0].l;
        break;

      case BraKet::xx_xs: assert(false && "this braket is not supported"); break;
      case BraKet::xs_xx: {
        /// lmax might be center dependent
        int ket_lmax = hard_lmax_;
        switch (deriv_order_) {
          case 0:
#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri
            ket_lmax = hard_default_lmax_;
#endif
            break;
          case 1:
#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri1
            ket_lmax = hard_default_lmax_;
#endif
            break;
          case 2:
#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri2
            ket_lmax = hard_default_lmax_;
#endif
            break;
// TODO 3rd and 4th derivatives not tested
//          case 3:
//#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri3
//            ket_lmax = hard_default_lmax_;
//#endif
//            break;
//          case 4:
//#ifdef LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri4
//            ket_lmax = hard_default_lmax_;
//#endif
//            break;
//        default:assert(false && "deriv_order>4 not yet supported");
          default:assert(false && "deriv_order>2 not yet supported");
        }
        buildfnidx =
            (bra1.contr[0].l * ket_lmax + ket1.contr[0].l) * ket_lmax +
                ket2.contr[0].l;
#ifdef ERI3_PURE_SH
        if (bra1.contr[0].l > 1)
          assert(bra1.contr[0].pure &&
                 "library assumes a solid harmonics shell in bra of a 3-center "
                 "2-body int, but a cartesian shell given");
#endif
      } break;

      case BraKet::xs_xs:
        buildfnidx = bra1.contr[0].l * hard_lmax_ + ket1.contr[0].l;
#ifdef ERI2_PURE_SH
        if (bra1.contr[0].l > 1)
          assert(bra1.contr[0].pure &&
                 "library assumes solid harmonics shells in a 2-center "
                 "2-body int, but a cartesian shell given in bra");
        if (ket1.contr[0].l > 1)
          assert(ket1.contr[0].pure &&
                 "library assumes solid harmonics shells in a 2-center "
                 "2-body int, but a cartesian shell given in bra");
#endif
        break;

      default:
        assert(false && "invalid braket");
    }

    assert(buildfnptrs_[buildfnidx] && "null build function ptr");
    buildfnptrs_[buildfnidx](&primdata_[0]);

#ifdef LIBINT2_ENGINE_TIMERS
    const auto t1 = timers.stop(1);
#ifdef LIBINT2_ENGINE_PROFILE_CLASS
#ifndef LIBINT2_PROFILE
    class_profiles[id].build_vrr += t1.count();
#else
    class_profiles[id].build_hrr += primdata_[0].timers->read(0) - t1_hrr_start;
    class_profiles[id].build_vrr += primdata_[0].timers->read(1) - t1_vrr_start;
#endif
#endif
#endif

#ifdef LIBINT2_ENGINE_TIMERS
    timers.start(2);
#endif

    const auto ntargets = nshellsets();

    // if needed, permute and transform
    if (use_scratch) {
      constexpr auto using_scalar_real = sizeof(value_type) == sizeof(scalar_type);
      static_assert(using_scalar_real,
                    "Libint2 C++11 API only supports scalar real types");
      typedef Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic,
                            Eigen::RowMajor>
          Matrix;

      // a 2-d view of the 4-d source tensor
      const auto nr1_cart = bra1.cartesian_size();
      const auto nr2_cart = bra2.cartesian_size();
      const auto nc1_cart = ket1.cartesian_size();
      const auto nc2_cart = ket2.cartesian_size();
      const auto ncol_cart = nc1_cart * nc2_cart;
      const auto n1234_cart = nr1_cart * nr2_cart * ncol_cart;
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
      const auto n_tgt = nr1_tgt * nr2_tgt * ncol_tgt;

      auto hotscr = &scratch_[0];  // points to the hot scratch

      // transform to solid harmonics first, then unpermute, if necessary
      for (auto s = 0; s != ntargets; ++s) {
        // when permuting derivatives may need to permute shellsets also, not
        // just integrals
        // within shellsets; this will point to where source shellset s should end up
        auto s_target = s;

        auto source =
            primdata_[0].targets[s];  // points to the most recent result
        auto target = hotscr;

        if (bra1.contr[0].pure) {
          libint2::solidharmonics::transform_first(
              bra1.contr[0].l, nr2_cart * ncol_cart, source, target);
          std::swap(source, target);
        }
        if (bra2.contr[0].pure) {
          libint2::solidharmonics::transform_inner(bra1.size(), bra2.contr[0].l,
                                                   ncol_cart, source, target);
          std::swap(source, target);
        }
        if (ket1.contr[0].pure) {
          libint2::solidharmonics::transform_inner(nrow, ket1.contr[0].l,
                                                   nc2_cart, source, target);
          std::swap(source, target);
        }
        if (ket2.contr[0].pure) {
          libint2::solidharmonics::transform_last(
              bra1.size() * bra2.size() * ket1.size(), ket2.contr[0].l, source,
              target);
          std::swap(source, target);
        }

        // need to permute?
        if (permute) {
          // loop over rows of the source matrix
          const auto* src_row_ptr = source;
          auto tgt_ptr = target;

          // if permuting derivatives ints must update their derivative index
          // There is only one mapDerivIndex, and it is always the same number of dimensions for all deriv_orders.
          // TODO are there other BraKet types that need to be added here? Would require adding support to function 'generate_deriv_index_map'
          if (deriv_order > 0){
            switch(braket_) {
              case BraKet::xx_xx: {
                s_target = (*mapDerivIndex)[swap_braket][swap_tbra][swap_tket][s];
              }
                break;
              case BraKet::xs_xx: {
                assert(swap_bra == false);
                assert(swap_braket == false);
                s_target = (*mapDerivIndex)[0][0][swap_tket][s];
              }
                break;
              case BraKet::xs_xs: {
                assert(swap_bra == false);
                assert(swap_ket == false);
                assert(swap_braket == false);
                s_target = s;
              }
                break;

              default:
                assert(false && "This braket type not yet supported for geometric derivatives");
            }
          }
          // if permuting derivatives ints must update their derivative index
          //switch (deriv_order) {
          //  case 0:
          //    break;  // nothing to do

          //  case 1: {
          //    switch(braket_) {
          //      case BraKet::xx_xx: {
          //        const unsigned mapDerivIndex1_xxxx[2][2][2][12] = {
          //            {{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
          //              {0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8}},
          //             {{3, 4, 5, 0, 1, 2, 6, 7, 8, 9, 10, 11},
          //              {3, 4, 5, 0, 1, 2, 9, 10, 11, 6, 7, 8}}},
          //            {{{6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5},
          //              {9, 10, 11, 6, 7, 8, 0, 1, 2, 3, 4, 5}},
          //             {{6, 7, 8, 9, 10, 11, 3, 4, 5, 0, 1, 2},
          //              {9, 10, 11, 6, 7, 8, 3, 4, 5, 0, 1, 2}}}};
          //        s_target = mapDerivIndex1_xxxx[swap_braket][swap_tbra][swap_tket][s];
          //      }
          //        break;

          //      case BraKet::xs_xx: {
          //        assert(swap_bra == false);
          //        assert(swap_braket == false);
          //        const unsigned mapDerivIndex1_xsxx[2][9] = {
          //            {0,1,2,3,4,5,6,7,8},
          //            {0,1,2,6,7,8,3,4,5}
          //        };
          //        s_target = mapDerivIndex1_xsxx[swap_tket][s];
          //      }
          //        break;

          //      case BraKet::xs_xs: {
          //        assert(swap_bra == false);
          //        assert(swap_ket == false);
          //        assert(swap_braket == false);
          //        s_target = s;
          //      }
          //        break;

          //      default:
          //        assert(false && "this backet type not yet supported for 1st geometric derivatives");
          //    }
          //  } break;

          //  case 2: {
          //    switch(braket_) {
          //      case BraKet::xx_xx: {
          //        const unsigned mapDerivIndex2_xxxx[2][2][2][78] = {
          //            {{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
          //               13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
          //               26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
          //               39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
          //               52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64,
          //               65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77},
          //              {0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8, 12,
          //               13, 14, 15, 16, 20, 21, 22, 17, 18, 19, 23, 24, 25,
          //               26, 30, 31, 32, 27, 28, 29, 33, 34, 35, 39, 40, 41,
          //               36, 37, 38, 42, 43, 47, 48, 49, 44, 45, 46, 50, 54,
          //               55, 56, 51, 52, 53, 72, 73, 74, 60, 65, 69, 75, 76,
          //               61, 66, 70, 77, 62, 67, 71, 57, 58, 59, 63, 64, 68}},
          //             {{33, 34, 35, 3, 14, 24, 36, 37, 38, 39, 40, 41, 42,
          //               43, 4, 15, 25, 44, 45, 46, 47, 48, 49, 50, 5, 16,
          //               26, 51, 52, 53, 54, 55, 56, 0, 1, 2, 6, 7, 8,
          //               9, 10, 11, 12, 13, 17, 18, 19, 20, 21, 22, 23, 27,
          //               28, 29, 30, 31, 32, 57, 58, 59, 60, 61, 62, 63, 64,
          //               65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77},
          //              {33, 34, 35, 3, 14, 24, 39, 40, 41, 36, 37, 38, 42,
          //               43, 4, 15, 25, 47, 48, 49, 44, 45, 46, 50, 5, 16,
          //               26, 54, 55, 56, 51, 52, 53, 0, 1, 2, 9, 10, 11,
          //               6, 7, 8, 12, 13, 20, 21, 22, 17, 18, 19, 23, 30,
          //               31, 32, 27, 28, 29, 72, 73, 74, 60, 65, 69, 75, 76,
          //               61, 66, 70, 77, 62, 67, 71, 57, 58, 59, 63, 64, 68}}},
          //            {{{57, 58, 59, 60, 61, 62, 6, 17, 27, 36, 44, 51, 63,
          //               64, 65, 66, 67, 7, 18, 28, 37, 45, 52, 68, 69, 70,
          //               71, 8, 19, 29, 38, 46, 53, 72, 73, 74, 9, 20, 30,
          //               39, 47, 54, 75, 76, 10, 21, 31, 40, 48, 55, 77, 11,
          //               22, 32, 41, 49, 56, 0, 1, 2, 3, 4, 5, 12, 13,
          //               14, 15, 16, 23, 24, 25, 26, 33, 34, 35, 42, 43, 50},
          //              {72, 73, 74, 60, 65, 69, 9, 20, 30, 39, 47, 54, 75,
          //               76, 61, 66, 70, 10, 21, 31, 40, 48, 55, 77, 62, 67,
          //               71, 11, 22, 32, 41, 49, 56, 57, 58, 59, 6, 17, 27,
          //               36, 44, 51, 63, 64, 7, 18, 28, 37, 45, 52, 68, 8,
          //               19, 29, 38, 46, 53, 0, 1, 2, 3, 4, 5, 12, 13,
          //               14, 15, 16, 23, 24, 25, 26, 33, 34, 35, 42, 43, 50}},
          //             {{57, 58, 59, 60, 61, 62, 36, 44, 51, 6, 17, 27, 63,
          //               64, 65, 66, 67, 37, 45, 52, 7, 18, 28, 68, 69, 70,
          //               71, 38, 46, 53, 8, 19, 29, 72, 73, 74, 39, 47, 54,
          //               9, 20, 30, 75, 76, 40, 48, 55, 10, 21, 31, 77, 41,
          //               49, 56, 11, 22, 32, 33, 34, 35, 3, 14, 24, 42, 43,
          //               4, 15, 25, 50, 5, 16, 26, 0, 1, 2, 12, 13, 23},
          //              {72, 73, 74, 60, 65, 69, 39, 47, 54, 9, 20, 30, 75,
          //               76, 61, 66, 70, 40, 48, 55, 10, 21, 31, 77, 62, 67,
          //               71, 41, 49, 56, 11, 22, 32, 57, 58, 59, 36, 44, 51,
          //               6, 17, 27, 63, 64, 37, 45, 52, 7, 18, 28, 68, 38,
          //               46, 53, 8, 19, 29, 33, 34, 35, 3, 14, 24, 42, 43,
          //               4, 15, 25, 50, 5, 16, 26, 0, 1, 2, 12, 13, 23}}}};
          //        s_target = mapDerivIndex2_xxxx[swap_braket][swap_tbra][swap_tket][s];
          //      }
          //        break;

          //      case BraKet::xs_xx: {
          //        assert(swap_bra == false);
          //        assert(swap_braket == false);
          //        const unsigned mapDerivIndex2_xsxx[2][45] = {
          //            {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
          //             12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
          //             24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
          //             36, 37, 38, 39, 40, 41, 42, 43, 44},
          //            {0,  1,  2,  6,  7,  8,  3,  4,  5,  9,  10, 14,
          //             15, 16, 11, 12, 13, 17, 21, 22, 23, 18, 19, 20,
          //             39, 40, 41, 27, 32, 36, 42, 43, 28, 33, 37, 44,
          //             29, 34, 38, 24, 25, 26, 30, 31, 35}};
          //        s_target = mapDerivIndex2_xsxx[swap_tket][s];
          //      }
          //        break;

          //      case BraKet::xs_xs: {
          //        assert(swap_bra == false);
          //        assert(swap_ket == false);
          //        assert(swap_braket == false);
          //        s_target = s;
          //      }
          //        break;

          //      default:
          //        assert(false && "this backet type not yet supported for 2st geometric derivatives");
          //    }
          //  } break;
          //  // deriv_order == 3
          //  case 3: {
          //    switch(braket_) {
          //      case BraKet::xx_xx: {
          //        const unsigned mapDerivIndex3_xxxx[2][2][2][364] = {
          //           {{{ 0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
          //               13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
          //               26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
          //               39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
          //               52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
          //               65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
          //               78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
          //               91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
          //              104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
          //              117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
          //              130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
          //              143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
          //              156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
          //              169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
          //              182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
          //              195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
          //              208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
          //              221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
          //              234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
          //              247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,
          //              260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
          //              273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285,
          //              286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298,
          //              299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311,
          //              312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324,
          //              325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337,
          //              338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350,
          //              351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363},
          //             {  0,   1,   2,   3,   4,   5,   9,  10,  11,   6,   7,   8,  12,
          //               13,  14,  15,  16,  20,  21,  22,  17,  18,  19,  23,  24,  25,
          //               26,  30,  31,  32,  27,  28,  29,  33,  34,  35,  39,  40,  41,
          //               36,  37,  38,  42,  43,  47,  48,  49,  44,  45,  46,  50,  54,
          //               55,  56,  51,  52,  53,  72,  73,  74,  60,  65,  69,  75,  76,
          //               61,  66,  70,  77,  62,  67,  71,  57,  58,  59,  63,  64,  68,
          //               78,  79,  80,  81,  82,  86,  87,  88,  83,  84,  85,  89,  90,
          //               91,  92,  96,  97,  98,  93,  94,  95,  99, 100, 101, 105, 106,
          //              107, 102, 103, 104, 108, 109, 113, 114, 115, 110, 111, 112, 116,
          //              120, 121, 122, 117, 118, 119, 138, 139, 140, 126, 131, 135, 141,
          //              142, 127, 132, 136, 143, 128, 133, 137, 123, 124, 125, 129, 130,
          //              134, 144, 145, 146, 147, 151, 152, 153, 148, 149, 150, 154, 155,
          //              156, 160, 161, 162, 157, 158, 159, 163, 164, 168, 169, 170, 165,
          //              166, 167, 171, 175, 176, 177, 172, 173, 174, 193, 194, 195, 181,
          //              186, 190, 196, 197, 182, 187, 191, 198, 183, 188, 192, 178, 179,
          //              180, 184, 185, 189, 199, 200, 201, 205, 206, 207, 202, 203, 204,
          //              208, 209, 213, 214, 215, 210, 211, 212, 216, 220, 221, 222, 217,
          //              218, 219, 238, 239, 240, 226, 231, 235, 241, 242, 227, 232, 236,
          //              243, 228, 233, 237, 223, 224, 225, 229, 230, 234, 244, 245, 249,
          //              250, 251, 246, 247, 248, 252, 256, 257, 258, 253, 254, 255, 274,
          //              275, 276, 262, 267, 271, 277, 278, 263, 268, 272, 279, 264, 269,
          //              273, 259, 260, 261, 265, 266, 270, 280, 284, 285, 286, 281, 282,
          //              283, 302, 303, 304, 290, 295, 299, 305, 306, 291, 296, 300, 307,
          //              292, 297, 301, 287, 288, 289, 293, 294, 298, 354, 355, 356, 323,
          //              338, 348, 357, 358, 324, 339, 349, 359, 325, 340, 350, 311, 316,
          //              320, 331, 335, 345, 360, 361, 326, 341, 351, 362, 327, 342, 352,
          //              312, 317, 321, 332, 336, 346, 363, 328, 343, 353, 313, 318, 322,
          //              333, 337, 347, 308, 309, 310, 314, 315, 319, 329, 330, 334, 344}},
          //            {{199, 200, 201,  33,  99, 154, 202, 203, 204, 205, 206, 207, 208,
          //              209,  34, 100, 155, 210, 211, 212, 213, 214, 215, 216,  35, 101,
          //              156, 217, 218, 219, 220, 221, 222,   3,  14,  24,  36,  37,  38,
          //               39,  40,  41,  80,  90, 102, 103, 104, 105, 106, 107, 145, 157,
          //              158, 159, 160, 161, 162, 223, 224, 225, 226, 227, 228, 229, 230,
          //              231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243,
          //              244, 245,  42, 108, 163, 246, 247, 248, 249, 250, 251, 252,  43,
          //              109, 164, 253, 254, 255, 256, 257, 258,   4,  15,  25,  44,  45,
          //               46,  47,  48,  49,  81,  91, 110, 111, 112, 113, 114, 115, 146,
          //              165, 166, 167, 168, 169, 170, 259, 260, 261, 262, 263, 264, 265,
          //              266, 267, 268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278,
          //              279, 280,  50, 116, 171, 281, 282, 283, 284, 285, 286,   5,  16,
          //               26,  51,  52,  53,  54,  55,  56,  82,  92, 117, 118, 119, 120,
          //              121, 122, 147, 172, 173, 174, 175, 176, 177, 287, 288, 289, 290,
          //              291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 301, 302, 303,
          //              304, 305, 306, 307,   0,   1,   2,   6,   7,   8,   9,  10,  11,
          //               12,  13,  17,  18,  19,  20,  21,  22,  23,  27,  28,  29,  30,
          //               31,  32,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,
          //               68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  83,
          //               84,  85,  86,  87,  88,  89,  93,  94,  95,  96,  97,  98, 123,
          //              124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136,
          //              137, 138, 139, 140, 141, 142, 143, 144, 148, 149, 150, 151, 152,
          //              153, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189,
          //              190, 191, 192, 193, 194, 195, 196, 197, 198, 308, 309, 310, 311,
          //              312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324,
          //              325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337,
          //              338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350,
          //              351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363},
          //             {199, 200, 201,  33,  99, 154, 205, 206, 207, 202, 203, 204, 208,
          //              209,  34, 100, 155, 213, 214, 215, 210, 211, 212, 216,  35, 101,
          //              156, 220, 221, 222, 217, 218, 219,   3,  14,  24,  39,  40,  41,
          //               36,  37,  38,  80,  90, 105, 106, 107, 102, 103, 104, 145, 160,
          //              161, 162, 157, 158, 159, 238, 239, 240, 226, 231, 235, 241, 242,
          //              227, 232, 236, 243, 228, 233, 237, 223, 224, 225, 229, 230, 234,
          //              244, 245,  42, 108, 163, 249, 250, 251, 246, 247, 248, 252,  43,
          //              109, 164, 256, 257, 258, 253, 254, 255,   4,  15,  25,  47,  48,
          //               49,  44,  45,  46,  81,  91, 113, 114, 115, 110, 111, 112, 146,
          //              168, 169, 170, 165, 166, 167, 274, 275, 276, 262, 267, 271, 277,
          //              278, 263, 268, 272, 279, 264, 269, 273, 259, 260, 261, 265, 266,
          //              270, 280,  50, 116, 171, 284, 285, 286, 281, 282, 283,   5,  16,
          //               26,  54,  55,  56,  51,  52,  53,  82,  92, 120, 121, 122, 117,
          //              118, 119, 147, 175, 176, 177, 172, 173, 174, 302, 303, 304, 290,
          //              295, 299, 305, 306, 291, 296, 300, 307, 292, 297, 301, 287, 288,
          //              289, 293, 294, 298,   0,   1,   2,   9,  10,  11,   6,   7,   8,
          //               12,  13,  20,  21,  22,  17,  18,  19,  23,  30,  31,  32,  27,
          //               28,  29,  72,  73,  74,  60,  65,  69,  75,  76,  61,  66,  70,
          //               77,  62,  67,  71,  57,  58,  59,  63,  64,  68,  78,  79,  86,
          //               87,  88,  83,  84,  85,  89,  96,  97,  98,  93,  94,  95, 138,
          //              139, 140, 126, 131, 135, 141, 142, 127, 132, 136, 143, 128, 133,
          //              137, 123, 124, 125, 129, 130, 134, 144, 151, 152, 153, 148, 149,
          //              150, 193, 194, 195, 181, 186, 190, 196, 197, 182, 187, 191, 198,
          //              183, 188, 192, 178, 179, 180, 184, 185, 189, 354, 355, 356, 323,
          //              338, 348, 357, 358, 324, 339, 349, 359, 325, 340, 350, 311, 316,
          //              320, 331, 335, 345, 360, 361, 326, 341, 351, 362, 327, 342, 352,
          //              312, 317, 321, 332, 336, 346, 363, 328, 343, 353, 313, 318, 322,
          //              333, 337, 347, 308, 309, 310, 314, 315, 319, 329, 330, 334, 344}}},
          //           {{{308, 309, 310, 311, 312, 313,  57, 123, 178, 223, 259, 287, 314,
          //              315, 316, 317, 318,  58, 124, 179, 224, 260, 288, 319, 320, 321,
          //              322,  59, 125, 180, 225, 261, 289, 323, 324, 325,  60, 126, 181,
          //              226, 262, 290, 326, 327,  61, 127, 182, 227, 263, 291, 328,  62,
          //              128, 183, 228, 264, 292,   6,  17,  27,  36,  44,  51,  83,  93,
          //              102, 110, 117, 148, 157, 165, 172, 202, 210, 217, 246, 253, 281,
          //              329, 330, 331, 332, 333,  63, 129, 184, 229, 265, 293, 334, 335,
          //              336, 337,  64, 130, 185, 230, 266, 294, 338, 339, 340,  65, 131,
          //              186, 231, 267, 295, 341, 342,  66, 132, 187, 232, 268, 296, 343,
          //               67, 133, 188, 233, 269, 297,   7,  18,  28,  37,  45,  52,  84,
          //               94, 103, 111, 118, 149, 158, 166, 173, 203, 211, 218, 247, 254,
          //              282, 344, 345, 346, 347,  68, 134, 189, 234, 270, 298, 348, 349,
          //              350,  69, 135, 190, 235, 271, 299, 351, 352,  70, 136, 191, 236,
          //              272, 300, 353,  71, 137, 192, 237, 273, 301,   8,  19,  29,  38,
          //               46,  53,  85,  95, 104, 112, 119, 150, 159, 167, 174, 204, 212,
          //              219, 248, 255, 283, 354, 355, 356,  72, 138, 193, 238, 274, 302,
          //              357, 358,  73, 139, 194, 239, 275, 303, 359,  74, 140, 195, 240,
          //              276, 304,   9,  20,  30,  39,  47,  54,  86,  96, 105, 113, 120,
          //              151, 160, 168, 175, 205, 213, 220, 249, 256, 284, 360, 361,  75,
          //              141, 196, 241, 277, 305, 362,  76, 142, 197, 242, 278, 306,  10,
          //               21,  31,  40,  48,  55,  87,  97, 106, 114, 121, 152, 161, 169,
          //              176, 206, 214, 221, 250, 257, 285, 363,  77, 143, 198, 243, 279,
          //              307,  11,  22,  32,  41,  49,  56,  88,  98, 107, 115, 122, 153,
          //              162, 170, 177, 207, 215, 222, 251, 258, 286,   0,   1,   2,   3,
          //                4,   5,  12,  13,  14,  15,  16,  23,  24,  25,  26,  33,  34,
          //               35,  42,  43,  50,  78,  79,  80,  81,  82,  89,  90,  91,  92,
          //               99, 100, 101, 108, 109, 116, 144, 145, 146, 147, 154, 155, 156,
          //              163, 164, 171, 199, 200, 201, 208, 209, 216, 244, 245, 252, 280},
          //             {354, 355, 356, 323, 338, 348,  72, 138, 193, 238, 274, 302, 357,
          //              358, 324, 339, 349,  73, 139, 194, 239, 275, 303, 359, 325, 340,
          //              350,  74, 140, 195, 240, 276, 304, 311, 316, 320,  60, 126, 181,
          //              226, 262, 290, 331, 335,  65, 131, 186, 231, 267, 295, 345,  69,
          //              135, 190, 235, 271, 299,   9,  20,  30,  39,  47,  54,  86,  96,
          //              105, 113, 120, 151, 160, 168, 175, 205, 213, 220, 249, 256, 284,
          //              360, 361, 326, 341, 351,  75, 141, 196, 241, 277, 305, 362, 327,
          //              342, 352,  76, 142, 197, 242, 278, 306, 312, 317, 321,  61, 127,
          //              182, 227, 263, 291, 332, 336,  66, 132, 187, 232, 268, 296, 346,
          //               70, 136, 191, 236, 272, 300,  10,  21,  31,  40,  48,  55,  87,
          //               97, 106, 114, 121, 152, 161, 169, 176, 206, 214, 221, 250, 257,
          //              285, 363, 328, 343, 353,  77, 143, 198, 243, 279, 307, 313, 318,
          //              322,  62, 128, 183, 228, 264, 292, 333, 337,  67, 133, 188, 233,
          //              269, 297, 347,  71, 137, 192, 237, 273, 301,  11,  22,  32,  41,
          //               49,  56,  88,  98, 107, 115, 122, 153, 162, 170, 177, 207, 215,
          //              222, 251, 258, 286, 308, 309, 310,  57, 123, 178, 223, 259, 287,
          //              314, 315,  58, 124, 179, 224, 260, 288, 319,  59, 125, 180, 225,
          //              261, 289,   6,  17,  27,  36,  44,  51,  83,  93, 102, 110, 117,
          //              148, 157, 165, 172, 202, 210, 217, 246, 253, 281, 329, 330,  63,
          //              129, 184, 229, 265, 293, 334,  64, 130, 185, 230, 266, 294,   7,
          //               18,  28,  37,  45,  52,  84,  94, 103, 111, 118, 149, 158, 166,
          //              173, 203, 211, 218, 247, 254, 282, 344,  68, 134, 189, 234, 270,
          //              298,   8,  19,  29,  38,  46,  53,  85,  95, 104, 112, 119, 150,
          //              159, 167, 174, 204, 212, 219, 248, 255, 283,   0,   1,   2,   3,
          //                4,   5,  12,  13,  14,  15,  16,  23,  24,  25,  26,  33,  34,
          //               35,  42,  43,  50,  78,  79,  80,  81,  82,  89,  90,  91,  92,
          //               99, 100, 101, 108, 109, 116, 144, 145, 146, 147, 154, 155, 156,
          //              163, 164, 171, 199, 200, 201, 208, 209, 216, 244, 245, 252, 280}},
          //            {{308, 309, 310, 311, 312, 313, 223, 259, 287,  57, 123, 178, 314,
          //              315, 316, 317, 318, 224, 260, 288,  58, 124, 179, 319, 320, 321,
          //              322, 225, 261, 289,  59, 125, 180, 323, 324, 325, 226, 262, 290,
          //               60, 126, 181, 326, 327, 227, 263, 291,  61, 127, 182, 328, 228,
          //              264, 292,  62, 128, 183, 202, 210, 217,  36, 102, 157, 246, 253,
          //               44, 110, 165, 281,  51, 117, 172,   6,  17,  27,  83,  93, 148,
          //              329, 330, 331, 332, 333, 229, 265, 293,  63, 129, 184, 334, 335,
          //              336, 337, 230, 266, 294,  64, 130, 185, 338, 339, 340, 231, 267,
          //              295,  65, 131, 186, 341, 342, 232, 268, 296,  66, 132, 187, 343,
          //              233, 269, 297,  67, 133, 188, 203, 211, 218,  37, 103, 158, 247,
          //              254,  45, 111, 166, 282,  52, 118, 173,   7,  18,  28,  84,  94,
          //              149, 344, 345, 346, 347, 234, 270, 298,  68, 134, 189, 348, 349,
          //              350, 235, 271, 299,  69, 135, 190, 351, 352, 236, 272, 300,  70,
          //              136, 191, 353, 237, 273, 301,  71, 137, 192, 204, 212, 219,  38,
          //              104, 159, 248, 255,  46, 112, 167, 283,  53, 119, 174,   8,  19,
          //               29,  85,  95, 150, 354, 355, 356, 238, 274, 302,  72, 138, 193,
          //              357, 358, 239, 275, 303,  73, 139, 194, 359, 240, 276, 304,  74,
          //              140, 195, 205, 213, 220,  39, 105, 160, 249, 256,  47, 113, 168,
          //              284,  54, 120, 175,   9,  20,  30,  86,  96, 151, 360, 361, 241,
          //              277, 305,  75, 141, 196, 362, 242, 278, 306,  76, 142, 197, 206,
          //              214, 221,  40, 106, 161, 250, 257,  48, 114, 169, 285,  55, 121,
          //              176,  10,  21,  31,  87,  97, 152, 363, 243, 279, 307,  77, 143,
          //              198, 207, 215, 222,  41, 107, 162, 251, 258,  49, 115, 170, 286,
          //               56, 122, 177,  11,  22,  32,  88,  98, 153, 199, 200, 201,  33,
          //               99, 154, 208, 209,  34, 100, 155, 216,  35, 101, 156,   3,  14,
          //               24,  80,  90, 145, 244, 245,  42, 108, 163, 252,  43, 109, 164,
          //                4,  15,  25,  81,  91, 146, 280,  50, 116, 171,   5,  16,  26,
          //               82,  92, 147,   0,   1,   2,  12,  13,  23,  78,  79,  89, 144},
          //             {354, 355, 356, 323, 338, 348, 238, 274, 302,  72, 138, 193, 357,
          //              358, 324, 339, 349, 239, 275, 303,  73, 139, 194, 359, 325, 340,
          //              350, 240, 276, 304,  74, 140, 195, 311, 316, 320, 226, 262, 290,
          //               60, 126, 181, 331, 335, 231, 267, 295,  65, 131, 186, 345, 235,
          //              271, 299,  69, 135, 190, 205, 213, 220,  39, 105, 160, 249, 256,
          //               47, 113, 168, 284,  54, 120, 175,   9,  20,  30,  86,  96, 151,
          //              360, 361, 326, 341, 351, 241, 277, 305,  75, 141, 196, 362, 327,
          //              342, 352, 242, 278, 306,  76, 142, 197, 312, 317, 321, 227, 263,
          //              291,  61, 127, 182, 332, 336, 232, 268, 296,  66, 132, 187, 346,
          //              236, 272, 300,  70, 136, 191, 206, 214, 221,  40, 106, 161, 250,
          //              257,  48, 114, 169, 285,  55, 121, 176,  10,  21,  31,  87,  97,
          //              152, 363, 328, 343, 353, 243, 279, 307,  77, 143, 198, 313, 318,
          //              322, 228, 264, 292,  62, 128, 183, 333, 337, 233, 269, 297,  67,
          //              133, 188, 347, 237, 273, 301,  71, 137, 192, 207, 215, 222,  41,
          //              107, 162, 251, 258,  49, 115, 170, 286,  56, 122, 177,  11,  22,
          //               32,  88,  98, 153, 308, 309, 310, 223, 259, 287,  57, 123, 178,
          //              314, 315, 224, 260, 288,  58, 124, 179, 319, 225, 261, 289,  59,
          //              125, 180, 202, 210, 217,  36, 102, 157, 246, 253,  44, 110, 165,
          //              281,  51, 117, 172,   6,  17,  27,  83,  93, 148, 329, 330, 229,
          //              265, 293,  63, 129, 184, 334, 230, 266, 294,  64, 130, 185, 203,
          //              211, 218,  37, 103, 158, 247, 254,  45, 111, 166, 282,  52, 118,
          //              173,   7,  18,  28,  84,  94, 149, 344, 234, 270, 298,  68, 134,
          //              189, 204, 212, 219,  38, 104, 159, 248, 255,  46, 112, 167, 283,
          //               53, 119, 174,   8,  19,  29,  85,  95, 150, 199, 200, 201,  33,
          //               99, 154, 208, 209,  34, 100, 155, 216,  35, 101, 156,   3,  14,
          //               24,  80,  90, 145, 244, 245,  42, 108, 163, 252,  43, 109, 164,
          //                4,  15,  25,  81,  91, 146, 280,  50, 116, 171,   5,  16,  26,
          //               82,  92, 147,   0,   1,   2,  12,  13,  23,  78,  79,  89, 144}}}};
          //        s_target = mapDerivIndex3_xxxx[swap_braket][swap_tbra][swap_tket][s];
          //      }
          //        break;

          //      // TODO 3center and 2center -- unable to test these, so commented out, but the array should be right. 
          //      // Requires definitions of LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri3/_3eri4 to work above
          //      //case BraKet::xs_xx: {
          //      //  assert(swap_bra == false);
          //      //  assert(swap_braket == false);
          //      //  const unsigned mapDerivIndex3_xsxx[2][165] = {
          //      //       {  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
          //      //         13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
          //      //         26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
          //      //         39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
          //      //         52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
          //      //         65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
          //      //         78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
          //      //         91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
          //      //        104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
          //      //        117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
          //      //        130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
          //      //        143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
          //      //        156, 157, 158, 159, 160, 161, 162, 163, 164},
          //      //       {  0,   1,   2,   6,   7,   8,   3,   4,   5,   9,  10,  14,  15,
          //      //         16,  11,  12,  13,  17,  21,  22,  23,  18,  19,  20,  39,  40,
          //      //         41,  27,  32,  36,  42,  43,  28,  33,  37,  44,  29,  34,  38,
          //      //         24,  25,  26,  30,  31,  35,  45,  46,  50,  51,  52,  47,  48,
          //      //         49,  53,  57,  58,  59,  54,  55,  56,  75,  76,  77,  63,  68,
          //      //         72,  78,  79,  64,  69,  73,  80,  65,  70,  74,  60,  61,  62,
          //      //         66,  67,  71,  81,  85,  86,  87,  82,  83,  84, 103, 104, 105,
          //      //         91,  96, 100, 106, 107,  92,  97, 101, 108,  93,  98, 102,  88,
          //      //         89,  90,  94,  95,  99, 155, 156, 157, 124, 139, 149, 158, 159,
          //      //        125, 140, 150, 160, 126, 141, 151, 112, 117, 121, 132, 136, 146,
          //      //        161, 162, 127, 142, 152, 163, 128, 143, 153, 113, 118, 122, 133,
          //      //        137, 147, 164, 129, 144, 154, 114, 119, 123, 134, 138, 148, 109,
          //      //        110, 111, 115, 116, 120, 130, 131, 135, 145}};
          //      //  s_target = mapDerivIndex3_xsxx[swap_tket][s];
          //      //}
          //      //  break;

          //      //case BraKet::xs_xs: {
          //      //  assert(swap_bra == false);
          //      //  assert(swap_ket == false);
          //      //  assert(swap_braket == false);
          //      //  s_target = s;
          //      //}
          //      //  break;

          //      default:
          //        assert(false && "this backet type not yet supported for 3rd geometric derivatives");
          //    }
          //  } break;

          //  // deriv_order == 4
          //  case 4: {
          //    switch(braket_) {
          //      case BraKet::xx_xx: {
          //        const unsigned mapDerivIndex4_xxxx[2][2][2][1365] = {
          //            {{{   0,    1,    2,    3,    4,    5,    6,    7,    8,    9,   10,
          //                 11,   12,   13,   14,   15,   16,   17,   18,   19,   20,   21,
          //                 22,   23,   24,   25,   26,   27,   28,   29,   30,   31,   32,
          //                 33,   34,   35,   36,   37,   38,   39,   40,   41,   42,   43,
          //                 44,   45,   46,   47,   48,   49,   50,   51,   52,   53,   54,
          //                 55,   56,   57,   58,   59,   60,   61,   62,   63,   64,   65,
          //                 66,   67,   68,   69,   70,   71,   72,   73,   74,   75,   76,
          //                 77,   78,   79,   80,   81,   82,   83,   84,   85,   86,   87,
          //                 88,   89,   90,   91,   92,   93,   94,   95,   96,   97,   98,
          //                 99,  100,  101,  102,  103,  104,  105,  106,  107,  108,  109,
          //                110,  111,  112,  113,  114,  115,  116,  117,  118,  119,  120,
          //                121,  122,  123,  124,  125,  126,  127,  128,  129,  130,  131,
          //                132,  133,  134,  135,  136,  137,  138,  139,  140,  141,  142,
          //                143,  144,  145,  146,  147,  148,  149,  150,  151,  152,  153,
          //                154,  155,  156,  157,  158,  159,  160,  161,  162,  163,  164,
          //                165,  166,  167,  168,  169,  170,  171,  172,  173,  174,  175,
          //                176,  177,  178,  179,  180,  181,  182,  183,  184,  185,  186,
          //                187,  188,  189,  190,  191,  192,  193,  194,  195,  196,  197,
          //                198,  199,  200,  201,  202,  203,  204,  205,  206,  207,  208,
          //                209,  210,  211,  212,  213,  214,  215,  216,  217,  218,  219,
          //                220,  221,  222,  223,  224,  225,  226,  227,  228,  229,  230,
          //                231,  232,  233,  234,  235,  236,  237,  238,  239,  240,  241,
          //                242,  243,  244,  245,  246,  247,  248,  249,  250,  251,  252,
          //                253,  254,  255,  256,  257,  258,  259,  260,  261,  262,  263,
          //                264,  265,  266,  267,  268,  269,  270,  271,  272,  273,  274,
          //                275,  276,  277,  278,  279,  280,  281,  282,  283,  284,  285,
          //                286,  287,  288,  289,  290,  291,  292,  293,  294,  295,  296,
          //                297,  298,  299,  300,  301,  302,  303,  304,  305,  306,  307,
          //                308,  309,  310,  311,  312,  313,  314,  315,  316,  317,  318,
          //                319,  320,  321,  322,  323,  324,  325,  326,  327,  328,  329,
          //                330,  331,  332,  333,  334,  335,  336,  337,  338,  339,  340,
          //                341,  342,  343,  344,  345,  346,  347,  348,  349,  350,  351,
          //                352,  353,  354,  355,  356,  357,  358,  359,  360,  361,  362,
          //                363,  364,  365,  366,  367,  368,  369,  370,  371,  372,  373,
          //                374,  375,  376,  377,  378,  379,  380,  381,  382,  383,  384,
          //                385,  386,  387,  388,  389,  390,  391,  392,  393,  394,  395,
          //                396,  397,  398,  399,  400,  401,  402,  403,  404,  405,  406,
          //                407,  408,  409,  410,  411,  412,  413,  414,  415,  416,  417,
          //                418,  419,  420,  421,  422,  423,  424,  425,  426,  427,  428,
          //                429,  430,  431,  432,  433,  434,  435,  436,  437,  438,  439,
          //                440,  441,  442,  443,  444,  445,  446,  447,  448,  449,  450,
          //                451,  452,  453,  454,  455,  456,  457,  458,  459,  460,  461,
          //                462,  463,  464,  465,  466,  467,  468,  469,  470,  471,  472,
          //                473,  474,  475,  476,  477,  478,  479,  480,  481,  482,  483,
          //                484,  485,  486,  487,  488,  489,  490,  491,  492,  493,  494,
          //                495,  496,  497,  498,  499,  500,  501,  502,  503,  504,  505,
          //                506,  507,  508,  509,  510,  511,  512,  513,  514,  515,  516,
          //                517,  518,  519,  520,  521,  522,  523,  524,  525,  526,  527,
          //                528,  529,  530,  531,  532,  533,  534,  535,  536,  537,  538,
          //                539,  540,  541,  542,  543,  544,  545,  546,  547,  548,  549,
          //                550,  551,  552,  553,  554,  555,  556,  557,  558,  559,  560,
          //                561,  562,  563,  564,  565,  566,  567,  568,  569,  570,  571,
          //                572,  573,  574,  575,  576,  577,  578,  579,  580,  581,  582,
          //                583,  584,  585,  586,  587,  588,  589,  590,  591,  592,  593,
          //                594,  595,  596,  597,  598,  599,  600,  601,  602,  603,  604,
          //                605,  606,  607,  608,  609,  610,  611,  612,  613,  614,  615,
          //                616,  617,  618,  619,  620,  621,  622,  623,  624,  625,  626,
          //                627,  628,  629,  630,  631,  632,  633,  634,  635,  636,  637,
          //                638,  639,  640,  641,  642,  643,  644,  645,  646,  647,  648,
          //                649,  650,  651,  652,  653,  654,  655,  656,  657,  658,  659,
          //                660,  661,  662,  663,  664,  665,  666,  667,  668,  669,  670,
          //                671,  672,  673,  674,  675,  676,  677,  678,  679,  680,  681,
          //                682,  683,  684,  685,  686,  687,  688,  689,  690,  691,  692,
          //                693,  694,  695,  696,  697,  698,  699,  700,  701,  702,  703,
          //                704,  705,  706,  707,  708,  709,  710,  711,  712,  713,  714,
          //                715,  716,  717,  718,  719,  720,  721,  722,  723,  724,  725,
          //                726,  727,  728,  729,  730,  731,  732,  733,  734,  735,  736,
          //                737,  738,  739,  740,  741,  742,  743,  744,  745,  746,  747,
          //                748,  749,  750,  751,  752,  753,  754,  755,  756,  757,  758,
          //                759,  760,  761,  762,  763,  764,  765,  766,  767,  768,  769,
          //                770,  771,  772,  773,  774,  775,  776,  777,  778,  779,  780,
          //                781,  782,  783,  784,  785,  786,  787,  788,  789,  790,  791,
          //                792,  793,  794,  795,  796,  797,  798,  799,  800,  801,  802,
          //                803,  804,  805,  806,  807,  808,  809,  810,  811,  812,  813,
          //                814,  815,  816,  817,  818,  819,  820,  821,  822,  823,  824,
          //                825,  826,  827,  828,  829,  830,  831,  832,  833,  834,  835,
          //                836,  837,  838,  839,  840,  841,  842,  843,  844,  845,  846,
          //                847,  848,  849,  850,  851,  852,  853,  854,  855,  856,  857,
          //                858,  859,  860,  861,  862,  863,  864,  865,  866,  867,  868,
          //                869,  870,  871,  872,  873,  874,  875,  876,  877,  878,  879,
          //                880,  881,  882,  883,  884,  885,  886,  887,  888,  889,  890,
          //                891,  892,  893,  894,  895,  896,  897,  898,  899,  900,  901,
          //                902,  903,  904,  905,  906,  907,  908,  909,  910,  911,  912,
          //                913,  914,  915,  916,  917,  918,  919,  920,  921,  922,  923,
          //                924,  925,  926,  927,  928,  929,  930,  931,  932,  933,  934,
          //                935,  936,  937,  938,  939,  940,  941,  942,  943,  944,  945,
          //                946,  947,  948,  949,  950,  951,  952,  953,  954,  955,  956,
          //                957,  958,  959,  960,  961,  962,  963,  964,  965,  966,  967,
          //                968,  969,  970,  971,  972,  973,  974,  975,  976,  977,  978,
          //                979,  980,  981,  982,  983,  984,  985,  986,  987,  988,  989,
          //                990,  991,  992,  993,  994,  995,  996,  997,  998,  999, 1000,
          //               1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011,
          //               1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022,
          //               1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033,
          //               1034, 1035, 1036, 1037, 1038, 1039, 1040, 1041, 1042, 1043, 1044,
          //               1045, 1046, 1047, 1048, 1049, 1050, 1051, 1052, 1053, 1054, 1055,
          //               1056, 1057, 1058, 1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066,
          //               1067, 1068, 1069, 1070, 1071, 1072, 1073, 1074, 1075, 1076, 1077,
          //               1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086, 1087, 1088,
          //               1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097, 1098, 1099,
          //               1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109, 1110,
          //               1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120, 1121,
          //               1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131, 1132,
          //               1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142, 1143,
          //               1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153, 1154,
          //               1155, 1156, 1157, 1158, 1159, 1160, 1161, 1162, 1163, 1164, 1165,
          //               1166, 1167, 1168, 1169, 1170, 1171, 1172, 1173, 1174, 1175, 1176,
          //               1177, 1178, 1179, 1180, 1181, 1182, 1183, 1184, 1185, 1186, 1187,
          //               1188, 1189, 1190, 1191, 1192, 1193, 1194, 1195, 1196, 1197, 1198,
          //               1199, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1209,
          //               1210, 1211, 1212, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220,
          //               1221, 1222, 1223, 1224, 1225, 1226, 1227, 1228, 1229, 1230, 1231,
          //               1232, 1233, 1234, 1235, 1236, 1237, 1238, 1239, 1240, 1241, 1242,
          //               1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253,
          //               1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264,
          //               1265, 1266, 1267, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275,
          //               1276, 1277, 1278, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286,
          //               1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297,
          //               1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308,
          //               1309, 1310, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319,
          //               1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330,
          //               1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341,
          //               1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1351, 1352,
          //               1353, 1354, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 1362, 1363,
          //               1364},
          //              {   0,    1,    2,    3,    4,    5,    9,   10,   11,    6,    7,
          //                  8,   12,   13,   14,   15,   16,   20,   21,   22,   17,   18,
          //                 19,   23,   24,   25,   26,   30,   31,   32,   27,   28,   29,
          //                 33,   34,   35,   39,   40,   41,   36,   37,   38,   42,   43,
          //                 47,   48,   49,   44,   45,   46,   50,   54,   55,   56,   51,
          //                 52,   53,   72,   73,   74,   60,   65,   69,   75,   76,   61,
          //                 66,   70,   77,   62,   67,   71,   57,   58,   59,   63,   64,
          //                 68,   78,   79,   80,   81,   82,   86,   87,   88,   83,   84,
          //                 85,   89,   90,   91,   92,   96,   97,   98,   93,   94,   95,
          //                 99,  100,  101,  105,  106,  107,  102,  103,  104,  108,  109,
          //                113,  114,  115,  110,  111,  112,  116,  120,  121,  122,  117,
          //                118,  119,  138,  139,  140,  126,  131,  135,  141,  142,  127,
          //                132,  136,  143,  128,  133,  137,  123,  124,  125,  129,  130,
          //                134,  144,  145,  146,  147,  151,  152,  153,  148,  149,  150,
          //                154,  155,  156,  160,  161,  162,  157,  158,  159,  163,  164,
          //                168,  169,  170,  165,  166,  167,  171,  175,  176,  177,  172,
          //                173,  174,  193,  194,  195,  181,  186,  190,  196,  197,  182,
          //                187,  191,  198,  183,  188,  192,  178,  179,  180,  184,  185,
          //                189,  199,  200,  201,  205,  206,  207,  202,  203,  204,  208,
          //                209,  213,  214,  215,  210,  211,  212,  216,  220,  221,  222,
          //                217,  218,  219,  238,  239,  240,  226,  231,  235,  241,  242,
          //                227,  232,  236,  243,  228,  233,  237,  223,  224,  225,  229,
          //                230,  234,  244,  245,  249,  250,  251,  246,  247,  248,  252,
          //                256,  257,  258,  253,  254,  255,  274,  275,  276,  262,  267,
          //                271,  277,  278,  263,  268,  272,  279,  264,  269,  273,  259,
          //                260,  261,  265,  266,  270,  280,  284,  285,  286,  281,  282,
          //                283,  302,  303,  304,  290,  295,  299,  305,  306,  291,  296,
          //                300,  307,  292,  297,  301,  287,  288,  289,  293,  294,  298,
          //                354,  355,  356,  323,  338,  348,  357,  358,  324,  339,  349,
          //                359,  325,  340,  350,  311,  316,  320,  331,  335,  345,  360,
          //                361,  326,  341,  351,  362,  327,  342,  352,  312,  317,  321,
          //                332,  336,  346,  363,  328,  343,  353,  313,  318,  322,  333,
          //                337,  347,  308,  309,  310,  314,  315,  319,  329,  330,  334,
          //                344,  364,  365,  366,  367,  368,  372,  373,  374,  369,  370,
          //                371,  375,  376,  377,  378,  382,  383,  384,  379,  380,  381,
          //                385,  386,  387,  391,  392,  393,  388,  389,  390,  394,  395,
          //                399,  400,  401,  396,  397,  398,  402,  406,  407,  408,  403,
          //                404,  405,  424,  425,  426,  412,  417,  421,  427,  428,  413,
          //                418,  422,  429,  414,  419,  423,  409,  410,  411,  415,  416,
          //                420,  430,  431,  432,  433,  437,  438,  439,  434,  435,  436,
          //                440,  441,  442,  446,  447,  448,  443,  444,  445,  449,  450,
          //                454,  455,  456,  451,  452,  453,  457,  461,  462,  463,  458,
          //                459,  460,  479,  480,  481,  467,  472,  476,  482,  483,  468,
          //                473,  477,  484,  469,  474,  478,  464,  465,  466,  470,  471,
          //                475,  485,  486,  487,  491,  492,  493,  488,  489,  490,  494,
          //                495,  499,  500,  501,  496,  497,  498,  502,  506,  507,  508,
          //                503,  504,  505,  524,  525,  526,  512,  517,  521,  527,  528,
          //                513,  518,  522,  529,  514,  519,  523,  509,  510,  511,  515,
          //                516,  520,  530,  531,  535,  536,  537,  532,  533,  534,  538,
          //                542,  543,  544,  539,  540,  541,  560,  561,  562,  548,  553,
          //                557,  563,  564,  549,  554,  558,  565,  550,  555,  559,  545,
          //                546,  547,  551,  552,  556,  566,  570,  571,  572,  567,  568,
          //                569,  588,  589,  590,  576,  581,  585,  591,  592,  577,  582,
          //                586,  593,  578,  583,  587,  573,  574,  575,  579,  580,  584,
          //                640,  641,  642,  609,  624,  634,  643,  644,  610,  625,  635,
          //                645,  611,  626,  636,  597,  602,  606,  617,  621,  631,  646,
          //                647,  612,  627,  637,  648,  613,  628,  638,  598,  603,  607,
          //                618,  622,  632,  649,  614,  629,  639,  599,  604,  608,  619,
          //                623,  633,  594,  595,  596,  600,  601,  605,  615,  616,  620,
          //                630,  650,  651,  652,  653,  657,  658,  659,  654,  655,  656,
          //                660,  661,  662,  666,  667,  668,  663,  664,  665,  669,  670,
          //                674,  675,  676,  671,  672,  673,  677,  681,  682,  683,  678,
          //                679,  680,  699,  700,  701,  687,  692,  696,  702,  703,  688,
          //                693,  697,  704,  689,  694,  698,  684,  685,  686,  690,  691,
          //                695,  705,  706,  707,  711,  712,  713,  708,  709,  710,  714,
          //                715,  719,  720,  721,  716,  717,  718,  722,  726,  727,  728,
          //                723,  724,  725,  744,  745,  746,  732,  737,  741,  747,  748,
          //                733,  738,  742,  749,  734,  739,  743,  729,  730,  731,  735,
          //                736,  740,  750,  751,  755,  756,  757,  752,  753,  754,  758,
          //                762,  763,  764,  759,  760,  761,  780,  781,  782,  768,  773,
          //                777,  783,  784,  769,  774,  778,  785,  770,  775,  779,  765,
          //                766,  767,  771,  772,  776,  786,  790,  791,  792,  787,  788,
          //                789,  808,  809,  810,  796,  801,  805,  811,  812,  797,  802,
          //                806,  813,  798,  803,  807,  793,  794,  795,  799,  800,  804,
          //                860,  861,  862,  829,  844,  854,  863,  864,  830,  845,  855,
          //                865,  831,  846,  856,  817,  822,  826,  837,  841,  851,  866,
          //                867,  832,  847,  857,  868,  833,  848,  858,  818,  823,  827,
          //                838,  842,  852,  869,  834,  849,  859,  819,  824,  828,  839,
          //                843,  853,  814,  815,  816,  820,  821,  825,  835,  836,  840,
          //                850,  870,  871,  872,  876,  877,  878,  873,  874,  875,  879,
          //                880,  884,  885,  886,  881,  882,  883,  887,  891,  892,  893,
          //                888,  889,  890,  909,  910,  911,  897,  902,  906,  912,  913,
          //                898,  903,  907,  914,  899,  904,  908,  894,  895,  896,  900,
          //                901,  905,  915,  916,  920,  921,  922,  917,  918,  919,  923,
          //                927,  928,  929,  924,  925,  926,  945,  946,  947,  933,  938,
          //                942,  948,  949,  934,  939,  943,  950,  935,  940,  944,  930,
          //                931,  932,  936,  937,  941,  951,  955,  956,  957,  952,  953,
          //                954,  973,  974,  975,  961,  966,  970,  976,  977,  962,  967,
          //                971,  978,  963,  968,  972,  958,  959,  960,  964,  965,  969,
          //               1025, 1026, 1027,  994, 1009, 1019, 1028, 1029,  995, 1010, 1020,
          //               1030,  996, 1011, 1021,  982,  987,  991, 1002, 1006, 1016, 1031,
          //               1032,  997, 1012, 1022, 1033,  998, 1013, 1023,  983,  988,  992,
          //               1003, 1007, 1017, 1034,  999, 1014, 1024,  984,  989,  993, 1004,
          //               1008, 1018,  979,  980,  981,  985,  986,  990, 1000, 1001, 1005,
          //               1015, 1035, 1036, 1040, 1041, 1042, 1037, 1038, 1039, 1043, 1047,
          //               1048, 1049, 1044, 1045, 1046, 1065, 1066, 1067, 1053, 1058, 1062,
          //               1068, 1069, 1054, 1059, 1063, 1070, 1055, 1060, 1064, 1050, 1051,
          //               1052, 1056, 1057, 1061, 1071, 1075, 1076, 1077, 1072, 1073, 1074,
          //               1093, 1094, 1095, 1081, 1086, 1090, 1096, 1097, 1082, 1087, 1091,
          //               1098, 1083, 1088, 1092, 1078, 1079, 1080, 1084, 1085, 1089, 1145,
          //               1146, 1147, 1114, 1129, 1139, 1148, 1149, 1115, 1130, 1140, 1150,
          //               1116, 1131, 1141, 1102, 1107, 1111, 1122, 1126, 1136, 1151, 1152,
          //               1117, 1132, 1142, 1153, 1118, 1133, 1143, 1103, 1108, 1112, 1123,
          //               1127, 1137, 1154, 1119, 1134, 1144, 1104, 1109, 1113, 1124, 1128,
          //               1138, 1099, 1100, 1101, 1105, 1106, 1110, 1120, 1121, 1125, 1135,
          //               1155, 1159, 1160, 1161, 1156, 1157, 1158, 1177, 1178, 1179, 1165,
          //               1170, 1174, 1180, 1181, 1166, 1171, 1175, 1182, 1167, 1172, 1176,
          //               1162, 1163, 1164, 1168, 1169, 1173, 1229, 1230, 1231, 1198, 1213,
          //               1223, 1232, 1233, 1199, 1214, 1224, 1234, 1200, 1215, 1225, 1186,
          //               1191, 1195, 1206, 1210, 1220, 1235, 1236, 1201, 1216, 1226, 1237,
          //               1202, 1217, 1227, 1187, 1192, 1196, 1207, 1211, 1221, 1238, 1203,
          //               1218, 1228, 1188, 1193, 1197, 1208, 1212, 1222, 1183, 1184, 1185,
          //               1189, 1190, 1194, 1204, 1205, 1209, 1219, 1350, 1351, 1352, 1285,
          //               1320, 1340, 1353, 1354, 1286, 1321, 1341, 1355, 1287, 1322, 1342,
          //               1254, 1269, 1279, 1304, 1314, 1334, 1356, 1357, 1288, 1323, 1343,
          //               1358, 1289, 1324, 1344, 1255, 1270, 1280, 1305, 1315, 1335, 1359,
          //               1290, 1325, 1345, 1256, 1271, 1281, 1306, 1316, 1336, 1242, 1247,
          //               1251, 1262, 1266, 1276, 1297, 1301, 1311, 1331, 1360, 1361, 1291,
          //               1326, 1346, 1362, 1292, 1327, 1347, 1257, 1272, 1282, 1307, 1317,
          //               1337, 1363, 1293, 1328, 1348, 1258, 1273, 1283, 1308, 1318, 1338,
          //               1243, 1248, 1252, 1263, 1267, 1277, 1298, 1302, 1312, 1332, 1364,
          //               1294, 1329, 1349, 1259, 1274, 1284, 1309, 1319, 1339, 1244, 1249,
          //               1253, 1264, 1268, 1278, 1299, 1303, 1313, 1333, 1239, 1240, 1241,
          //               1245, 1246, 1250, 1260, 1261, 1265, 1275, 1295, 1296, 1300, 1310,
          //               1330}},
          //             {{ 870,  871,  872,  199,  485,  705,  873,  874,  875,  876,  877,
          //                878,  879,  880,  200,  486,  706,  881,  882,  883,  884,  885,
          //                886,  887,  201,  487,  707,  888,  889,  890,  891,  892,  893,
          //                 33,   99,  154,  202,  203,  204,  205,  206,  207,  385,  440,
          //                488,  489,  490,  491,  492,  493,  660,  708,  709,  710,  711,
          //                712,  713,  894,  895,  896,  897,  898,  899,  900,  901,  902,
          //                903,  904,  905,  906,  907,  908,  909,  910,  911,  912,  913,
          //                914,  915,  916,  208,  494,  714,  917,  918,  919,  920,  921,
          //                922,  923,  209,  495,  715,  924,  925,  926,  927,  928,  929,
          //                 34,  100,  155,  210,  211,  212,  213,  214,  215,  386,  441,
          //                496,  497,  498,  499,  500,  501,  661,  716,  717,  718,  719,
          //                720,  721,  930,  931,  932,  933,  934,  935,  936,  937,  938,
          //                939,  940,  941,  942,  943,  944,  945,  946,  947,  948,  949,
          //                950,  951,  216,  502,  722,  952,  953,  954,  955,  956,  957,
          //                 35,  101,  156,  217,  218,  219,  220,  221,  222,  387,  442,
          //                503,  504,  505,  506,  507,  508,  662,  723,  724,  725,  726,
          //                727,  728,  958,  959,  960,  961,  962,  963,  964,  965,  966,
          //                967,  968,  969,  970,  971,  972,  973,  974,  975,  976,  977,
          //                978,    3,   14,   24,   36,   37,   38,   39,   40,   41,   80,
          //                 90,  102,  103,  104,  105,  106,  107,  145,  157,  158,  159,
          //                160,  161,  162,  223,  224,  225,  226,  227,  228,  229,  230,
          //                231,  232,  233,  234,  235,  236,  237,  238,  239,  240,  241,
          //                242,  243,  366,  376,  388,  389,  390,  391,  392,  393,  431,
          //                443,  444,  445,  446,  447,  448,  509,  510,  511,  512,  513,
          //                514,  515,  516,  517,  518,  519,  520,  521,  522,  523,  524,
          //                525,  526,  527,  528,  529,  651,  663,  664,  665,  666,  667,
          //                668,  729,  730,  731,  732,  733,  734,  735,  736,  737,  738,
          //                739,  740,  741,  742,  743,  744,  745,  746,  747,  748,  749,
          //                979,  980,  981,  982,  983,  984,  985,  986,  987,  988,  989,
          //                990,  991,  992,  993,  994,  995,  996,  997,  998,  999, 1000,
          //               1001, 1002, 1003, 1004, 1005, 1006, 1007, 1008, 1009, 1010, 1011,
          //               1012, 1013, 1014, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022,
          //               1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033,
          //               1034, 1035, 1036,  244,  530,  750, 1037, 1038, 1039, 1040, 1041,
          //               1042, 1043,  245,  531,  751, 1044, 1045, 1046, 1047, 1048, 1049,
          //                 42,  108,  163,  246,  247,  248,  249,  250,  251,  394,  449,
          //                532,  533,  534,  535,  536,  537,  669,  752,  753,  754,  755,
          //                756,  757, 1050, 1051, 1052, 1053, 1054, 1055, 1056, 1057, 1058,
          //               1059, 1060, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1068, 1069,
          //               1070, 1071,  252,  538,  758, 1072, 1073, 1074, 1075, 1076, 1077,
          //                 43,  109,  164,  253,  254,  255,  256,  257,  258,  395,  450,
          //                539,  540,  541,  542,  543,  544,  670,  759,  760,  761,  762,
          //                763,  764, 1078, 1079, 1080, 1081, 1082, 1083, 1084, 1085, 1086,
          //               1087, 1088, 1089, 1090, 1091, 1092, 1093, 1094, 1095, 1096, 1097,
          //               1098,    4,   15,   25,   44,   45,   46,   47,   48,   49,   81,
          //                 91,  110,  111,  112,  113,  114,  115,  146,  165,  166,  167,
          //                168,  169,  170,  259,  260,  261,  262,  263,  264,  265,  266,
          //                267,  268,  269,  270,  271,  272,  273,  274,  275,  276,  277,
          //                278,  279,  367,  377,  396,  397,  398,  399,  400,  401,  432,
          //                451,  452,  453,  454,  455,  456,  545,  546,  547,  548,  549,
          //                550,  551,  552,  553,  554,  555,  556,  557,  558,  559,  560,
          //                561,  562,  563,  564,  565,  652,  671,  672,  673,  674,  675,
          //                676,  765,  766,  767,  768,  769,  770,  771,  772,  773,  774,
          //                775,  776,  777,  778,  779,  780,  781,  782,  783,  784,  785,
          //               1099, 1100, 1101, 1102, 1103, 1104, 1105, 1106, 1107, 1108, 1109,
          //               1110, 1111, 1112, 1113, 1114, 1115, 1116, 1117, 1118, 1119, 1120,
          //               1121, 1122, 1123, 1124, 1125, 1126, 1127, 1128, 1129, 1130, 1131,
          //               1132, 1133, 1134, 1135, 1136, 1137, 1138, 1139, 1140, 1141, 1142,
          //               1143, 1144, 1145, 1146, 1147, 1148, 1149, 1150, 1151, 1152, 1153,
          //               1154, 1155,  280,  566,  786, 1156, 1157, 1158, 1159, 1160, 1161,
          //                 50,  116,  171,  281,  282,  283,  284,  285,  286,  402,  457,
          //                567,  568,  569,  570,  571,  572,  677,  787,  788,  789,  790,
          //                791,  792, 1162, 1163, 1164, 1165, 1166, 1167, 1168, 1169, 1170,
          //               1171, 1172, 1173, 1174, 1175, 1176, 1177, 1178, 1179, 1180, 1181,
          //               1182,    5,   16,   26,   51,   52,   53,   54,   55,   56,   82,
          //                 92,  117,  118,  119,  120,  121,  122,  147,  172,  173,  174,
          //                175,  176,  177,  287,  288,  289,  290,  291,  292,  293,  294,
          //                295,  296,  297,  298,  299,  300,  301,  302,  303,  304,  305,
          //                306,  307,  368,  378,  403,  404,  405,  406,  407,  408,  433,
          //                458,  459,  460,  461,  462,  463,  573,  574,  575,  576,  577,
          //                578,  579,  580,  581,  582,  583,  584,  585,  586,  587,  588,
          //                589,  590,  591,  592,  593,  653,  678,  679,  680,  681,  682,
          //                683,  793,  794,  795,  796,  797,  798,  799,  800,  801,  802,
          //                803,  804,  805,  806,  807,  808,  809,  810,  811,  812,  813,
          //               1183, 1184, 1185, 1186, 1187, 1188, 1189, 1190, 1191, 1192, 1193,
          //               1194, 1195, 1196, 1197, 1198, 1199, 1200, 1201, 1202, 1203, 1204,
          //               1205, 1206, 1207, 1208, 1209, 1210, 1211, 1212, 1213, 1214, 1215,
          //               1216, 1217, 1218, 1219, 1220, 1221, 1222, 1223, 1224, 1225, 1226,
          //               1227, 1228, 1229, 1230, 1231, 1232, 1233, 1234, 1235, 1236, 1237,
          //               1238,    0,    1,    2,    6,    7,    8,    9,   10,   11,   12,
          //                 13,   17,   18,   19,   20,   21,   22,   23,   27,   28,   29,
          //                 30,   31,   32,   57,   58,   59,   60,   61,   62,   63,   64,
          //                 65,   66,   67,   68,   69,   70,   71,   72,   73,   74,   75,
          //                 76,   77,   78,   79,   83,   84,   85,   86,   87,   88,   89,
          //                 93,   94,   95,   96,   97,   98,  123,  124,  125,  126,  127,
          //                128,  129,  130,  131,  132,  133,  134,  135,  136,  137,  138,
          //                139,  140,  141,  142,  143,  144,  148,  149,  150,  151,  152,
          //                153,  178,  179,  180,  181,  182,  183,  184,  185,  186,  187,
          //                188,  189,  190,  191,  192,  193,  194,  195,  196,  197,  198,
          //                308,  309,  310,  311,  312,  313,  314,  315,  316,  317,  318,
          //                319,  320,  321,  322,  323,  324,  325,  326,  327,  328,  329,
          //                330,  331,  332,  333,  334,  335,  336,  337,  338,  339,  340,
          //                341,  342,  343,  344,  345,  346,  347,  348,  349,  350,  351,
          //                352,  353,  354,  355,  356,  357,  358,  359,  360,  361,  362,
          //                363,  364,  365,  369,  370,  371,  372,  373,  374,  375,  379,
          //                380,  381,  382,  383,  384,  409,  410,  411,  412,  413,  414,
          //                415,  416,  417,  418,  419,  420,  421,  422,  423,  424,  425,
          //                426,  427,  428,  429,  430,  434,  435,  436,  437,  438,  439,
          //                464,  465,  466,  467,  468,  469,  470,  471,  472,  473,  474,
          //                475,  476,  477,  478,  479,  480,  481,  482,  483,  484,  594,
          //                595,  596,  597,  598,  599,  600,  601,  602,  603,  604,  605,
          //                606,  607,  608,  609,  610,  611,  612,  613,  614,  615,  616,
          //                617,  618,  619,  620,  621,  622,  623,  624,  625,  626,  627,
          //                628,  629,  630,  631,  632,  633,  634,  635,  636,  637,  638,
          //                639,  640,  641,  642,  643,  644,  645,  646,  647,  648,  649,
          //                650,  654,  655,  656,  657,  658,  659,  684,  685,  686,  687,
          //                688,  689,  690,  691,  692,  693,  694,  695,  696,  697,  698,
          //                699,  700,  701,  702,  703,  704,  814,  815,  816,  817,  818,
          //                819,  820,  821,  822,  823,  824,  825,  826,  827,  828,  829,
          //                830,  831,  832,  833,  834,  835,  836,  837,  838,  839,  840,
          //                841,  842,  843,  844,  845,  846,  847,  848,  849,  850,  851,
          //                852,  853,  854,  855,  856,  857,  858,  859,  860,  861,  862,
          //                863,  864,  865,  866,  867,  868,  869, 1239, 1240, 1241, 1242,
          //               1243, 1244, 1245, 1246, 1247, 1248, 1249, 1250, 1251, 1252, 1253,
          //               1254, 1255, 1256, 1257, 1258, 1259, 1260, 1261, 1262, 1263, 1264,
          //               1265, 1266, 1267, 1268, 1269, 1270, 1271, 1272, 1273, 1274, 1275,
          //               1276, 1277, 1278, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286,
          //               1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297,
          //               1298, 1299, 1300, 1301, 1302, 1303, 1304, 1305, 1306, 1307, 1308,
          //               1309, 1310, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319,
          //               1320, 1321, 1322, 1323, 1324, 1325, 1326, 1327, 1328, 1329, 1330,
          //               1331, 1332, 1333, 1334, 1335, 1336, 1337, 1338, 1339, 1340, 1341,
          //               1342, 1343, 1344, 1345, 1346, 1347, 1348, 1349, 1350, 1351, 1352,
          //               1353, 1354, 1355, 1356, 1357, 1358, 1359, 1360, 1361, 1362, 1363,
          //               1364},
          //              { 870,  871,  872,  199,  485,  705,  876,  877,  878,  873,  874,
          //                875,  879,  880,  200,  486,  706,  884,  885,  886,  881,  882,
          //                883,  887,  201,  487,  707,  891,  892,  893,  888,  889,  890,
          //                 33,   99,  154,  205,  206,  207,  202,  203,  204,  385,  440,
          //                491,  492,  493,  488,  489,  490,  660,  711,  712,  713,  708,
          //                709,  710,  909,  910,  911,  897,  902,  906,  912,  913,  898,
          //                903,  907,  914,  899,  904,  908,  894,  895,  896,  900,  901,
          //                905,  915,  916,  208,  494,  714,  920,  921,  922,  917,  918,
          //                919,  923,  209,  495,  715,  927,  928,  929,  924,  925,  926,
          //                 34,  100,  155,  213,  214,  215,  210,  211,  212,  386,  441,
          //                499,  500,  501,  496,  497,  498,  661,  719,  720,  721,  716,
          //                717,  718,  945,  946,  947,  933,  938,  942,  948,  949,  934,
          //                939,  943,  950,  935,  940,  944,  930,  931,  932,  936,  937,
          //                941,  951,  216,  502,  722,  955,  956,  957,  952,  953,  954,
          //                 35,  101,  156,  220,  221,  222,  217,  218,  219,  387,  442,
          //                506,  507,  508,  503,  504,  505,  662,  726,  727,  728,  723,
          //                724,  725,  973,  974,  975,  961,  966,  970,  976,  977,  962,
          //                967,  971,  978,  963,  968,  972,  958,  959,  960,  964,  965,
          //                969,    3,   14,   24,   39,   40,   41,   36,   37,   38,   80,
          //                 90,  105,  106,  107,  102,  103,  104,  145,  160,  161,  162,
          //                157,  158,  159,  238,  239,  240,  226,  231,  235,  241,  242,
          //                227,  232,  236,  243,  228,  233,  237,  223,  224,  225,  229,
          //                230,  234,  366,  376,  391,  392,  393,  388,  389,  390,  431,
          //                446,  447,  448,  443,  444,  445,  524,  525,  526,  512,  517,
          //                521,  527,  528,  513,  518,  522,  529,  514,  519,  523,  509,
          //                510,  511,  515,  516,  520,  651,  666,  667,  668,  663,  664,
          //                665,  744,  745,  746,  732,  737,  741,  747,  748,  733,  738,
          //                742,  749,  734,  739,  743,  729,  730,  731,  735,  736,  740,
          //               1025, 1026, 1027,  994, 1009, 1019, 1028, 1029,  995, 1010, 1020,
          //               1030,  996, 1011, 1021,  982,  987,  991, 1002, 1006, 1016, 1031,
          //               1032,  997, 1012, 1022, 1033,  998, 1013, 1023,  983,  988,  992,
          //               1003, 1007, 1017, 1034,  999, 1014, 1024,  984,  989,  993, 1004,
          //               1008, 1018,  979,  980,  981,  985,  986,  990, 1000, 1001, 1005,
          //               1015, 1035, 1036,  244,  530,  750, 1040, 1041, 1042, 1037, 1038,
          //               1039, 1043,  245,  531,  751, 1047, 1048, 1049, 1044, 1045, 1046,
          //                 42,  108,  163,  249,  250,  251,  246,  247,  248,  394,  449,
          //                535,  536,  537,  532,  533,  534,  669,  755,  756,  757,  752,
          //                753,  754, 1065, 1066, 1067, 1053, 1058, 1062, 1068, 1069, 1054,
          //               1059, 1063, 1070, 1055, 1060, 1064, 1050, 1051, 1052, 1056, 1057,
          //               1061, 1071,  252,  538,  758, 1075, 1076, 1077, 1072, 1073, 1074,
          //                 43,  109,  164,  256,  257,  258,  253,  254,  255,  395,  450,
          //                542,  543,  544,  539,  540,  541,  670,  762,  763,  764,  759,
          //                760,  761, 1093, 1094, 1095, 1081, 1086, 1090, 1096, 1097, 1082,
          //               1087, 1091, 1098, 1083, 1088, 1092, 1078, 1079, 1080, 1084, 1085,
          //               1089,    4,   15,   25,   47,   48,   49,   44,   45,   46,   81,
          //                 91,  113,  114,  115,  110,  111,  112,  146,  168,  169,  170,
          //                165,  166,  167,  274,  275,  276,  262,  267,  271,  277,  278,
          //                263,  268,  272,  279,  264,  269,  273,  259,  260,  261,  265,
          //                266,  270,  367,  377,  399,  400,  401,  396,  397,  398,  432,
          //                454,  455,  456,  451,  452,  453,  560,  561,  562,  548,  553,
          //                557,  563,  564,  549,  554,  558,  565,  550,  555,  559,  545,
          //                546,  547,  551,  552,  556,  652,  674,  675,  676,  671,  672,
          //                673,  780,  781,  782,  768,  773,  777,  783,  784,  769,  774,
          //                778,  785,  770,  775,  779,  765,  766,  767,  771,  772,  776,
          //               1145, 1146, 1147, 1114, 1129, 1139, 1148, 1149, 1115, 1130, 1140,
          //               1150, 1116, 1131, 1141, 1102, 1107, 1111, 1122, 1126, 1136, 1151,
          //               1152, 1117, 1132, 1142, 1153, 1118, 1133, 1143, 1103, 1108, 1112,
          //               1123, 1127, 1137, 1154, 1119, 1134, 1144, 1104, 1109, 1113, 1124,
          //               1128, 1138, 1099, 1100, 1101, 1105, 1106, 1110, 1120, 1121, 1125,
          //               1135, 1155,  280,  566,  786, 1159, 1160, 1161, 1156, 1157, 1158,
          //                 50,  116,  171,  284,  285,  286,  281,  282,  283,  402,  457,
          //                570,  571,  572,  567,  568,  569,  677,  790,  791,  792,  787,
          //                788,  789, 1177, 1178, 1179, 1165, 1170, 1174, 1180, 1181, 1166,
          //               1171, 1175, 1182, 1167, 1172, 1176, 1162, 1163, 1164, 1168, 1169,
          //               1173,    5,   16,   26,   54,   55,   56,   51,   52,   53,   82,
          //                 92,  120,  121,  122,  117,  118,  119,  147,  175,  176,  177,
          //                172,  173,  174,  302,  303,  304,  290,  295,  299,  305,  306,
          //                291,  296,  300,  307,  292,  297,  301,  287,  288,  289,  293,
          //                294,  298,  368,  378,  406,  407,  408,  403,  404,  405,  433,
          //                461,  462,  463,  458,  459,  460,  588,  589,  590,  576,  581,
          //                585,  591,  592,  577,  582,  586,  593,  578,  583,  587,  573,
          //                574,  575,  579,  580,  584,  653,  681,  682,  683,  678,  679,
          //                680,  808,  809,  810,  796,  801,  805,  811,  812,  797,  802,
          //                806,  813,  798,  803,  807,  793,  794,  795,  799,  800,  804,
          //               1229, 1230, 1231, 1198, 1213, 1223, 1232, 1233, 1199, 1214, 1224,
          //               1234, 1200, 1215, 1225, 1186, 1191, 1195, 1206, 1210, 1220, 1235,
          //               1236, 1201, 1216, 1226, 1237, 1202, 1217, 1227, 1187, 1192, 1196,
          //               1207, 1211, 1221, 1238, 1203, 1218, 1228, 1188, 1193, 1197, 1208,
          //               1212, 1222, 1183, 1184, 1185, 1189, 1190, 1194, 1204, 1205, 1209,
          //               1219,    0,    1,    2,    9,   10,   11,    6,    7,    8,   12,
          //                 13,   20,   21,   22,   17,   18,   19,   23,   30,   31,   32,
          //                 27,   28,   29,   72,   73,   74,   60,   65,   69,   75,   76,
          //                 61,   66,   70,   77,   62,   67,   71,   57,   58,   59,   63,
          //                 64,   68,   78,   79,   86,   87,   88,   83,   84,   85,   89,
          //                 96,   97,   98,   93,   94,   95,  138,  139,  140,  126,  131,
          //                135,  141,  142,  127,  132,  136,  143,  128,  133,  137,  123,
          //                124,  125,  129,  130,  134,  144,  151,  152,  153,  148,  149,
          //                150,  193,  194,  195,  181,  186,  190,  196,  197,  182,  187,
          //                191,  198,  183,  188,  192,  178,  179,  180,  184,  185,  189,
          //                354,  355,  356,  323,  338,  348,  357,  358,  324,  339,  349,
          //                359,  325,  340,  350,  311,  316,  320,  331,  335,  345,  360,
          //                361,  326,  341,  351,  362,  327,  342,  352,  312,  317,  321,
          //                332,  336,  346,  363,  328,  343,  353,  313,  318,  322,  333,
          //                337,  347,  308,  309,  310,  314,  315,  319,  329,  330,  334,
          //                344,  364,  365,  372,  373,  374,  369,  370,  371,  375,  382,
          //                383,  384,  379,  380,  381,  424,  425,  426,  412,  417,  421,
          //                427,  428,  413,  418,  422,  429,  414,  419,  423,  409,  410,
          //                411,  415,  416,  420,  430,  437,  438,  439,  434,  435,  436,
          //                479,  480,  481,  467,  472,  476,  482,  483,  468,  473,  477,
          //                484,  469,  474,  478,  464,  465,  466,  470,  471,  475,  640,
          //                641,  642,  609,  624,  634,  643,  644,  610,  625,  635,  645,
          //                611,  626,  636,  597,  602,  606,  617,  621,  631,  646,  647,
          //                612,  627,  637,  648,  613,  628,  638,  598,  603,  607,  618,
          //                622,  632,  649,  614,  629,  639,  599,  604,  608,  619,  623,
          //                633,  594,  595,  596,  600,  601,  605,  615,  616,  620,  630,
          //                650,  657,  658,  659,  654,  655,  656,  699,  700,  701,  687,
          //                692,  696,  702,  703,  688,  693,  697,  704,  689,  694,  698,
          //                684,  685,  686,  690,  691,  695,  860,  861,  862,  829,  844,
          //                854,  863,  864,  830,  845,  855,  865,  831,  846,  856,  817,
          //                822,  826,  837,  841,  851,  866,  867,  832,  847,  857,  868,
          //                833,  848,  858,  818,  823,  827,  838,  842,  852,  869,  834,
          //                849,  859,  819,  824,  828,  839,  843,  853,  814,  815,  816,
          //                820,  821,  825,  835,  836,  840,  850, 1350, 1351, 1352, 1285,
          //               1320, 1340, 1353, 1354, 1286, 1321, 1341, 1355, 1287, 1322, 1342,
          //               1254, 1269, 1279, 1304, 1314, 1334, 1356, 1357, 1288, 1323, 1343,
          //               1358, 1289, 1324, 1344, 1255, 1270, 1280, 1305, 1315, 1335, 1359,
          //               1290, 1325, 1345, 1256, 1271, 1281, 1306, 1316, 1336, 1242, 1247,
          //               1251, 1262, 1266, 1276, 1297, 1301, 1311, 1331, 1360, 1361, 1291,
          //               1326, 1346, 1362, 1292, 1327, 1347, 1257, 1272, 1282, 1307, 1317,
          //               1337, 1363, 1293, 1328, 1348, 1258, 1273, 1283, 1308, 1318, 1338,
          //               1243, 1248, 1252, 1263, 1267, 1277, 1298, 1302, 1312, 1332, 1364,
          //               1294, 1329, 1349, 1259, 1274, 1284, 1309, 1319, 1339, 1244, 1249,
          //               1253, 1264, 1268, 1278, 1299, 1303, 1313, 1333, 1239, 1240, 1241,
          //               1245, 1246, 1250, 1260, 1261, 1265, 1275, 1295, 1296, 1300, 1310,
          //               1330}}},
          //            {{{1239, 1240, 1241, 1242, 1243, 1244,  308,  594,  814,  979, 1099,
          //               1183, 1245, 1246, 1247, 1248, 1249,  309,  595,  815,  980, 1100,
          //               1184, 1250, 1251, 1252, 1253,  310,  596,  816,  981, 1101, 1185,
          //               1254, 1255, 1256,  311,  597,  817,  982, 1102, 1186, 1257, 1258,
          //                312,  598,  818,  983, 1103, 1187, 1259,  313,  599,  819,  984,
          //               1104, 1188,   57,  123,  178,  223,  259,  287,  409,  464,  509,
          //                545,  573,  684,  729,  765,  793,  894,  930,  958, 1050, 1078,
          //               1162, 1260, 1261, 1262, 1263, 1264,  314,  600,  820,  985, 1105,
          //               1189, 1265, 1266, 1267, 1268,  315,  601,  821,  986, 1106, 1190,
          //               1269, 1270, 1271,  316,  602,  822,  987, 1107, 1191, 1272, 1273,
          //                317,  603,  823,  988, 1108, 1192, 1274,  318,  604,  824,  989,
          //               1109, 1193,   58,  124,  179,  224,  260,  288,  410,  465,  510,
          //                546,  574,  685,  730,  766,  794,  895,  931,  959, 1051, 1079,
          //               1163, 1275, 1276, 1277, 1278,  319,  605,  825,  990, 1110, 1194,
          //               1279, 1280, 1281,  320,  606,  826,  991, 1111, 1195, 1282, 1283,
          //                321,  607,  827,  992, 1112, 1196, 1284,  322,  608,  828,  993,
          //               1113, 1197,   59,  125,  180,  225,  261,  289,  411,  466,  511,
          //                547,  575,  686,  731,  767,  795,  896,  932,  960, 1052, 1080,
          //               1164, 1285, 1286, 1287,  323,  609,  829,  994, 1114, 1198, 1288,
          //               1289,  324,  610,  830,  995, 1115, 1199, 1290,  325,  611,  831,
          //                996, 1116, 1200,   60,  126,  181,  226,  262,  290,  412,  467,
          //                512,  548,  576,  687,  732,  768,  796,  897,  933,  961, 1053,
          //               1081, 1165, 1291, 1292,  326,  612,  832,  997, 1117, 1201, 1293,
          //                327,  613,  833,  998, 1118, 1202,   61,  127,  182,  227,  263,
          //                291,  413,  468,  513,  549,  577,  688,  733,  769,  797,  898,
          //                934,  962, 1054, 1082, 1166, 1294,  328,  614,  834,  999, 1119,
          //               1203,   62,  128,  183,  228,  264,  292,  414,  469,  514,  550,
          //                578,  689,  734,  770,  798,  899,  935,  963, 1055, 1083, 1167,
          //                  6,   17,   27,   36,   44,   51,   83,   93,  102,  110,  117,
          //                148,  157,  165,  172,  202,  210,  217,  246,  253,  281,  369,
          //                379,  388,  396,  403,  434,  443,  451,  458,  488,  496,  503,
          //                532,  539,  567,  654,  663,  671,  678,  708,  716,  723,  752,
          //                759,  787,  873,  881,  888,  917,  924,  952, 1037, 1044, 1072,
          //               1156, 1295, 1296, 1297, 1298, 1299,  329,  615,  835, 1000, 1120,
          //               1204, 1300, 1301, 1302, 1303,  330,  616,  836, 1001, 1121, 1205,
          //               1304, 1305, 1306,  331,  617,  837, 1002, 1122, 1206, 1307, 1308,
          //                332,  618,  838, 1003, 1123, 1207, 1309,  333,  619,  839, 1004,
          //               1124, 1208,   63,  129,  184,  229,  265,  293,  415,  470,  515,
          //                551,  579,  690,  735,  771,  799,  900,  936,  964, 1056, 1084,
          //               1168, 1310, 1311, 1312, 1313,  334,  620,  840, 1005, 1125, 1209,
          //               1314, 1315, 1316,  335,  621,  841, 1006, 1126, 1210, 1317, 1318,
          //                336,  622,  842, 1007, 1127, 1211, 1319,  337,  623,  843, 1008,
          //               1128, 1212,   64,  130,  185,  230,  266,  294,  416,  471,  516,
          //                552,  580,  691,  736,  772,  800,  901,  937,  965, 1057, 1085,
          //               1169, 1320, 1321, 1322,  338,  624,  844, 1009, 1129, 1213, 1323,
          //               1324,  339,  625,  845, 1010, 1130, 1214, 1325,  340,  626,  846,
          //               1011, 1131, 1215,   65,  131,  186,  231,  267,  295,  417,  472,
          //                517,  553,  581,  692,  737,  773,  801,  902,  938,  966, 1058,
          //               1086, 1170, 1326, 1327,  341,  627,  847, 1012, 1132, 1216, 1328,
          //                342,  628,  848, 1013, 1133, 1217,   66,  132,  187,  232,  268,
          //                296,  418,  473,  518,  554,  582,  693,  738,  774,  802,  903,
          //                939,  967, 1059, 1087, 1171, 1329,  343,  629,  849, 1014, 1134,
          //               1218,   67,  133,  188,  233,  269,  297,  419,  474,  519,  555,
          //                583,  694,  739,  775,  803,  904,  940,  968, 1060, 1088, 1172,
          //                  7,   18,   28,   37,   45,   52,   84,   94,  103,  111,  118,
          //                149,  158,  166,  173,  203,  211,  218,  247,  254,  282,  370,
          //                380,  389,  397,  404,  435,  444,  452,  459,  489,  497,  504,
          //                533,  540,  568,  655,  664,  672,  679,  709,  717,  724,  753,
          //                760,  788,  874,  882,  889,  918,  925,  953, 1038, 1045, 1073,
          //               1157, 1330, 1331, 1332, 1333,  344,  630,  850, 1015, 1135, 1219,
          //               1334, 1335, 1336,  345,  631,  851, 1016, 1136, 1220, 1337, 1338,
          //                346,  632,  852, 1017, 1137, 1221, 1339,  347,  633,  853, 1018,
          //               1138, 1222,   68,  134,  189,  234,  270,  298,  420,  475,  520,
          //                556,  584,  695,  740,  776,  804,  905,  941,  969, 1061, 1089,
          //               1173, 1340, 1341, 1342,  348,  634,  854, 1019, 1139, 1223, 1343,
          //               1344,  349,  635,  855, 1020, 1140, 1224, 1345,  350,  636,  856,
          //               1021, 1141, 1225,   69,  135,  190,  235,  271,  299,  421,  476,
          //                521,  557,  585,  696,  741,  777,  805,  906,  942,  970, 1062,
          //               1090, 1174, 1346, 1347,  351,  637,  857, 1022, 1142, 1226, 1348,
          //                352,  638,  858, 1023, 1143, 1227,   70,  136,  191,  236,  272,
          //                300,  422,  477,  522,  558,  586,  697,  742,  778,  806,  907,
          //                943,  971, 1063, 1091, 1175, 1349,  353,  639,  859, 1024, 1144,
          //               1228,   71,  137,  192,  237,  273,  301,  423,  478,  523,  559,
          //                587,  698,  743,  779,  807,  908,  944,  972, 1064, 1092, 1176,
          //                  8,   19,   29,   38,   46,   53,   85,   95,  104,  112,  119,
          //                150,  159,  167,  174,  204,  212,  219,  248,  255,  283,  371,
          //                381,  390,  398,  405,  436,  445,  453,  460,  490,  498,  505,
          //                534,  541,  569,  656,  665,  673,  680,  710,  718,  725,  754,
          //                761,  789,  875,  883,  890,  919,  926,  954, 1039, 1046, 1074,
          //               1158, 1350, 1351, 1352,  354,  640,  860, 1025, 1145, 1229, 1353,
          //               1354,  355,  641,  861, 1026, 1146, 1230, 1355,  356,  642,  862,
          //               1027, 1147, 1231,   72,  138,  193,  238,  274,  302,  424,  479,
          //                524,  560,  588,  699,  744,  780,  808,  909,  945,  973, 1065,
          //               1093, 1177, 1356, 1357,  357,  643,  863, 1028, 1148, 1232, 1358,
          //                358,  644,  864, 1029, 1149, 1233,   73,  139,  194,  239,  275,
          //                303,  425,  480,  525,  561,  589,  700,  745,  781,  809,  910,
          //                946,  974, 1066, 1094, 1178, 1359,  359,  645,  865, 1030, 1150,
          //               1234,   74,  140,  195,  240,  276,  304,  426,  481,  526,  562,
          //                590,  701,  746,  782,  810,  911,  947,  975, 1067, 1095, 1179,
          //                  9,   20,   30,   39,   47,   54,   86,   96,  105,  113,  120,
          //                151,  160,  168,  175,  205,  213,  220,  249,  256,  284,  372,
          //                382,  391,  399,  406,  437,  446,  454,  461,  491,  499,  506,
          //                535,  542,  570,  657,  666,  674,  681,  711,  719,  726,  755,
          //                762,  790,  876,  884,  891,  920,  927,  955, 1040, 1047, 1075,
          //               1159, 1360, 1361,  360,  646,  866, 1031, 1151, 1235, 1362,  361,
          //                647,  867, 1032, 1152, 1236,   75,  141,  196,  241,  277,  305,
          //                427,  482,  527,  563,  591,  702,  747,  783,  811,  912,  948,
          //                976, 1068, 1096, 1180, 1363,  362,  648,  868, 1033, 1153, 1237,
          //                 76,  142,  197,  242,  278,  306,  428,  483,  528,  564,  592,
          //                703,  748,  784,  812,  913,  949,  977, 1069, 1097, 1181,   10,
          //                 21,   31,   40,   48,   55,   87,   97,  106,  114,  121,  152,
          //                161,  169,  176,  206,  214,  221,  250,  257,  285,  373,  383,
          //                392,  400,  407,  438,  447,  455,  462,  492,  500,  507,  536,
          //                543,  571,  658,  667,  675,  682,  712,  720,  727,  756,  763,
          //                791,  877,  885,  892,  921,  928,  956, 1041, 1048, 1076, 1160,
          //               1364,  363,  649,  869, 1034, 1154, 1238,   77,  143,  198,  243,
          //                279,  307,  429,  484,  529,  565,  593,  704,  749,  785,  813,
          //                914,  950,  978, 1070, 1098, 1182,   11,   22,   32,   41,   49,
          //                 56,   88,   98,  107,  115,  122,  153,  162,  170,  177,  207,
          //                215,  222,  251,  258,  286,  374,  384,  393,  401,  408,  439,
          //                448,  456,  463,  493,  501,  508,  537,  544,  572,  659,  668,
          //                676,  683,  713,  721,  728,  757,  764,  792,  878,  886,  893,
          //                922,  929,  957, 1042, 1049, 1077, 1161,    0,    1,    2,    3,
          //                  4,    5,   12,   13,   14,   15,   16,   23,   24,   25,   26,
          //                 33,   34,   35,   42,   43,   50,   78,   79,   80,   81,   82,
          //                 89,   90,   91,   92,   99,  100,  101,  108,  109,  116,  144,
          //                145,  146,  147,  154,  155,  156,  163,  164,  171,  199,  200,
          //                201,  208,  209,  216,  244,  245,  252,  280,  364,  365,  366,
          //                367,  368,  375,  376,  377,  378,  385,  386,  387,  394,  395,
          //                402,  430,  431,  432,  433,  440,  441,  442,  449,  450,  457,
          //                485,  486,  487,  494,  495,  502,  530,  531,  538,  566,  650,
          //                651,  652,  653,  660,  661,  662,  669,  670,  677,  705,  706,
          //                707,  714,  715,  722,  750,  751,  758,  786,  870,  871,  872,
          //                879,  880,  887,  915,  916,  923,  951, 1035, 1036, 1043, 1071,
          //               1155},
          //              {1350, 1351, 1352, 1285, 1320, 1340,  354,  640,  860, 1025, 1145,
          //               1229, 1353, 1354, 1286, 1321, 1341,  355,  641,  861, 1026, 1146,
          //               1230, 1355, 1287, 1322, 1342,  356,  642,  862, 1027, 1147, 1231,
          //               1254, 1269, 1279,  323,  609,  829,  994, 1114, 1198, 1304, 1314,
          //                338,  624,  844, 1009, 1129, 1213, 1334,  348,  634,  854, 1019,
          //               1139, 1223,   72,  138,  193,  238,  274,  302,  424,  479,  524,
          //                560,  588,  699,  744,  780,  808,  909,  945,  973, 1065, 1093,
          //               1177, 1356, 1357, 1288, 1323, 1343,  357,  643,  863, 1028, 1148,
          //               1232, 1358, 1289, 1324, 1344,  358,  644,  864, 1029, 1149, 1233,
          //               1255, 1270, 1280,  324,  610,  830,  995, 1115, 1199, 1305, 1315,
          //                339,  625,  845, 1010, 1130, 1214, 1335,  349,  635,  855, 1020,
          //               1140, 1224,   73,  139,  194,  239,  275,  303,  425,  480,  525,
          //                561,  589,  700,  745,  781,  809,  910,  946,  974, 1066, 1094,
          //               1178, 1359, 1290, 1325, 1345,  359,  645,  865, 1030, 1150, 1234,
          //               1256, 1271, 1281,  325,  611,  831,  996, 1116, 1200, 1306, 1316,
          //                340,  626,  846, 1011, 1131, 1215, 1336,  350,  636,  856, 1021,
          //               1141, 1225,   74,  140,  195,  240,  276,  304,  426,  481,  526,
          //                562,  590,  701,  746,  782,  810,  911,  947,  975, 1067, 1095,
          //               1179, 1242, 1247, 1251,  311,  597,  817,  982, 1102, 1186, 1262,
          //               1266,  316,  602,  822,  987, 1107, 1191, 1276,  320,  606,  826,
          //                991, 1111, 1195,   60,  126,  181,  226,  262,  290,  412,  467,
          //                512,  548,  576,  687,  732,  768,  796,  897,  933,  961, 1053,
          //               1081, 1165, 1297, 1301,  331,  617,  837, 1002, 1122, 1206, 1311,
          //                335,  621,  841, 1006, 1126, 1210,   65,  131,  186,  231,  267,
          //                295,  417,  472,  517,  553,  581,  692,  737,  773,  801,  902,
          //                938,  966, 1058, 1086, 1170, 1331,  345,  631,  851, 1016, 1136,
          //               1220,   69,  135,  190,  235,  271,  299,  421,  476,  521,  557,
          //                585,  696,  741,  777,  805,  906,  942,  970, 1062, 1090, 1174,
          //                  9,   20,   30,   39,   47,   54,   86,   96,  105,  113,  120,
          //                151,  160,  168,  175,  205,  213,  220,  249,  256,  284,  372,
          //                382,  391,  399,  406,  437,  446,  454,  461,  491,  499,  506,
          //                535,  542,  570,  657,  666,  674,  681,  711,  719,  726,  755,
          //                762,  790,  876,  884,  891,  920,  927,  955, 1040, 1047, 1075,
          //               1159, 1360, 1361, 1291, 1326, 1346,  360,  646,  866, 1031, 1151,
          //               1235, 1362, 1292, 1327, 1347,  361,  647,  867, 1032, 1152, 1236,
          //               1257, 1272, 1282,  326,  612,  832,  997, 1117, 1201, 1307, 1317,
          //                341,  627,  847, 1012, 1132, 1216, 1337,  351,  637,  857, 1022,
          //               1142, 1226,   75,  141,  196,  241,  277,  305,  427,  482,  527,
          //                563,  591,  702,  747,  783,  811,  912,  948,  976, 1068, 1096,
          //               1180, 1363, 1293, 1328, 1348,  362,  648,  868, 1033, 1153, 1237,
          //               1258, 1273, 1283,  327,  613,  833,  998, 1118, 1202, 1308, 1318,
          //                342,  628,  848, 1013, 1133, 1217, 1338,  352,  638,  858, 1023,
          //               1143, 1227,   76,  142,  197,  242,  278,  306,  428,  483,  528,
          //                564,  592,  703,  748,  784,  812,  913,  949,  977, 1069, 1097,
          //               1181, 1243, 1248, 1252,  312,  598,  818,  983, 1103, 1187, 1263,
          //               1267,  317,  603,  823,  988, 1108, 1192, 1277,  321,  607,  827,
          //                992, 1112, 1196,   61,  127,  182,  227,  263,  291,  413,  468,
          //                513,  549,  577,  688,  733,  769,  797,  898,  934,  962, 1054,
          //               1082, 1166, 1298, 1302,  332,  618,  838, 1003, 1123, 1207, 1312,
          //                336,  622,  842, 1007, 1127, 1211,   66,  132,  187,  232,  268,
          //                296,  418,  473,  518,  554,  582,  693,  738,  774,  802,  903,
          //                939,  967, 1059, 1087, 1171, 1332,  346,  632,  852, 1017, 1137,
          //               1221,   70,  136,  191,  236,  272,  300,  422,  477,  522,  558,
          //                586,  697,  742,  778,  806,  907,  943,  971, 1063, 1091, 1175,
          //                 10,   21,   31,   40,   48,   55,   87,   97,  106,  114,  121,
          //                152,  161,  169,  176,  206,  214,  221,  250,  257,  285,  373,
          //                383,  392,  400,  407,  438,  447,  455,  462,  492,  500,  507,
          //                536,  543,  571,  658,  667,  675,  682,  712,  720,  727,  756,
          //                763,  791,  877,  885,  892,  921,  928,  956, 1041, 1048, 1076,
          //               1160, 1364, 1294, 1329, 1349,  363,  649,  869, 1034, 1154, 1238,
          //               1259, 1274, 1284,  328,  614,  834,  999, 1119, 1203, 1309, 1319,
          //                343,  629,  849, 1014, 1134, 1218, 1339,  353,  639,  859, 1024,
          //               1144, 1228,   77,  143,  198,  243,  279,  307,  429,  484,  529,
          //                565,  593,  704,  749,  785,  813,  914,  950,  978, 1070, 1098,
          //               1182, 1244, 1249, 1253,  313,  599,  819,  984, 1104, 1188, 1264,
          //               1268,  318,  604,  824,  989, 1109, 1193, 1278,  322,  608,  828,
          //                993, 1113, 1197,   62,  128,  183,  228,  264,  292,  414,  469,
          //                514,  550,  578,  689,  734,  770,  798,  899,  935,  963, 1055,
          //               1083, 1167, 1299, 1303,  333,  619,  839, 1004, 1124, 1208, 1313,
          //                337,  623,  843, 1008, 1128, 1212,   67,  133,  188,  233,  269,
          //                297,  419,  474,  519,  555,  583,  694,  739,  775,  803,  904,
          //                940,  968, 1060, 1088, 1172, 1333,  347,  633,  853, 1018, 1138,
          //               1222,   71,  137,  192,  237,  273,  301,  423,  478,  523,  559,
          //                587,  698,  743,  779,  807,  908,  944,  972, 1064, 1092, 1176,
          //                 11,   22,   32,   41,   49,   56,   88,   98,  107,  115,  122,
          //                153,  162,  170,  177,  207,  215,  222,  251,  258,  286,  374,
          //                384,  393,  401,  408,  439,  448,  456,  463,  493,  501,  508,
          //                537,  544,  572,  659,  668,  676,  683,  713,  721,  728,  757,
          //                764,  792,  878,  886,  893,  922,  929,  957, 1042, 1049, 1077,
          //               1161, 1239, 1240, 1241,  308,  594,  814,  979, 1099, 1183, 1245,
          //               1246,  309,  595,  815,  980, 1100, 1184, 1250,  310,  596,  816,
          //                981, 1101, 1185,   57,  123,  178,  223,  259,  287,  409,  464,
          //                509,  545,  573,  684,  729,  765,  793,  894,  930,  958, 1050,
          //               1078, 1162, 1260, 1261,  314,  600,  820,  985, 1105, 1189, 1265,
          //                315,  601,  821,  986, 1106, 1190,   58,  124,  179,  224,  260,
          //                288,  410,  465,  510,  546,  574,  685,  730,  766,  794,  895,
          //                931,  959, 1051, 1079, 1163, 1275,  319,  605,  825,  990, 1110,
          //               1194,   59,  125,  180,  225,  261,  289,  411,  466,  511,  547,
          //                575,  686,  731,  767,  795,  896,  932,  960, 1052, 1080, 1164,
          //                  6,   17,   27,   36,   44,   51,   83,   93,  102,  110,  117,
          //                148,  157,  165,  172,  202,  210,  217,  246,  253,  281,  369,
          //                379,  388,  396,  403,  434,  443,  451,  458,  488,  496,  503,
          //                532,  539,  567,  654,  663,  671,  678,  708,  716,  723,  752,
          //                759,  787,  873,  881,  888,  917,  924,  952, 1037, 1044, 1072,
          //               1156, 1295, 1296,  329,  615,  835, 1000, 1120, 1204, 1300,  330,
          //                616,  836, 1001, 1121, 1205,   63,  129,  184,  229,  265,  293,
          //                415,  470,  515,  551,  579,  690,  735,  771,  799,  900,  936,
          //                964, 1056, 1084, 1168, 1310,  334,  620,  840, 1005, 1125, 1209,
          //                 64,  130,  185,  230,  266,  294,  416,  471,  516,  552,  580,
          //                691,  736,  772,  800,  901,  937,  965, 1057, 1085, 1169,    7,
          //                 18,   28,   37,   45,   52,   84,   94,  103,  111,  118,  149,
          //                158,  166,  173,  203,  211,  218,  247,  254,  282,  370,  380,
          //                389,  397,  404,  435,  444,  452,  459,  489,  497,  504,  533,
          //                540,  568,  655,  664,  672,  679,  709,  717,  724,  753,  760,
          //                788,  874,  882,  889,  918,  925,  953, 1038, 1045, 1073, 1157,
          //               1330,  344,  630,  850, 1015, 1135, 1219,   68,  134,  189,  234,
          //                270,  298,  420,  475,  520,  556,  584,  695,  740,  776,  804,
          //                905,  941,  969, 1061, 1089, 1173,    8,   19,   29,   38,   46,
          //                 53,   85,   95,  104,  112,  119,  150,  159,  167,  174,  204,
          //                212,  219,  248,  255,  283,  371,  381,  390,  398,  405,  436,
          //                445,  453,  460,  490,  498,  505,  534,  541,  569,  656,  665,
          //                673,  680,  710,  718,  725,  754,  761,  789,  875,  883,  890,
          //                919,  926,  954, 1039, 1046, 1074, 1158,    0,    1,    2,    3,
          //                  4,    5,   12,   13,   14,   15,   16,   23,   24,   25,   26,
          //                 33,   34,   35,   42,   43,   50,   78,   79,   80,   81,   82,
          //                 89,   90,   91,   92,   99,  100,  101,  108,  109,  116,  144,
          //                145,  146,  147,  154,  155,  156,  163,  164,  171,  199,  200,
          //                201,  208,  209,  216,  244,  245,  252,  280,  364,  365,  366,
          //                367,  368,  375,  376,  377,  378,  385,  386,  387,  394,  395,
          //                402,  430,  431,  432,  433,  440,  441,  442,  449,  450,  457,
          //                485,  486,  487,  494,  495,  502,  530,  531,  538,  566,  650,
          //                651,  652,  653,  660,  661,  662,  669,  670,  677,  705,  706,
          //                707,  714,  715,  722,  750,  751,  758,  786,  870,  871,  872,
          //                879,  880,  887,  915,  916,  923,  951, 1035, 1036, 1043, 1071,
          //               1155}},
          //     
          //             {{1239, 1240, 1241, 1242, 1243, 1244,  979, 1099, 1183,  308,  594,
          //                814, 1245, 1246, 1247, 1248, 1249,  980, 1100, 1184,  309,  595,
          //                815, 1250, 1251, 1252, 1253,  981, 1101, 1185,  310,  596,  816,
          //               1254, 1255, 1256,  982, 1102, 1186,  311,  597,  817, 1257, 1258,
          //                983, 1103, 1187,  312,  598,  818, 1259,  984, 1104, 1188,  313,
          //                599,  819,  894,  930,  958,  223,  509,  729, 1050, 1078,  259,
          //                545,  765, 1162,  287,  573,  793,   57,  123,  178,  409,  464,
          //                684, 1260, 1261, 1262, 1263, 1264,  985, 1105, 1189,  314,  600,
          //                820, 1265, 1266, 1267, 1268,  986, 1106, 1190,  315,  601,  821,
          //               1269, 1270, 1271,  987, 1107, 1191,  316,  602,  822, 1272, 1273,
          //                988, 1108, 1192,  317,  603,  823, 1274,  989, 1109, 1193,  318,
          //                604,  824,  895,  931,  959,  224,  510,  730, 1051, 1079,  260,
          //                546,  766, 1163,  288,  574,  794,   58,  124,  179,  410,  465,
          //                685, 1275, 1276, 1277, 1278,  990, 1110, 1194,  319,  605,  825,
          //               1279, 1280, 1281,  991, 1111, 1195,  320,  606,  826, 1282, 1283,
          //                992, 1112, 1196,  321,  607,  827, 1284,  993, 1113, 1197,  322,
          //                608,  828,  896,  932,  960,  225,  511,  731, 1052, 1080,  261,
          //                547,  767, 1164,  289,  575,  795,   59,  125,  180,  411,  466,
          //                686, 1285, 1286, 1287,  994, 1114, 1198,  323,  609,  829, 1288,
          //               1289,  995, 1115, 1199,  324,  610,  830, 1290,  996, 1116, 1200,
          //                325,  611,  831,  897,  933,  961,  226,  512,  732, 1053, 1081,
          //                262,  548,  768, 1165,  290,  576,  796,   60,  126,  181,  412,
          //                467,  687, 1291, 1292,  997, 1117, 1201,  326,  612,  832, 1293,
          //                998, 1118, 1202,  327,  613,  833,  898,  934,  962,  227,  513,
          //                733, 1054, 1082,  263,  549,  769, 1166,  291,  577,  797,   61,
          //                127,  182,  413,  468,  688, 1294,  999, 1119, 1203,  328,  614,
          //                834,  899,  935,  963,  228,  514,  734, 1055, 1083,  264,  550,
          //                770, 1167,  292,  578,  798,   62,  128,  183,  414,  469,  689,
          //                873,  881,  888,  202,  488,  708,  917,  924,  210,  496,  716,
          //                952,  217,  503,  723,   36,  102,  157,  388,  443,  663, 1037,
          //               1044,  246,  532,  752, 1072,  253,  539,  759,   44,  110,  165,
          //                396,  451,  671, 1156,  281,  567,  787,   51,  117,  172,  403,
          //                458,  678,    6,   17,   27,   83,   93,  148,  369,  379,  434,
          //                654, 1295, 1296, 1297, 1298, 1299, 1000, 1120, 1204,  329,  615,
          //                835, 1300, 1301, 1302, 1303, 1001, 1121, 1205,  330,  616,  836,
          //               1304, 1305, 1306, 1002, 1122, 1206,  331,  617,  837, 1307, 1308,
          //               1003, 1123, 1207,  332,  618,  838, 1309, 1004, 1124, 1208,  333,
          //                619,  839,  900,  936,  964,  229,  515,  735, 1056, 1084,  265,
          //                551,  771, 1168,  293,  579,  799,   63,  129,  184,  415,  470,
          //                690, 1310, 1311, 1312, 1313, 1005, 1125, 1209,  334,  620,  840,
          //               1314, 1315, 1316, 1006, 1126, 1210,  335,  621,  841, 1317, 1318,
          //               1007, 1127, 1211,  336,  622,  842, 1319, 1008, 1128, 1212,  337,
          //                623,  843,  901,  937,  965,  230,  516,  736, 1057, 1085,  266,
          //                552,  772, 1169,  294,  580,  800,   64,  130,  185,  416,  471,
          //                691, 1320, 1321, 1322, 1009, 1129, 1213,  338,  624,  844, 1323,
          //               1324, 1010, 1130, 1214,  339,  625,  845, 1325, 1011, 1131, 1215,
          //                340,  626,  846,  902,  938,  966,  231,  517,  737, 1058, 1086,
          //                267,  553,  773, 1170,  295,  581,  801,   65,  131,  186,  417,
          //                472,  692, 1326, 1327, 1012, 1132, 1216,  341,  627,  847, 1328,
          //               1013, 1133, 1217,  342,  628,  848,  903,  939,  967,  232,  518,
          //                738, 1059, 1087,  268,  554,  774, 1171,  296,  582,  802,   66,
          //                132,  187,  418,  473,  693, 1329, 1014, 1134, 1218,  343,  629,
          //                849,  904,  940,  968,  233,  519,  739, 1060, 1088,  269,  555,
          //                775, 1172,  297,  583,  803,   67,  133,  188,  419,  474,  694,
          //                874,  882,  889,  203,  489,  709,  918,  925,  211,  497,  717,
          //                953,  218,  504,  724,   37,  103,  158,  389,  444,  664, 1038,
          //               1045,  247,  533,  753, 1073,  254,  540,  760,   45,  111,  166,
          //                397,  452,  672, 1157,  282,  568,  788,   52,  118,  173,  404,
          //                459,  679,    7,   18,   28,   84,   94,  149,  370,  380,  435,
          //                655, 1330, 1331, 1332, 1333, 1015, 1135, 1219,  344,  630,  850,
          //               1334, 1335, 1336, 1016, 1136, 1220,  345,  631,  851, 1337, 1338,
          //               1017, 1137, 1221,  346,  632,  852, 1339, 1018, 1138, 1222,  347,
          //                633,  853,  905,  941,  969,  234,  520,  740, 1061, 1089,  270,
          //                556,  776, 1173,  298,  584,  804,   68,  134,  189,  420,  475,
          //                695, 1340, 1341, 1342, 1019, 1139, 1223,  348,  634,  854, 1343,
          //               1344, 1020, 1140, 1224,  349,  635,  855, 1345, 1021, 1141, 1225,
          //                350,  636,  856,  906,  942,  970,  235,  521,  741, 1062, 1090,
          //                271,  557,  777, 1174,  299,  585,  805,   69,  135,  190,  421,
          //                476,  696, 1346, 1347, 1022, 1142, 1226,  351,  637,  857, 1348,
          //               1023, 1143, 1227,  352,  638,  858,  907,  943,  971,  236,  522,
          //                742, 1063, 1091,  272,  558,  778, 1175,  300,  586,  806,   70,
          //                136,  191,  422,  477,  697, 1349, 1024, 1144, 1228,  353,  639,
          //                859,  908,  944,  972,  237,  523,  743, 1064, 1092,  273,  559,
          //                779, 1176,  301,  587,  807,   71,  137,  192,  423,  478,  698,
          //                875,  883,  890,  204,  490,  710,  919,  926,  212,  498,  718,
          //                954,  219,  505,  725,   38,  104,  159,  390,  445,  665, 1039,
          //               1046,  248,  534,  754, 1074,  255,  541,  761,   46,  112,  167,
          //                398,  453,  673, 1158,  283,  569,  789,   53,  119,  174,  405,
          //                460,  680,    8,   19,   29,   85,   95,  150,  371,  381,  436,
          //                656, 1350, 1351, 1352, 1025, 1145, 1229,  354,  640,  860, 1353,
          //               1354, 1026, 1146, 1230,  355,  641,  861, 1355, 1027, 1147, 1231,
          //                356,  642,  862,  909,  945,  973,  238,  524,  744, 1065, 1093,
          //                274,  560,  780, 1177,  302,  588,  808,   72,  138,  193,  424,
          //                479,  699, 1356, 1357, 1028, 1148, 1232,  357,  643,  863, 1358,
          //               1029, 1149, 1233,  358,  644,  864,  910,  946,  974,  239,  525,
          //                745, 1066, 1094,  275,  561,  781, 1178,  303,  589,  809,   73,
          //                139,  194,  425,  480,  700, 1359, 1030, 1150, 1234,  359,  645,
          //                865,  911,  947,  975,  240,  526,  746, 1067, 1095,  276,  562,
          //                782, 1179,  304,  590,  810,   74,  140,  195,  426,  481,  701,
          //                876,  884,  891,  205,  491,  711,  920,  927,  213,  499,  719,
          //                955,  220,  506,  726,   39,  105,  160,  391,  446,  666, 1040,
          //               1047,  249,  535,  755, 1075,  256,  542,  762,   47,  113,  168,
          //                399,  454,  674, 1159,  284,  570,  790,   54,  120,  175,  406,
          //                461,  681,    9,   20,   30,   86,   96,  151,  372,  382,  437,
          //                657, 1360, 1361, 1031, 1151, 1235,  360,  646,  866, 1362, 1032,
          //               1152, 1236,  361,  647,  867,  912,  948,  976,  241,  527,  747,
          //               1068, 1096,  277,  563,  783, 1180,  305,  591,  811,   75,  141,
          //                196,  427,  482,  702, 1363, 1033, 1153, 1237,  362,  648,  868,
          //                913,  949,  977,  242,  528,  748, 1069, 1097,  278,  564,  784,
          //               1181,  306,  592,  812,   76,  142,  197,  428,  483,  703,  877,
          //                885,  892,  206,  492,  712,  921,  928,  214,  500,  720,  956,
          //                221,  507,  727,   40,  106,  161,  392,  447,  667, 1041, 1048,
          //                250,  536,  756, 1076,  257,  543,  763,   48,  114,  169,  400,
          //                455,  675, 1160,  285,  571,  791,   55,  121,  176,  407,  462,
          //                682,   10,   21,   31,   87,   97,  152,  373,  383,  438,  658,
          //               1364, 1034, 1154, 1238,  363,  649,  869,  914,  950,  978,  243,
          //                529,  749, 1070, 1098,  279,  565,  785, 1182,  307,  593,  813,
          //                 77,  143,  198,  429,  484,  704,  878,  886,  893,  207,  493,
          //                713,  922,  929,  215,  501,  721,  957,  222,  508,  728,   41,
          //                107,  162,  393,  448,  668, 1042, 1049,  251,  537,  757, 1077,
          //                258,  544,  764,   49,  115,  170,  401,  456,  676, 1161,  286,
          //                572,  792,   56,  122,  177,  408,  463,  683,   11,   22,   32,
          //                 88,   98,  153,  374,  384,  439,  659,  870,  871,  872,  199,
          //                485,  705,  879,  880,  200,  486,  706,  887,  201,  487,  707,
          //                 33,   99,  154,  385,  440,  660,  915,  916,  208,  494,  714,
          //                923,  209,  495,  715,   34,  100,  155,  386,  441,  661,  951,
          //                216,  502,  722,   35,  101,  156,  387,  442,  662,    3,   14,
          //                 24,   80,   90,  145,  366,  376,  431,  651, 1035, 1036,  244,
          //                530,  750, 1043,  245,  531,  751,   42,  108,  163,  394,  449,
          //                669, 1071,  252,  538,  758,   43,  109,  164,  395,  450,  670,
          //                  4,   15,   25,   81,   91,  146,  367,  377,  432,  652, 1155,
          //                280,  566,  786,   50,  116,  171,  402,  457,  677,    5,   16,
          //                 26,   82,   92,  147,  368,  378,  433,  653,    0,    1,    2,
          //                 12,   13,   23,   78,   79,   89,  144,  364,  365,  375,  430,
          //                650},
          //              {1350, 1351, 1352, 1285, 1320, 1340, 1025, 1145, 1229,  354,  640,
          //                860, 1353, 1354, 1286, 1321, 1341, 1026, 1146, 1230,  355,  641,
          //                861, 1355, 1287, 1322, 1342, 1027, 1147, 1231,  356,  642,  862,
          //               1254, 1269, 1279,  994, 1114, 1198,  323,  609,  829, 1304, 1314,
          //               1009, 1129, 1213,  338,  624,  844, 1334, 1019, 1139, 1223,  348,
          //                634,  854,  909,  945,  973,  238,  524,  744, 1065, 1093,  274,
          //                560,  780, 1177,  302,  588,  808,   72,  138,  193,  424,  479,
          //                699, 1356, 1357, 1288, 1323, 1343, 1028, 1148, 1232,  357,  643,
          //                863, 1358, 1289, 1324, 1344, 1029, 1149, 1233,  358,  644,  864,
          //               1255, 1270, 1280,  995, 1115, 1199,  324,  610,  830, 1305, 1315,
          //               1010, 1130, 1214,  339,  625,  845, 1335, 1020, 1140, 1224,  349,
          //                635,  855,  910,  946,  974,  239,  525,  745, 1066, 1094,  275,
          //                561,  781, 1178,  303,  589,  809,   73,  139,  194,  425,  480,
          //                700, 1359, 1290, 1325, 1345, 1030, 1150, 1234,  359,  645,  865,
          //               1256, 1271, 1281,  996, 1116, 1200,  325,  611,  831, 1306, 1316,
          //               1011, 1131, 1215,  340,  626,  846, 1336, 1021, 1141, 1225,  350,
          //                636,  856,  911,  947,  975,  240,  526,  746, 1067, 1095,  276,
          //                562,  782, 1179,  304,  590,  810,   74,  140,  195,  426,  481,
          //                701, 1242, 1247, 1251,  982, 1102, 1186,  311,  597,  817, 1262,
          //               1266,  987, 1107, 1191,  316,  602,  822, 1276,  991, 1111, 1195,
          //                320,  606,  826,  897,  933,  961,  226,  512,  732, 1053, 1081,
          //                262,  548,  768, 1165,  290,  576,  796,   60,  126,  181,  412,
          //                467,  687, 1297, 1301, 1002, 1122, 1206,  331,  617,  837, 1311,
          //               1006, 1126, 1210,  335,  621,  841,  902,  938,  966,  231,  517,
          //                737, 1058, 1086,  267,  553,  773, 1170,  295,  581,  801,   65,
          //                131,  186,  417,  472,  692, 1331, 1016, 1136, 1220,  345,  631,
          //                851,  906,  942,  970,  235,  521,  741, 1062, 1090,  271,  557,
          //                777, 1174,  299,  585,  805,   69,  135,  190,  421,  476,  696,
          //                876,  884,  891,  205,  491,  711,  920,  927,  213,  499,  719,
          //                955,  220,  506,  726,   39,  105,  160,  391,  446,  666, 1040,
          //               1047,  249,  535,  755, 1075,  256,  542,  762,   47,  113,  168,
          //                399,  454,  674, 1159,  284,  570,  790,   54,  120,  175,  406,
          //                461,  681,    9,   20,   30,   86,   96,  151,  372,  382,  437,
          //                657, 1360, 1361, 1291, 1326, 1346, 1031, 1151, 1235,  360,  646,
          //                866, 1362, 1292, 1327, 1347, 1032, 1152, 1236,  361,  647,  867,
          //               1257, 1272, 1282,  997, 1117, 1201,  326,  612,  832, 1307, 1317,
          //               1012, 1132, 1216,  341,  627,  847, 1337, 1022, 1142, 1226,  351,
          //                637,  857,  912,  948,  976,  241,  527,  747, 1068, 1096,  277,
          //                563,  783, 1180,  305,  591,  811,   75,  141,  196,  427,  482,
          //                702, 1363, 1293, 1328, 1348, 1033, 1153, 1237,  362,  648,  868,
          //               1258, 1273, 1283,  998, 1118, 1202,  327,  613,  833, 1308, 1318,
          //               1013, 1133, 1217,  342,  628,  848, 1338, 1023, 1143, 1227,  352,
          //                638,  858,  913,  949,  977,  242,  528,  748, 1069, 1097,  278,
          //                564,  784, 1181,  306,  592,  812,   76,  142,  197,  428,  483,
          //                703, 1243, 1248, 1252,  983, 1103, 1187,  312,  598,  818, 1263,
          //               1267,  988, 1108, 1192,  317,  603,  823, 1277,  992, 1112, 1196,
          //                321,  607,  827,  898,  934,  962,  227,  513,  733, 1054, 1082,
          //                263,  549,  769, 1166,  291,  577,  797,   61,  127,  182,  413,
          //                468,  688, 1298, 1302, 1003, 1123, 1207,  332,  618,  838, 1312,
          //               1007, 1127, 1211,  336,  622,  842,  903,  939,  967,  232,  518,
          //                738, 1059, 1087,  268,  554,  774, 1171,  296,  582,  802,   66,
          //                132,  187,  418,  473,  693, 1332, 1017, 1137, 1221,  346,  632,
          //                852,  907,  943,  971,  236,  522,  742, 1063, 1091,  272,  558,
          //                778, 1175,  300,  586,  806,   70,  136,  191,  422,  477,  697,
          //                877,  885,  892,  206,  492,  712,  921,  928,  214,  500,  720,
          //                956,  221,  507,  727,   40,  106,  161,  392,  447,  667, 1041,
          //               1048,  250,  536,  756, 1076,  257,  543,  763,   48,  114,  169,
          //                400,  455,  675, 1160,  285,  571,  791,   55,  121,  176,  407,
          //                462,  682,   10,   21,   31,   87,   97,  152,  373,  383,  438,
          //                658, 1364, 1294, 1329, 1349, 1034, 1154, 1238,  363,  649,  869,
          //               1259, 1274, 1284,  999, 1119, 1203,  328,  614,  834, 1309, 1319,
          //               1014, 1134, 1218,  343,  629,  849, 1339, 1024, 1144, 1228,  353,
          //                639,  859,  914,  950,  978,  243,  529,  749, 1070, 1098,  279,
          //                565,  785, 1182,  307,  593,  813,   77,  143,  198,  429,  484,
          //                704, 1244, 1249, 1253,  984, 1104, 1188,  313,  599,  819, 1264,
          //               1268,  989, 1109, 1193,  318,  604,  824, 1278,  993, 1113, 1197,
          //                322,  608,  828,  899,  935,  963,  228,  514,  734, 1055, 1083,
          //                264,  550,  770, 1167,  292,  578,  798,   62,  128,  183,  414,
          //                469,  689, 1299, 1303, 1004, 1124, 1208,  333,  619,  839, 1313,
          //               1008, 1128, 1212,  337,  623,  843,  904,  940,  968,  233,  519,
          //                739, 1060, 1088,  269,  555,  775, 1172,  297,  583,  803,   67,
          //                133,  188,  419,  474,  694, 1333, 1018, 1138, 1222,  347,  633,
          //                853,  908,  944,  972,  237,  523,  743, 1064, 1092,  273,  559,
          //                779, 1176,  301,  587,  807,   71,  137,  192,  423,  478,  698,
          //                878,  886,  893,  207,  493,  713,  922,  929,  215,  501,  721,
          //                957,  222,  508,  728,   41,  107,  162,  393,  448,  668, 1042,
          //               1049,  251,  537,  757, 1077,  258,  544,  764,   49,  115,  170,
          //                401,  456,  676, 1161,  286,  572,  792,   56,  122,  177,  408,
          //                463,  683,   11,   22,   32,   88,   98,  153,  374,  384,  439,
          //                659, 1239, 1240, 1241,  979, 1099, 1183,  308,  594,  814, 1245,
          //               1246,  980, 1100, 1184,  309,  595,  815, 1250,  981, 1101, 1185,
          //                310,  596,  816,  894,  930,  958,  223,  509,  729, 1050, 1078,
          //                259,  545,  765, 1162,  287,  573,  793,   57,  123,  178,  409,
          //                464,  684, 1260, 1261,  985, 1105, 1189,  314,  600,  820, 1265,
          //                986, 1106, 1190,  315,  601,  821,  895,  931,  959,  224,  510,
          //                730, 1051, 1079,  260,  546,  766, 1163,  288,  574,  794,   58,
          //                124,  179,  410,  465,  685, 1275,  990, 1110, 1194,  319,  605,
          //                825,  896,  932,  960,  225,  511,  731, 1052, 1080,  261,  547,
          //                767, 1164,  289,  575,  795,   59,  125,  180,  411,  466,  686,
          //                873,  881,  888,  202,  488,  708,  917,  924,  210,  496,  716,
          //                952,  217,  503,  723,   36,  102,  157,  388,  443,  663, 1037,
          //               1044,  246,  532,  752, 1072,  253,  539,  759,   44,  110,  165,
          //                396,  451,  671, 1156,  281,  567,  787,   51,  117,  172,  403,
          //                458,  678,    6,   17,   27,   83,   93,  148,  369,  379,  434,
          //                654, 1295, 1296, 1000, 1120, 1204,  329,  615,  835, 1300, 1001,
          //               1121, 1205,  330,  616,  836,  900,  936,  964,  229,  515,  735,
          //               1056, 1084,  265,  551,  771, 1168,  293,  579,  799,   63,  129,
          //                184,  415,  470,  690, 1310, 1005, 1125, 1209,  334,  620,  840,
          //                901,  937,  965,  230,  516,  736, 1057, 1085,  266,  552,  772,
          //               1169,  294,  580,  800,   64,  130,  185,  416,  471,  691,  874,
          //                882,  889,  203,  489,  709,  918,  925,  211,  497,  717,  953,
          //                218,  504,  724,   37,  103,  158,  389,  444,  664, 1038, 1045,
          //                247,  533,  753, 1073,  254,  540,  760,   45,  111,  166,  397,
          //                452,  672, 1157,  282,  568,  788,   52,  118,  173,  404,  459,
          //                679,    7,   18,   28,   84,   94,  149,  370,  380,  435,  655,
          //               1330, 1015, 1135, 1219,  344,  630,  850,  905,  941,  969,  234,
          //                520,  740, 1061, 1089,  270,  556,  776, 1173,  298,  584,  804,
          //                 68,  134,  189,  420,  475,  695,  875,  883,  890,  204,  490,
          //                710,  919,  926,  212,  498,  718,  954,  219,  505,  725,   38,
          //                104,  159,  390,  445,  665, 1039, 1046,  248,  534,  754, 1074,
          //                255,  541,  761,   46,  112,  167,  398,  453,  673, 1158,  283,
          //                569,  789,   53,  119,  174,  405,  460,  680,    8,   19,   29,
          //                 85,   95,  150,  371,  381,  436,  656,  870,  871,  872,  199,
          //                485,  705,  879,  880,  200,  486,  706,  887,  201,  487,  707,
          //                 33,   99,  154,  385,  440,  660,  915,  916,  208,  494,  714,
          //                923,  209,  495,  715,   34,  100,  155,  386,  441,  661,  951,
          //                216,  502,  722,   35,  101,  156,  387,  442,  662,    3,   14,
          //                 24,   80,   90,  145,  366,  376,  431,  651, 1035, 1036,  244,
          //                530,  750, 1043,  245,  531,  751,   42,  108,  163,  394,  449,
          //                669, 1071,  252,  538,  758,   43,  109,  164,  395,  450,  670,
          //                  4,   15,   25,   81,   91,  146,  367,  377,  432,  652, 1155,
          //                280,  566,  786,   50,  116,  171,  402,  457,  677,    5,   16,
          //                 26,   82,   92,  147,  368,  378,  433,  653,    0,    1,    2,
          //                 12,   13,   23,   78,   79,   89,  144,  364,  365,  375,  430,
          //                650}}}};
          //        s_target = mapDerivIndex4_xxxx[swap_braket][swap_tbra][swap_tket][s];
          //      }
          //        break;
          //      // TODO 3center and 2center -- unable to test these, so commented out. Array should be right.
          //      // Requires definitions of LIBINT2_CENTER_DEPENDENT_MAX_AM_3eri3/_3eri4 to work above
          //      //case BraKet::xs_xx: {
          //      //  assert(swap_bra == false);
          //      //  assert(swap_braket == false);
          //      //  const unsigned mapDerivIndex4_xsxx[2][495] = {
          //      //      {  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
          //      //         13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
          //      //         26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
          //      //         39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
          //      //         52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
          //      //         65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
          //      //         78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
          //      //         91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
          //      //        104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
          //      //        117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
          //      //        130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
          //      //        143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
          //      //        156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
          //      //        169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
          //      //        182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
          //      //        195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
          //      //        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
          //      //        221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
          //      //        234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
          //      //        247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,
          //      //        260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
          //      //        273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285,
          //      //        286, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 297, 298,
          //      //        299, 300, 301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311,
          //      //        312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324,
          //      //        325, 326, 327, 328, 329, 330, 331, 332, 333, 334, 335, 336, 337,
          //      //        338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350,
          //      //        351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 363,
          //      //        364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376,
          //      //        377, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389,
          //      //        390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402,
          //      //        403, 404, 405, 406, 407, 408, 409, 410, 411, 412, 413, 414, 415,
          //      //        416, 417, 418, 419, 420, 421, 422, 423, 424, 425, 426, 427, 428,
          //      //        429, 430, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441,
          //      //        442, 443, 444, 445, 446, 447, 448, 449, 450, 451, 452, 453, 454,
          //      //        455, 456, 457, 458, 459, 460, 461, 462, 463, 464, 465, 466, 467,
          //      //        468, 469, 470, 471, 472, 473, 474, 475, 476, 477, 478, 479, 480,
          //      //        481, 482, 483, 484, 485, 486, 487, 488, 489, 490, 491, 492, 493,
          //      //        494},
          //      //       {  0,   1,   2,   6,   7,   8,   3,   4,   5,   9,  10,  14,  15,
          //      //         16,  11,  12,  13,  17,  21,  22,  23,  18,  19,  20,  39,  40,
          //      //         41,  27,  32,  36,  42,  43,  28,  33,  37,  44,  29,  34,  38,
          //      //         24,  25,  26,  30,  31,  35,  45,  46,  50,  51,  52,  47,  48,
          //      //         49,  53,  57,  58,  59,  54,  55,  56,  75,  76,  77,  63,  68,
          //      //         72,  78,  79,  64,  69,  73,  80,  65,  70,  74,  60,  61,  62,
          //      //         66,  67,  71,  81,  85,  86,  87,  82,  83,  84, 103, 104, 105,
          //      //         91,  96, 100, 106, 107,  92,  97, 101, 108,  93,  98, 102,  88,
          //      //         89,  90,  94,  95,  99, 155, 156, 157, 124, 139, 149, 158, 159,
          //      //        125, 140, 150, 160, 126, 141, 151, 112, 117, 121, 132, 136, 146,
          //      //        161, 162, 127, 142, 152, 163, 128, 143, 153, 113, 118, 122, 133,
          //      //        137, 147, 164, 129, 144, 154, 114, 119, 123, 134, 138, 148, 109,
          //      //        110, 111, 115, 116, 120, 130, 131, 135, 145, 165, 166, 170, 171,
          //      //        172, 167, 168, 169, 173, 177, 178, 179, 174, 175, 176, 195, 196,
          //      //        197, 183, 188, 192, 198, 199, 184, 189, 193, 200, 185, 190, 194,
          //      //        180, 181, 182, 186, 187, 191, 201, 205, 206, 207, 202, 203, 204,
          //      //        223, 224, 225, 211, 216, 220, 226, 227, 212, 217, 221, 228, 213,
          //      //        218, 222, 208, 209, 210, 214, 215, 219, 275, 276, 277, 244, 259,
          //      //        269, 278, 279, 245, 260, 270, 280, 246, 261, 271, 232, 237, 241,
          //      //        252, 256, 266, 281, 282, 247, 262, 272, 283, 248, 263, 273, 233,
          //      //        238, 242, 253, 257, 267, 284, 249, 264, 274, 234, 239, 243, 254,
          //      //        258, 268, 229, 230, 231, 235, 236, 240, 250, 251, 255, 265, 285,
          //      //        289, 290, 291, 286, 287, 288, 307, 308, 309, 295, 300, 304, 310,
          //      //        311, 296, 301, 305, 312, 297, 302, 306, 292, 293, 294, 298, 299,
          //      //        303, 359, 360, 361, 328, 343, 353, 362, 363, 329, 344, 354, 364,
          //      //        330, 345, 355, 316, 321, 325, 336, 340, 350, 365, 366, 331, 346,
          //      //        356, 367, 332, 347, 357, 317, 322, 326, 337, 341, 351, 368, 333,
          //      //        348, 358, 318, 323, 327, 338, 342, 352, 313, 314, 315, 319, 320,
          //      //        324, 334, 335, 339, 349, 480, 481, 482, 415, 450, 470, 483, 484,
          //      //        416, 451, 471, 485, 417, 452, 472, 384, 399, 409, 434, 444, 464,
          //      //        486, 487, 418, 453, 473, 488, 419, 454, 474, 385, 400, 410, 435,
          //      //        445, 465, 489, 420, 455, 475, 386, 401, 411, 436, 446, 466, 372,
          //      //        377, 381, 392, 396, 406, 427, 431, 441, 461, 490, 491, 421, 456,
          //      //        476, 492, 422, 457, 477, 387, 402, 412, 437, 447, 467, 493, 423,
          //      //        458, 478, 388, 403, 413, 438, 448, 468, 373, 378, 382, 393, 397,
          //      //        407, 428, 432, 442, 462, 494, 424, 459, 479, 389, 404, 414, 439,
          //      //        449, 469, 374, 379, 383, 394, 398, 408, 429, 433, 443, 463, 369,
          //      //        370, 371, 375, 376, 380, 390, 391, 395, 405, 425, 426, 430, 440,
          //      //        460}};
          //      //  s_target = mapDerivIndex4_xsxx[swap_tket][s];
          //      //}
          //      //  break;

          //      //case BraKet::xs_xs: {
          //      //  assert(swap_bra == false);
          //      //  assert(swap_ket == false);
          //      //  assert(swap_braket == false);
          //      //  s_target = s;
          //      //}
          //      //  break;

          //      default:
          //        assert(false && "this backet type not yet supported for 4th geometric derivatives");
          //    }
          //  } break;

          //  default:
          //    assert(false &&
          //           "5-th and higher derivatives not yet generalized");
          //}

          for (auto r1 = 0; r1 != nr1; ++r1) {
            for (auto r2 = 0; r2 != nr2; ++r2, src_row_ptr += ncol) {
              typedef Eigen::Map<const Matrix> ConstMap;
              typedef Eigen::Map<Matrix> Map;
              typedef Eigen::Map<Matrix, Eigen::Unaligned,
                                 Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>
                  StridedMap;

              // represent this source row as a matrix
              ConstMap src_blk_mat(src_row_ptr, nc1, nc2);

              // and copy to the block of the target matrix
              if (swap_braket) {
                // if swapped bra and ket, a row of source becomes a column
                // of
                // target
                // source row {r1,r2} is mapped to target column {r1,r2} if
                // !swap_tket, else to {r2,r1}
                const auto tgt_col_idx =
                    !swap_tket ? r1 * nr2 + r2 : r2 * nr1 + r1;
                StridedMap tgt_blk_mat(
                    tgt_ptr + tgt_col_idx, nr1_tgt, nr2_tgt,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>(
                        nr2_tgt * ncol_tgt, ncol_tgt));
                if (swap_tbra)
                  tgt_blk_mat = src_blk_mat.transpose();
                else
                  tgt_blk_mat = src_blk_mat;
              } else {
                // source row {r1,r2} is mapped to target row {r1,r2} if
                // !swap_tbra, else to {r2,r1}
                const auto tgt_row_idx =
                    !swap_tbra ? r1 * nr2 + r2 : r2 * nr1 + r1;
                Map tgt_blk_mat(tgt_ptr + tgt_row_idx * ncol, nc1_tgt, nc2_tgt);
                if (swap_tket)
                  tgt_blk_mat = src_blk_mat.transpose();
                else
                  tgt_blk_mat = src_blk_mat;
              }
            }  // end of loop
          }    // over rows of source
          std::swap(source, target);
        }  // need to permute?

        // if the integrals ended up in scratch_, keep them there, update the
        // hot buffer
        // to the next available scratch space, and update targets_
        if (source != primdata_[0].targets[s]) {
          hotscr += n1234_cart;
          if (s != s_target)
            assert(set_targets_ && "logic error");  // mess if targets_ points
                                                    // to primdata_[0].targets
          targets_[s_target] = source;
        } else {
          // only needed if permuted derivs or set_targets_ is true
          // for simplicity always set targets_
          if (s != s_target)
            assert(set_targets_ && "logic error");  // mess if targets_ points
                                                    // to primdata_[0].targets
          targets_[s_target] = source;
        }
      }     // loop over shellsets
    }       // if need_scratch => needed to transpose and/or tform
    else {  // did not use scratch? may still need to update targets_
      if (set_targets_) {
        for (auto s = 0; s != ntargets; ++s)
          targets_[s] = primdata_[0].targets[s];
      }
    }

#ifdef LIBINT2_ENGINE_TIMERS
    const auto t2 = timers.stop(2);
#ifdef LIBINT2_ENGINE_PROFILE_CLASS
    class_profiles[id].tform += t2.count();
#endif
#endif
  }  // not (ss|ss)

  if (cartesian_shell_normalization() == CartesianShellNormalization::uniform) {
    std::array<std::reference_wrapper<const Shell>, 4> shells{bra1, bra2, ket1, ket2};
    for (auto s = 0ul; s != targets_.size(); ++s) {
      uniform_normalize_cartesian_shells(const_cast<value_type*>(targets_[s]), shells);
    }
  }

  return targets_;
}

#undef BOOST_PP_NBODY_OPERATOR_LIST
#undef BOOST_PP_NBODY_OPERATOR_INDEX_TUPLE
#undef BOOST_PP_NBODY_OPERATOR_INDEX_LIST
#undef BOOST_PP_NBODY_BRAKET_INDEX_TUPLE
#undef BOOST_PP_NBODY_BRAKET_INDEX_LIST
#undef BOOST_PP_NBODY_DERIV_ORDER_TUPLE
#undef BOOST_PP_NBODY_DERIV_ORDER_LIST
#undef BOOST_PP_NBODYENGINE_MCR3
#undef BOOST_PP_NBODYENGINE_MCR3_ncenter
#undef BOOST_PP_NBODYENGINE_MCR3_default_ncenter
#undef BOOST_PP_NBODYENGINE_MCR3_NCENTER
#undef BOOST_PP_NBODYENGINE_MCR3_OPER
#undef BOOST_PP_NBODYENGINE_MCR3_DERIV
#undef BOOST_PP_NBODYENGINE_MCR3_task
#undef BOOST_PP_NBODYENGINE_MCR3_TASK
#undef BOOST_PP_NBODYENGINE_MCR4
#undef BOOST_PP_NBODYENGINE_MCR5
#undef BOOST_PP_NBODYENGINE_MCR6
#undef BOOST_PP_NBODYENGINE_MCR7

#ifdef LIBINT2_DOES_NOT_INLINE_ENGINE
template any Engine::enforce_params_type<Engine::empty_pod>(
    Operator oper, const Engine::empty_pod& params, bool throw_if_wrong_type);

template const Engine::target_ptr_vec& Engine::compute<Shell>(
    const Shell& first_shell, const Shell&);

template const Engine::target_ptr_vec& Engine::compute<Shell, Shell>(
    const Shell& first_shell, const Shell&, const Shell&);

template const Engine::target_ptr_vec& Engine::compute<Shell, Shell, Shell>(
    const Shell& first_shell, const Shell&, const Shell&, const Shell&);
#endif

}  // namespace libint2

#endif /* _libint2_src_lib_libint_engineimpl_h_ */
