/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_include_derivmap_h_
#define _libint2_include_derivmap_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
#error "The simple Libint API requires C++11 support"
#endif

#include <libint2/braket.h>
#include <libint2/tensor.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <vector>

namespace libint2 {
/** This class statically initializes all index permutation maps
 * for each BraKet type which requires them in Engine. A multidimensional array
 * is created for all deriv_orders up to LIBINT2_MAX_DERIV_ORDER for each BraKet
 * type. Each of these index permutation maps have 4 dimensions of size 2,2,2
 * and nderivs, where nderivs is the total number of geometric derivatives for
 * this shell set at the particular order of differentiation. The 3 dimensions
 * of size 2 each indicate whether the permutation has the effect of 1. swapping
 * the bra and ket
 * 2. swapping the bra and 3. swapping the ket.
 * These maps are created when libint2::initialize() is called.
 */
class DerivMapGenerator {
 public:
  // Statically generate all mapDerivIndex arrays up to LIBINT2_MAX_DERIV_ORDER
  // for both BraKet::xx_xx and BraKet::xs_xx and store in a vector
  static void initialize() {
    int max_deriv_order = 0;
#if defined(LIBINT2_MAX_DERIV_ORDER)
    max_deriv_order = LIBINT2_MAX_DERIV_ORDER;
#else  // defined(LIBINT2_MAX_DERIV_ORDER)
#ifdef INCLUDE_ERI
    if (INCLUDE_ERI > max_deriv_order) max_deriv_order = INCLUDE_ERI;
#endif  // INCLUDE_ERI
#ifdef INCLUDE_ERI3
    if (INCLUDE_ERI3 > max_deriv_order) max_deriv_order = INCLUDE_ERI3;
#endif  // INCLUDE_ERI3
#endif  // defined(LIBINT2_MAX_DERIV_ORDER)
    for (int d = 1; d <= max_deriv_order; ++d) {
      braket_xx_xx().push_back(
          DerivMapGenerator::generate_deriv_index_map(d, BraKet::xx_xx));
      braket_xs_xx().push_back(
          DerivMapGenerator::generate_deriv_index_map(d, BraKet::xs_xx));
    }
  }

  // Helper functions for declaring/building/calling the maps in both
  // initialize() and instance()
  static std::vector<Tensor<size_t>>& braket_xx_xx() {
    static std::vector<Tensor<size_t>> braket_xx_xx_maps;
    return braket_xx_xx_maps;
  }

  static std::vector<Tensor<size_t>>& braket_xs_xx() {
    static std::vector<Tensor<size_t>> braket_xs_xx_maps;
    return braket_xs_xx_maps;
  }

  static Tensor<size_t>& instance(int deriv_order_, BraKet braket_) {
    // access and return the one that is needed
    switch (braket_) {
      case BraKet::xx_xx: {
        return braket_xx_xx()[deriv_order_ - 1];
      }
      case BraKet::xs_xx: {
        return braket_xs_xx()[deriv_order_ - 1];
      }
      default:
        assert(false &&
               "Derivative index permutation maps for this braket type not yet "
               "supported.");
    }
    std::abort();
  }

 private:
  // Functions required for generating the maps.
  // Combinations with repetition.
  // Combinations of size 'k' from 'n' elements stored in vector 'inp'.
  // Requires instantiating vector 'out' and vector of vectors 'result', which
  // stores every combination.
  static void cwr_recursion(std::vector<size_t> inp, std::vector<size_t>& out,
                            std::vector<std::vector<size_t>>& result, size_t k,
                            size_t i, size_t n) {
    // base case: if combination size is k, add to result
    if (out.size() == k) {
      result.push_back(out);
      return;
    }
    for (size_t j = i; j < n; j++) {
      out.push_back(inp[j]);
      cwr_recursion(inp, out, result, k, j, n);
      // backtrack - remove current element from solution
      out.pop_back();
    }
  }

  // How many shell set derivatives
  // k is number of centers, n is deriv order.
  // Assumes coordinate-independent operator, 3 coordinates per center
  static size_t how_many_derivs(size_t k, size_t n) {
    size_t val = 1;
    size_t factorial = 1;
    for (size_t i = 0; i < n; i++) {
      val *= (3 * k + i);
      factorial *= i + 1;
    }
    val /= factorial;
    return val;
  }

  // Create array which is of size (nderivs, deriv_order)
  // which when indexed with a 1d buffer index (the flattened generalized upper
  // triangle index) it maps it to the corresponding to multidimensional index
  // tuple (i,j,k,...) where i<=j<=k<=... The multidimensional indices each
  // correspond to a partial derivative operator, and the value of each index
  // indicates which parameter is being differentiated. For example, if there
  // are 6 differentiable parameters, and deriv_order = 2, calling
  // generate_multi_index_lookup(6,2) yields 'combos' such that combos[13]
  // yields the pair [2,4] which matches the value and index pair in the below
  // array: 0  1  2  3  4  5
  //    6  7  8  9 10
  //      11 12 13 14
  //         15 16 17
  //            18 19
  //               20
  static std::vector<std::vector<size_t>> generate_multi_index_lookup(
      size_t nparams, size_t deriv_order) {
    using namespace std;
    // Generate vector of indices 0 through nparams-1
    vector<size_t> inp;
    for (size_t i = 0; i < nparams; i++) {
      inp.push_back(i);
    }
    // Generate all possible combinations with repitition.
    // These are upper triangle indices, and the length of them is the total
    // number of derivatives
    vector<size_t> out;
    vector<vector<size_t>> combos;
    cwr_recursion(inp, out, combos, deriv_order, 0, nparams);
    return combos;
  }

  // Function which computes mapDerivIndex array
  static Tensor<size_t> generate_deriv_index_map(size_t deriv_order,
                                                 BraKet braket_) {
    using namespace std;
    // Check BraKet type to determine permutations, number of centers
    vector<size_t> swap_braket_perm;
    vector<size_t> swap_bra_perm;
    vector<size_t> swap_ket_perm;
    vector<vector<size_t>> swap_combos;
    size_t ncenters = 0;
    switch (braket_) {
      case BraKet::xx_xx: {
        swap_braket_perm = {6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 5};
        swap_bra_perm = {3, 4, 5, 0, 1, 2, 6, 7, 8, 9, 10, 11};
        swap_ket_perm = {0, 1, 2, 3, 4, 5, 9, 10, 11, 6, 7, 8};
        // All possible on/off combinations of swap_braket, swap_bra, and
        // swap_ket
        swap_combos = {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1},
                       {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}};
        ncenters = 4;
      } break;
      case BraKet::xs_xx: {
        swap_braket_perm = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        swap_bra_perm = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        swap_ket_perm = {0, 1, 2, 6, 7, 8, 3, 4, 5};
        swap_combos = {{0, 0, 0}, {0, 0, 1}};
        ncenters = 3;
      } break;

      default:
        assert(false &&
               "Derivative index maps for this BraKet type not implemented.");
    }

    // Number of possible derivatives
    auto nderivs = how_many_derivs(
        ncenters, deriv_order);  // e.g. for 4 center: 12,78,364,1365
    // Number of differentiable parameters in a shell set (assumes 3 cartesian
    // components each center)
    auto nparams = ncenters * 3;

    // Initialize mapDerivIndex. First three dimensions are either 0 or 1
    // acting as a bool for whether to swap braket, swap bra, or swap ket.
    // Last index is the integer map.
    // Note that BraKet::xx_xx uses the whole thing, BraKet::xs_xx, only need
    // [0][0][:][:] slice and BraKet::xx_sx, only need [0][:][0][:] slice
    size_t swap_dim = 2;
    Tensor<size_t> mapDerivIndex{swap_dim, swap_dim, swap_dim, nderivs};

    // Get lookup which maps flattened upper triangle index to the
    // multidimensional index in terms of full array axes. Each axis of this
    // multidimensional array represents a different partial derivative of the
    // shell set
    vector<vector<size_t>> lookup =
        generate_multi_index_lookup(nparams, deriv_order);

    // Now loop over every combination of swap_*
    for (size_t swap = 0; swap < swap_combos.size(); swap++) {
      auto swap_braket = swap_combos[swap][0];
      auto swap_bra = swap_combos[swap][1];
      auto swap_ket = swap_combos[swap][2];
      // For every single derivative index, lookup its multidimensional index
      // and apply the permutation rules for this Braket type
      for (size_t i = 0; i < nderivs; i++) {
        vector<size_t> multi_idx = lookup[i];
        vector<size_t> permuted_indices;
        for (size_t j = 0; j < multi_idx.size(); j++) {
          size_t idx = multi_idx[j];
          if (swap_braket == 1) idx = swap_braket_perm[idx];
          if (swap_bra == 1) idx = swap_bra_perm[idx];
          if (swap_ket == 1) idx = swap_ket_perm[idx];
          permuted_indices.push_back(idx);
        }
        // Sort permuted indices into order of upper triangle indices, i <= j <=
        // k...
        sort(permuted_indices.begin(), permuted_indices.end());

        // Since the vector of vectors 'lookup' generated by
        // generate_multi_index_lookup is sorted such that each vector is
        // elementwise <= the next vector, we can use binary search to find the
        // new 1d index from the permuted multidimensional index
        size_t new_idx = 0;
        auto it = lower_bound(lookup.begin(), lookup.end(), permuted_indices);
        if (it != lookup.end()) new_idx = it - lookup.begin();
        (*mapDerivIndex.data(swap_braket, swap_bra, swap_ket, i)) = new_idx;
      }
    }
    return mapDerivIndex;
  }

};  // class definition

}  // namespace libint2

#endif /* _libint2_include_derivmap_h_*/
