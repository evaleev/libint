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


#ifndef _libint2_include_derivmap_h_
#define _libint2_include_derivmap_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "The simple Libint API requires C++11 support"
#endif

#include <stdlib.h>
#include <vector>
#include <algorithm>

namespace libint2 {

  namespace derivmap {

  // Singleton class to statically initialize all mapDerivIndex arrays up to 
  // the order that the library supports     
    class DerivMapGenerator { 
      public:
        // TODO here generate all arrays statically up to LIBINT2_MAX_DERIV_ORDER

      private:
        // Functions required for generating the maps.
        // Combinations with repitition. 
        // Combinations of size 'k' from 'n' elements stored in vector 'inp'.
        // Requires instantiating vector 'out' and vector of vectors 'result', which stores every combination.
        void cwr_recursion(std::vector<int> inp,
                                  std::vector<int> &out,
                                  std::vector<std::vector<int>> &result,
                                  int k, int i, int n)
        {   
            // base case: if combination size is k, add to result 
            if (out.size() == k){
                result.push_back(out);
                return;
            }
            for (int j = i; j < n; j++){
                out.push_back(inp[j]);
                cwr_recursion(inp, out, result, k, j, n);
                // backtrack - remove current element from solution
                out.pop_back();
            }
        }

        // How many shell set derivatives
        // k is number of centers, n is deriv order.
        // Assumes coordinate-independent operator, 3 coordinates per center
        int how_many_derivs(int k, int n) {
            int val = 1;
            int factorial = 1;
            for (int i=0; i < n; i++) {
                val *= (3 * k + i); 
                factorial *= i + 1;
            }   
            val /= factorial;
            return val;
        }

        // Cartesian product of a series of vectors
        std::vector<std::vector<int>> cartesian_product (const std::vector<std::vector<int>>& v) {
            std::vector<std::vector<int>> s = {{}};
            for (const auto& u : v) {
                std::vector<std::vector<int>> r;
                for (const auto& x : s) {
                    for (const auto y : u) {
                        r.push_back(x);
                        r.back().push_back(y);
                    }
                }
                s = std::move(r);
            }   
            return s;
        }

        // Create array which is of size (nderivs, deriv_order)
        // which when indexed with a 1d buffer index (the flattened generalized upper triangle index)
        // it maps it to the corresponding to multidimensional index tuple (i,j,k,...) where i<=j<=k<=...
        // The multidimensional indices each correspond to a partial derivative operator, 
        // and the value of each index indicates which parameter is being differentiated.
        // For example, if there are 6 differentiable parameters, and deriv_order = 2,
        // calling generate_multi_index_lookup(6,2) yields 'combos' such that combos[13] yields the pair [2,4]
        // which matches the value and index pair in the below array:
        // 0  1  2  3  4  5
        //    6  7  8  9 10
        //      11 12 13 14
        //         15 16 17
        //            18 19
        //               20
        std::vector<std::vector<int>> generate_multi_index_lookup(int nparams, int deriv_order) {
            using namespace std;
            // Generate vector of indices 0 through nparams-1
            vector<int> inp;
            for (int i = 0; i < nparams; i++) {
                inp.push_back(i);
            }
            // Generate all possible combinations with repitition. 
            // These are upper triangle indices, and the length of them is the total number of derivatives
            vector<int> out;
            vector<vector<int>> combos;
            cwr_recursion(inp, out, combos, deriv_order, 0, nparams);
            return combos;
        }

        // Function which computes mapDerivIndex 
        std::vector<std::vector<std::vector<std::vector<int>>>> 
        generate_deriv_index_map(int deriv_order, BraKet braket_)
        {
            using namespace std;
            // Number of possible derivatives
            int nderivs = how_many_derivs(ncenters, deriv_order); // e.g. for 4 center: 12,78,364,1365
            // Number of differentiable parameters in a shell set (assumes 3 cartesian components each center)
            int nparams = ncenters * 3;
        
            // The BraKet type determines what the permutations are.
            vector<int> swap_braket_perm;
            vector<int> swap_bra_perm;
            vector<int> swap_ket_perm;
            vector<vector<int>> swap_combos;
            // TODO switch to checking BraKet types
            if (ncenters == 4){
                swap_braket_perm = {6,7,8,9,10,11,0,1,2,3,4,5};
                swap_bra_perm = {3,4,5,0,1,2,6,7,8,9,10,11};
                swap_ket_perm = {0,1,2,3,4,5,9,10,11,6,7,8};
                // All possible on/off combinations of swap_braket, swap_bra, and swap_ket
                swap_combos = {{0, 0, 0},
                               {0, 0, 1},
                               {0, 1, 0},
                               {0, 1, 1},
                               {1, 0, 0},
                               {1, 0, 1},
                               {1, 1, 0},
                               {1, 1, 1}};
            }
            // TODO switch to checking BraKet types
            // TODO should I add BraKet::xx_xs?
            if (ncenters == 3){
                swap_braket_perm = {0,1,2,3,4,5,6,7,8};
                swap_bra_perm = {0,1,2,3,4,5,6,7,8};
                swap_ket_perm = {0,1,2,6,7,8,3,4,5};
                swap_combos = {{0,0,0}, {0,0,1}};
            }
            // Initialize mapDerivIndex. First three dimensions are either 0 or 1
            // acting as a bool for whether to swap braket, swap bra, or swap ket.
            // Last index is the integer map.
            // Note that BraKet::xx_xx uses the whole thing, BraKet::xs_xx, only need [0][0][:][:] slice
            // and BraKet::xx_sx, only need [0][:][0][:] slice
            vector<vector<vector<vector<int>>>> mapDerivIndex(2, vector<vector<vector<int>>>
                                                             (2, vector<vector<int>>
                                                             (2, vector<int>
                                                             (nderivs, 0))));

            // Get lookup which maps flattened upper triangle index to the multidimensional index 
            // in terms of full array axes. Each axis of this multidimensional array represents
            // a different partial derivative of the shell set
            vector<vector<int>> lookup = generate_multi_index_lookup(nparams, deriv_order);
        
            // Now loop over every combination of swap_* 
            for (int swap = 0; swap < swap_combos.size(); swap++){
                int swap_braket = swap_combos[swap][0];
                int swap_bra = swap_combos[swap][1];
                int swap_ket = swap_combos[swap][2];
                // For every single derivative index, lookup its multidimensional index
                // and apply the permutation rules for this Braket::*
                for (int i = 0; i < nderivs; i++){
                    vector<int> multi_idx = lookup[i];
                    vector<int> permuted_indices;
                    for (int j=0; j < multi_idx.size(); j++){
                        int idx = multi_idx[j];
                        if (swap_braket == 1) idx = swap_braket_perm[idx];
                        if (swap_bra == 1) idx = swap_bra_perm[idx];
                        if (swap_ket == 1) idx = swap_ket_perm[idx];
                        permuted_indices.push_back(idx);
                    }
                    // Sort permuted indices into order of upper triangle indices, i <= j <= k...
                    sort(permuted_indices.begin(), permuted_indices.end());

                    // Now we need to map back to the 1d index.
                    // There are two easy ways to do this.
                    // The first is to construct a (deriv_order)-dimensional array, where each element is a 1d buffer index.
                    // Then we can just index the array with 'permuted_indices' to find the 1d index
                    // This is very fast, but has large memory requirement at higher orders for storing the array.
                    // Instead, we can loop through 'lookup' until 'permuted_indices' matches an entry, and then that entry is the desired 1d buffer index.
                    // This requires many more loops, but is cheaper memory wise.
                    int new_idx = 0;
                    auto it = lower_bound(lookup.begin(), lookup.end(), permuted_indices);
                    if (it != lookup.end()) new_idx = it - lookup.begin();
                    mapDerivIndex[swap_braket][swap_bra][swap_ket][i] = new_idx;
                }
            }
            return mapDerivIndex;
        }

    } // class definition

  } // namespace libint2::derivmap

} // namespace libint2

#endif  /* _libint2_include_solidharmonics_h_ */
