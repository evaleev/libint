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

#ifndef _libint2_src_bin_libint_deriviter_h_
#define _libint2_src_bin_libint_deriviter_h_

#include <numeric>
#include <string>
#include <stdexcept>

namespace libint2 {

  /// Iterates over unique derivative indices
  template <unsigned int NCenters>
  struct DerivIndexIterator {
    public:
      DerivIndexIterator(unsigned int deriv_order) : deriv_order_(deriv_order) {
        assert(NCenters != 0);
        std::fill(deriv_index_, deriv_index_+NCenters*3, 0u);
        deriv_index_[0] = deriv_order_;
      }

      unsigned int range_rank() const {
        unsigned int result = 1;
        for(unsigned int d=1; d<=deriv_order_; ++d) {
          result *= (NCenters*3+d-1); result /= d;
        }
        return result;
      }

      unsigned int value(unsigned int i) const {
        assert(i < NCenters*3);
        return deriv_index_[i];
      }

      const unsigned int* values() const {
        return deriv_index_;
      }

      ///
      bool last() const {
        return last(const_cast<unsigned int*>(deriv_index_), NCenters*3);
      }
      /// will throw if last() == true
      void next() {
        next(deriv_index_, NCenters*3);
      }

    private:
      unsigned int deriv_order_;
      unsigned int deriv_index_[NCenters*3];

      static void
      first(unsigned int* deriv_index, unsigned int n) {
        assert(n != 0);
        const unsigned int deriv_order = std::accumulate(deriv_index, deriv_index+n, 0u);
        std::fill(deriv_index, deriv_index+n, 0u);
        deriv_index[0] = deriv_order;
      }
      static bool
      last(unsigned int* deriv_index, unsigned int n) {
        const unsigned int deriv_order = std::accumulate(deriv_index, deriv_index+n, 0u);
        return deriv_index[n-1] == deriv_order;
      }
      static void
      next(unsigned int* deriv_index, unsigned int n) {
        if (n == 1) return;
        if (last(deriv_index+1, n-1)) {
          assert(deriv_index[0]!=0u);
          --deriv_index[0];
          ++deriv_index[1];
          first(deriv_index+1, n-1);
        }
        else
          next(deriv_index+1, n-1);
      }
  };
}

#endif /* header guard */
