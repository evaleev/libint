#ifndef _libint2_src_bin_libint_util_h_
#define _libint2_src_bin_libint_util_h_

#include <numeric>
#include <string>
#include <stdexcept>
#include <smart_ptr.h>
#include <util_types.h>

namespace libint2 {
  std::string to_string(FunctionPosition pos);
  
  template <class Target, class Source> SafePtr<Target> require_dynamic_cast(const SafePtr<Source>& s) {
    const SafePtr<Target> t = dynamic_pointer_cast<Target,Source>(s);
    if (t == 0)
      throw std::runtime_error("require_dynamic_cast: dynamic case failed");
    return t;
  }
  template <class Target, class Source> const Target* require_dynamic_cast(const Source* s) {
    const Target* t = dynamic_cast<Target*>(s);
    if (t == 0)
      throw std::runtime_error("require_dynamic_cast: dynamic case failed");
    return t;
  }
  
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
