#ifndef _libint2_src_bin_libint_boostutil_h_
#define _libint2_src_bin_libint_boostutil_h_

#include <boost/tuple/tuple.hpp>

namespace libint2 {

  /// vector of N elements of type T
  template <typename T, int N>
  class VectorN {
    public:
    /// Default is vector of zeroes
    VectorN() { for(int i=0; i<N; ++i) data_[i] = T(); }
    VectorN(const VectorN& a) { for(int i=0; i<N; ++i) data_[i] = a.data_[i]; }

    const VectorN& operator+=(const VectorN& a) {
      for(int i=0; i<N; ++i)
        data_[i] += a.data_[i];
    }
    const VectorN& operator-=(const VectorN& a) {
      for(int i=0; i<N; ++i)
        data_[i] -= a.data_[i];
    }
    
    /// 1-norm
    T norm1() const { T result(0); for(int i=0; i<N; ++i) result += abs(data_[i]); return result; }
    
    T& operator[](int i) {
      assert(i>=0 && i<N);
      return data_[i];
    }
    const T& operator[](int i) const {
      assert(i>=0 && i<N);
      return data_[i];
    }
    private:
      T data_[N];
  };
  
  template<typename T, int N>
  VectorN<T,N> operator+(const VectorN<T,N>& a,
                         const VectorN<T,N>& b) {
    VectorN<T,N> result(a);
    result += b;
    return result;
  }
  template<typename T, int N>
  VectorN<T,N> operator-(const VectorN<T,N>& a,
                         const VectorN<T,N>& b) {
    VectorN<T,N> result(a);
    result -= b;
    return result;
  }
  
  template <typename T, int N>
  inline VectorN<T,N> unit_vector(int i) {
    assert(i >= 0 && i < N);
    VectorN<T,N> result;
    result[i] += 1;
    return result;
  }
  
  // useful typedefs
  typedef VectorN<int,3> IntVec3;
  inline IntVec3 unit_intvec3(int i) {
    return unit_vector<int,3>(i);
  }
  /// return true if has elements < 0
  inline bool ltzero(const IntVec3& a) {
    for(int xyz=0; xyz<3; ++xyz)
      if (a[xyz] < 0)
        return true;
    return false;
  }
  
}

#endif /* header guard */
