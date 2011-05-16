
#ifndef _libint2_src_bin_libint_hashable_h_
#define _libint2_src_bin_libint_hashable_h_

#include <libint2_intrinsic_types.h>

namespace libint2 {

  /** KeyManagePolicy defines whether to compute+cache, compute, or key is just an object */
  typedef enum { CacheKey, ComputeKey, ReferToKey} KeyManagePolicy;
  
  /// use OwnKey to figure out whether the key should be stored in Hashable
  template <KeyManagePolicy KeyManage> struct OwnKey { enum {result=false}; };
  template<> struct OwnKey<CacheKey> { enum {result=true}; };

  /** If OwnsKey is true then KeyStore<T> has the key of type T, otherwise it's empty */
  template <class T, bool HasAKey>
    struct KeyStore;
  template <class T>
    struct KeyStore<T,true> {
      mutable T value;
    };
  template <class T>
    struct KeyStore<T,false> {
    };

  /** KeyTraits<T> describes following properties of type T:
      1) how to return objects of type T

      By default, objects are returned by value.
  */
  template <typename T>
    struct KeyTraits
    {
      typedef T ReturnType;
    };
  /// std::string should be returned by const reference
  template <>
    struct KeyTraits<std::string>
    {
      typedef const std::string& ReturnType;
    };
    
  /** Objects of Hashable<T> class provide hashing function key() which computes keys of type KeyType.
      key() returns KeyTraits<KeyType>::ReturnType.
   */
  template <typename KeyType, KeyManagePolicy KeyMP>
    class Hashable
    {
    public:
      typedef typename KeyTraits<KeyType>::ReturnType KeyReturnType;
      Hashable() {}
      virtual ~Hashable() {}

      //// Computes key
      virtual KeyReturnType key() const =0;

    protected:
      KeyStore<KeyType,OwnKey<KeyMP>::result> key_;
    };

  /// FNVStringHash uses Fowler/Noll/Vo algorithm to hash a char string to a 64-bit integer
  class FNVStringHash {
  public:
    /// The type of key computed using this hash
    typedef LIBINT2_UINT_LEAST64 KeyType;
    FNVStringHash():
      hval_(hval_init) {}
    ~FNVStringHash() {}

    /// Returns 64-bit integer hash of S
    inline LIBINT2_UINT_LEAST64 hash(const std::string& S);

  private:
#if __WORDSIZE == 64
    static const LIBINT2_UINT_LEAST64 hval_init = 0xcbf29ce484222325UL;
    static const LIBINT2_UINT_LEAST64 fnv_prime64 = 0x100000001b3UL;
#else
    static const LIBINT2_UINT_LEAST64 hval_init = 0xcbf29ce484222325ULL;
    static const LIBINT2_UINT_LEAST64 fnv_prime64 = 0x100000001b3ULL;
#endif
    LIBINT2_UINT_LEAST64 hval_;
  };

  LIBINT2_UINT_LEAST64
  FNVStringHash::hash(const std::string& S) {
    const unsigned char* cS = reinterpret_cast<const unsigned char*>(S.c_str());
    const unsigned char* cptr = cS;
    while (*cptr) {
      hval_ ^= (LIBINT2_UINT_LEAST64)*cptr;
      cptr++;
      hval_ *= fnv_prime64;
    }

    return hval_;
  }
  
};

#endif
