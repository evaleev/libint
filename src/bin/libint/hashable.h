
#include <libint2_types.h>

#ifndef _libint2_src_bin_libint_hashable_h_
#define _libint2_src_bin_libint_hashable_h_

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
  */
  template <typename T>
    struct KeyTraits;
  template <>
    struct KeyTraits<unsigned>
    {
      typedef unsigned ReturnType;
    };
  template <>
    struct KeyTraits<LIBINT2_UINT_LEAST64>
    {
      typedef LIBINT2_UINT_LEAST64 ReturnType;
    };
  template <>
    struct KeyTraits<std::string>
    {
      typedef const std::string& ReturnType;
    };
  template <>
    struct KeyTraits<double>
    {
      typedef const double& ReturnType;
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
      ~Hashable() {}

      //// Computes key
      virtual KeyReturnType key() const =0;

    protected:
      KeyStore<KeyType,OwnKey<KeyMP>::result> key_;
    };

};

#endif
