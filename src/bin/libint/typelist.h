
#ifndef _libint2_src_bin_libint_typelist_h_
#define _libint2_src_bin_libint_typelist_h_


namespace libint2 {

////////////////////////////////////////////////////////////////////////////////
// class NullType
// Used as a placeholder for "no type here"
// Useful as an end marker in typelists 
////////////////////////////////////////////////////////////////////////////////

  template <class Param>
    struct NullType {
      typedef Param ParamType;
    };

////////////////////////////////////////////////////////////////////////////////
// class template Int2Type
// Converts each integral constant into a unique type
// Invocation: Int2Type<v> where v is a compile-time constant integral
// Defines 'value', an enum that evaluates to v
////////////////////////////////////////////////////////////////////////////////

    template <int v>
    struct Int2Type
    {
        enum { value = v };
    };
    
////////////////////////////////////////////////////////////////////////////////
// class template Type2Type
// Converts each type into a unique, insipid type
// Invocation Type2Type<T> where T is a type
// Defines the type OriginalType which maps back to T
////////////////////////////////////////////////////////////////////////////////

    template <typename T>
    struct Type2Type
    {
        typedef T OriginalType;
    };
    
////////////////////////////////////////////////////////////////////////////////
// class template Select
// Selects one of two types based upon a boolean constant
// Invocation: Select<flag, T, U>::Result
// where:
// flag is a compile-time boolean constant
// T and U are types
// Result evaluates to T if flag is true, and to U otherwise.
////////////////////////////////////////////////////////////////////////////////

    template <bool flag, typename T, typename U>
    struct Select
    {
        typedef T Result;
    };
    template <typename T, typename U>
    struct Select<false, T, U>
    {
        typedef U Result;
    };
    
////////////////////////////////////////////////////////////////////////////////
// class template IsSameType
// Return true iff two given types are the same
// Invocation: SameType<T, U>::value
// where:
// T and U are types
// Result evaluates to true iff U == T (types equal)
////////////////////////////////////////////////////////////////////////////////

    template <typename T, typename U>
    struct IsSameType
    {
        enum { value = false };
    };
    
    template <typename T>
    struct IsSameType<T,T>
    {
        enum { value = true };
    };

////////////////////////////////////////////////////////////////////////////////
// Helper types Small and Big - guarantee that sizeof(Small) < sizeof(Big)
////////////////////////////////////////////////////////////////////////////////

    namespace Private
    {
        template <class T, class U>
        struct ConversionHelper
        {
            typedef char Small;
            struct Big { char dummy[2]; };
            static Big   Test(...);
            static Small Test(U);
            static T MakeT();
        };
    }

////////////////////////////////////////////////////////////////////////////////
// class template Conversion
// Figures out the conversion relationships between two types
// Invocations (T and U are types):
// a) Conversion<T, U>::exists
// returns (at compile time) true if there is an implicit conversion from T
// to U (example: Derived to Base)
// b) Conversion<T, U>::exists2Way
// returns (at compile time) true if there are both conversions from T
// to U and from U to T (example: int to char and back)
// c) Conversion<T, U>::sameType
// returns (at compile time) true if T and U represent the same type
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

    template <class T, class U>
    struct Conversion
    {
        typedef Private::ConversionHelper<T, U> H;
#ifndef __MWERKS__
        enum { exists = sizeof(typename H::Small) == sizeof((H::Test(H::MakeT()))) };
#else
        enum { exists = false };
#endif
        enum { exists2Way = exists && Conversion<U, T>::exists };
        enum { sameType = false };
    };
    
    template <class T>
    struct Conversion<T, T>    
    {
        enum { exists = 1, exists2Way = 1, sameType = 1 };
    };
    
    template <class T>
    struct Conversion<void, T>    
    {
        enum { exists = 0, exists2Way = 0, sameType = 0 };
    };
    
    template <class T>
    struct Conversion<T, void>    
    {
        enum { exists = 0, exists2Way = 0, sameType = 0 };
    };
    
    template <>
    struct Conversion<void, void>    
    {
    public:
        enum { exists = 1, exists2Way = 1, sameType = 1 };
    };

////////////////////////////////////////////////////////////////////////////////
// class template SuperSubclass
// Invocation: SuperSubclass<B, D>::value where B and D are types. 
// Returns true if B is a public base of D, or if B and D are aliases of the 
// same type.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

template <class T, class U>
struct SuperSubclass
{
  enum { value = (::libint2::Conversion<const volatile U*, const volatile T*>::exists &&
                  !::libint2::Conversion<const volatile T*, const volatile void*>::sameType) };
};

////////////////////////////////////////////////////////////////////////////////
// class template SuperSubclassStrict
// Invocation: SuperSubclassStrict<B, D>::value where B and D are types. 
// Returns true if B is a public base of D.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
////////////////////////////////////////////////////////////////////////////////

template<class T,class U>
struct SuperSubclassStrict
{
  enum { value = (::libint2::Conversion<const volatile U*, const volatile T*>::exists &&
                 !::libint2::Conversion<const volatile T*, const volatile void*>::sameType &&
                 !::libint2::Conversion<const volatile T*, const volatile U*>::sameType) };
};

}   // namespace libint2

////////////////////////////////////////////////////////////////////////////////
// macro SUPERSUBCLASS
// Invocation: SUPERSUBCLASS(B, D) where B and D are types. 
// Returns true if B is a public base of D, or if B and D are aliases of the 
// same type.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
// Deprecated: Use SuperSubclass class template instead.
////////////////////////////////////////////////////////////////////////////////

#define SUPERSUBCLASS(T, U) \
    ::libint2::SuperSubclass<T,U>::value

////////////////////////////////////////////////////////////////////////////////
// macro SUPERSUBCLASS_STRICT
// Invocation: SUPERSUBCLASS(B, D) where B and D are types. 
// Returns true if B is a public base of D.
//
// Caveat: might not work if T and U are in a private inheritance hierarchy.
// Deprecated: Use SuperSubclassStrict class template instead.
////////////////////////////////////////////////////////////////////////////////

#define SUPERSUBCLASS_STRICT(T, U) \
    ::libint2::SuperSubclassStrict<T,U>::value



#define PTYPELIST_1(param, T1) ::libint2::Typelist<param, T1, ::libint2::NullType<param> >

#define PTYPELIST_2(param, T1, T2) ::libint2::Typelist<param, T1, PTYPELIST_1(param, T2) >

#define PTYPELIST_3(param, T1, T2, T3) ::libint2::Typelist<param, T1, PTYPELIST_2(param, T2, T3) >


namespace libint2 {

  template <class Param, class T, class U>
    struct Typelist
    {
      typedef Param ParamType;
      typedef T Head;
      typedef U Tail;

      /// body_ and tail_ are here to be able to create instances of Typelist objects
      Head head_;
      Tail tail_;
    };

  /**
  How to figure out the length of a typelist
   */
  template <class TList> struct Length;
  template <class Param> struct Length< NullType<Param> >
  {
    enum { value = 0 };
  };
        
  template <class Param, class T, class U>
    struct Length< Typelist<Param, T, U> >
    {
      enum { value = 1 + Length<U>::value };
    };

  /**
    Return reference to value of one of the members of a typelist object
   */
  template <class TList, unsigned int i> struct TypeAt;

  template <class Param, class T, class U>
    struct TypeAt<Typelist<Param, T, U>, 0>
  {
    typedef typename Typelist<Param, T, U>::Head Result;
  };

  template <class Param, class T, class U, unsigned int i>
    struct TypeAt<Typelist<Param, T, U>, i>
  {
    typedef typename TypeAt<U,i-1>::Result Result;
  };
  
  /**
    Return reference to value of one of the members of a typelist object
   */
  template <class TList, unsigned int i> struct ValueAt;

  template <class Param, class T, class U>
    struct ValueAt<Typelist<Param, T, U>, 0>
    {
      typedef typename Typelist<Param, T, U>::Head ResultType;
      static const ResultType& typelist_member(const Typelist<Param, T, U>& tlist)
      {
        return tlist.head_;
      };
    };

  template <class Param, class T, class U, unsigned int i>
    struct ValueAt<Typelist<Param, T, U>, i>
    {
      typedef typename TypeAt<Typelist<Param, T, U>, i>::Result ResultType;
      static const ResultType& typelist_member(const Typelist<Param, T, U>& tlist)
      {
        return ValueAt<U, i-1>::typelist_member(tlist.tail_);
      };
    };

  
}

#endif
