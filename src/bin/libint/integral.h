
#ifndef _libint2_src_bin_libint_integral_h_
#define _libint2_src_bin_libint_integral_h_

namespace libint2 {

  /**
     O is an operator
     TList is a typelist TList1 that consists 2 second-level typelists: TList2 for bra and TList3 for ket
     TList2 and TList3 consist of O::np types which specify basis functions for each particle in bra/ket
  */
  
  // Tags for first and second typelists
  typedef Int2Type<1> Tag1;
  typedef Int2Type<2> Tag2;
  typedef Int2Type<2> Tag3;

#if 0
  template < class O > struct IntegralImpl
  {
    BFSet* bra_[O::np];
    BFSet* ket_[O::np];
  };

  // General Integral declaration takes Operator, list of bra function types and list of key functions types
  template < class O, class BraList, class KetList > class Integral;
  
  // unrolling BraList
  template < class O, class T, class U, class KetList>
    class Integral< O, PTYPELIST_2(Tag1,T,U), KetList> :
    public Integral<O, T, KetList>,
    public Integral<O, U, KetList>,
    virtual public IntegralImpl<O>
    {
    public:
      typedef PTYPELIST_2(Tag1,T,U) BraList;

      BFSet* bra[O::np];
      BFSet* ket[O::np];

      Integral(const T& bra_fn)

      ~Integral()
        {
          // This must be a list for bra and ket
          if (Length<BraKetList>::value != 2)
            assert(false);
        }
    };

  // KetList
  template < class O, class T>
    class Integral< O, PTYPELIST_2(Tag1,T,NullType<Tag1>) > :
    public Integral<O, T>
    {
    public:
      typedef PTYPELIST_2(Tag1,T,NullType<Tag1>) KetList;
    };

  // Bra list
  template < class O, class T, class U >
    class Integral< O, PTYPELIST_2(Tag2,T,U) > :
    public Integral<O, T>,
    public Integral<O, U>
    {
    public:
      typedef PTYPELIST_2(Tag2,T,U) BraFuncList;

      T* 
      T bra[O::np];

      ~Integral()
        {
          // This must be a list for bra and ket functions
          if (Length<BraFuncList>::value != O::np)
            assert(false);
        }
    };

  // specialization only exists when TList is tagged with Tag1
  template < class O, class T>
    class Integral< O, PTYPELIST_2(Tag2,T,NullType<Tag2>) > :
    public Integral<O, T>
    {
    public:
      typedef PTYPELIST_2(Tag2,T,NullType<Tag2>) KetFuncList;

      T ket[O::np];

      ~Integral()
        {
          // This must be a list for ket functions
          if (Length<KetFuncList>::value != O::np)
            assert(false);
        }
    };

#endif

};

#endif
