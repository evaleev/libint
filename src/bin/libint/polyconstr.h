

#ifndef _libint2_src_bin_libint_polyconstr_h_
#define _libint2_src_bin_libint_polyconstr_h_


namespace libint2 {

  /** ConstructablePolymorphically is a base for all objects
      which can be constructed using a SafePtr to a base or a
      SafePtr to ConstructablePolymorphically.
  */
  class ConstructablePolymorphically {
  protected:
    ConstructablePolymorphically() {}
    virtual ~ConstructablePolymorphically() {}
  };

};

#endif
