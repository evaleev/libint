
#ifndef _libint2_src_bin_libint_contract_h_
#define _libint2_src_bin_libint_contract_h_

namespace libint2 {

  /// use this as a base to add to Derived a "contracted()" attribute
  template <typename Derived> class Contractable {
    public:
      Contractable() : value_(default_value_) {}
      Contractable(const Contractable& source) : value_(source.value_) {}
      Contractable& operator=(const Contractable& source) {
        value_ = source.value_;
        return *this;
      }
      bool contracted() const { return value_; }
      void uncontract() { value_ = false; }
      static void set_contracted_default_value(bool dv) { default_value_ = dv; }
    private:
      bool value_;
      static bool default_value_;
  };
  template <typename Derived>
  bool Contractable<Derived>::default_value_ = false;

};

#endif
