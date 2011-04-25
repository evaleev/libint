
#ifndef _libint2_src_bin_libint_bfset_h_
#define _libint2_src_bin_libint_bfset_h_

#include <iostream>
#include <string>
#include <smart_ptr.h>
#include <polyconstr.h>
#include <hashable.h>
#include <contractable.h>
#include <global_macros.h>

namespace libint2 {

  /** Set of basis functions. Sets must be constructable using
      SafePtr<BFSet> or SafePtr<ConstructablePolymorphically>.
  */
  class BFSet : public ConstructablePolymorphically {

  public:
    virtual ~BFSet() {}
    virtual unsigned int num_bf() const =0;
    virtual const std::string label() const =0;

  protected:
    BFSet() {}

  };

  /** Set of basis functions with incrementable/decrementable quantum numbers.
      Sets must be constructable using SafePtr<BFSet> or SafePtr<ConstructablePolymorphically>.

      Call to dec() may invalidate the object. No further
      modification of such object's state is possible.
      Assignment of such object to another preserves the invalid state,
      but copy construction of a valid object is possible.
  */
  class IncableBFSet : public BFSet {

  public:
    virtual ~IncableBFSet() {}

    /// Add c quanta along xyz.
    virtual void inc(unsigned int xyz, unsigned int c = 1u) =0;
    /// Subtract c quanta along xyz. If impossible, invalidate the object, but do not change its quanta!
    virtual void dec(unsigned int xyz, unsigned int c = 1u) =0;
    /// Returns the norm of the quantum numbers
    virtual unsigned int norm() const =0;
    /// norm() == 0
    bool zero() const { return norm() == 0; }
    /// Return false if this object is invalid
    bool valid() const { return valid_; }

  protected:
    IncableBFSet() : valid_(true) {}

    /// make this object invalid
    void invalidate() { valid_ = false; }

  private:
    bool valid_;
  };

  /// BF with unit quantum along X. F must behave like IncableBFSet
  template <typename F>
  F unit(unsigned int X) { F tmp; tmp.inc(X,1); return tmp; }
  /// Return true if A is valid
  inline bool exists(const IncableBFSet& A) { return A.valid(); }

  /** Represents cartesian derivatives of atom-centered basis functions
  */
  class OriginDerivative : public Hashable<unsigned,ReferToKey> {

  public:
    OriginDerivative() : valid_(true) {
      d_[0] = d_[1] = d_[2] = 0u;
    }
    OriginDerivative(const OriginDerivative& other) : valid_(true) {
      std::copy(other.d_, other.d_ + 3, d_);
    }
    OriginDerivative& operator=(const OriginDerivative& other) {
      valid_ = other.valid_;
      std::copy(other.d_, other.d_ + 3, d_);
      return *this;
    }
    OriginDerivative& operator+=(const OriginDerivative& other) {
      assert(valid_);
      for(int xyz=0; xyz<3; ++xyz)
        d_[xyz] += other.d_[xyz];
      return *this;
    }
    OriginDerivative& operator-=(const OriginDerivative& other) {
      assert(valid_);
      for(int xyz=0; xyz<3; ++xyz)
        d_[xyz] -= other.d_[xyz];
      return *this;
    }

    /// returns the number of quanta along xyz
    unsigned int d(unsigned int xyz) const {
      return d_[xyz];
    }
    /// Add c quanta along xyz.
    void inc(unsigned int xyz, unsigned int c = 1u) {
      assert(valid_);
      d_[xyz] += c;
    }
    /// Subtract c quanta along xyz. If impossible, invalidate the object, but do not change its quanta!
    void dec(unsigned int xyz, unsigned int c = 1u) {
      //assert(valid_);
      if (d_[xyz] >= c)
        d_[xyz] -= c;
      else
        valid_ = false;
    }
    /// Returns the norm of the quantum numbers
    unsigned int norm() const {
      return d_[0] + d_[1] + d_[2];
    }
    /// norm() == 0
    bool zero() const { return norm() == 0; }
    /// Return false if this object is invalid
    bool valid() const { return valid_; }
    /// Implements Hashable<unsigned>::key()
    unsigned key() const {
      unsigned nxy = d_[1] + d_[2];
      unsigned l = nxy + d_[0];
      unsigned key = nxy*(nxy+1)/2 + d_[2];
      return key + key_l_offset[l];
    }
    /// Return a compact label
    const std::string label() const {
      char result[] = "000";
      for(unsigned int xyz=0; xyz<3; ++xyz)
        result[xyz] += d_[xyz];
      return std::string(result);
    }

    const static unsigned max_deriv = 4;
    /// The range of keys is [0,max_key). The formula is easily derived by summing (L+1)(L+2)/2 up to max_deriv
    const static unsigned max_key = (1 + max_deriv)*(2 + max_deriv)*(3 + max_deriv)/6;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  protected:

    /// make this object invalid
    void invalidate() { valid_ = false; }

  private:
    unsigned int d_[3];
    bool valid_;  // indicates valid/invalid state
    /// key_l_offset[L] is the number of all possible derivatives of order up to L
    static unsigned key_l_offset[max_deriv+1];

  };

  OriginDerivative operator-(const OriginDerivative& A, const OriginDerivative& B);
  bool operator==(const OriginDerivative& A, const OriginDerivative& B);
  /// Return true if A is valid
  inline bool exists(const OriginDerivative& A) { return A.valid(); }

  class CGF;

  /// Cartesian Gaussian Shell
  class CGShell : public IncableBFSet, public Hashable<unsigned,ReferToKey>,
                  public Contractable<CGShell> {

    unsigned int qn_[1];
    OriginDerivative deriv_;

    friend CGShell operator+(const CGShell& A, const CGShell& B);
    friend CGShell operator-(const CGShell& A, const CGShell& B);

  public:
    /// As far as SetIterator is concerned, CGShell is a set of one CGF
    typedef CGF iter_type;
    typedef IncableBFSet parent_type;

    /// Default constructor creates an s-type shell
    CGShell();
    CGShell(unsigned int qn);
    CGShell(const CGShell&);
    virtual ~CGShell();
    CGShell& operator=(const CGShell&);

    const OriginDerivative& deriv() const { return deriv_; }
    OriginDerivative& deriv() { return deriv_; }

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set
    unsigned int num_bf() const { return (qn_[0]+1)*(qn_[0]+2)/2; };
    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz=0) const { return qn_[0]; }

    /// Comparison operator
    bool operator==(const CGShell&) const;

    /// Implementation of IncableBFSet::inc().
    void inc(unsigned int xyz, unsigned int c = 1u);
    /// Implementation of IncableBFSet::dec().
    void dec(unsigned int xyz, unsigned int c = 1u);
    /// Implements IncableBFSet::norm()
    unsigned int norm() const;
    /// Implements Hashable<unsigned>::key()
    unsigned key() const { return (deriv().key() * 2 + (contracted() ? 1 : 0)) * (max_qn+1) + qn_[0]; }
    const static unsigned max_qn = LIBINT_CARTGAUSS_MAX_AM;
    // The range of keys is [0,max_key]
    const static unsigned max_key = 2 * (max_qn + 1) * OriginDerivative::max_key;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  };

  CGShell operator+(const CGShell& A, const CGShell& B);
  CGShell operator-(const CGShell& A, const CGShell& B);

  /// Cartesian Gaussian Function
  class CGF : public IncableBFSet, public Hashable<unsigned,ComputeKey>,
              public Contractable<CGF> {

    unsigned int qn_[3];
    OriginDerivative deriv_;

    friend CGF operator+(const CGF& A, const CGF& B);
    friend CGF operator-(const CGF& A, const CGF& B);

  public:
    /// As far as SetIterator is concerned, CGF is a set of one CGF
    typedef CGF iter_type;
    typedef IncableBFSet parent_type;
    /// How to return key
    //typedef typename Hashable<unsigned,ComputeKey>::KeyReturnType KeyReturnType;

    /// Default constructor makes an s-type Gaussian
    CGF();
    CGF(unsigned int qn[3]);
    CGF(const CGF&);
    CGF(const ConstructablePolymorphically&);
    virtual ~CGF();
    /// assignment
    CGF& operator=(const CGF&);

    const OriginDerivative& deriv() const { return deriv_; }
    OriginDerivative& deriv() { return deriv_; }

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set (always 1)
    unsigned int num_bf() const { return 1; };
    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz) const;

    /// Comparison operator
    bool operator==(const CGF&) const;

    /// Implementation of IncableBFSet::inc().
    void inc(unsigned int xyz, unsigned int c = 1u);
    /// Implementation of IncableBFSet::dec().
    void dec(unsigned int xyz, unsigned int c = 1u);
    /// Implements IncableBFSet::norm()
    unsigned int norm() const;
    /// Implements Hashable<unsigned>::key()
    unsigned key() const {
      unsigned nxy = qn_[1] + qn_[2];
      unsigned l = nxy + qn_[0];
      unsigned key = nxy*(nxy+1)/2 + qn_[2];
      return ( deriv().key() * 2 + (contracted() ? 1 : 0)) * max_num_qn + key + key_l_offset[l];
    }
    /// The range of keys is [0,max_key). The formula is easily derived by summing (L+1)(L+2)/2 up to CGShell::max_key
    /// The factor of 2 to account for contracted vs. uncontracted basis functions
    /// The factor of OriginDerivative::max_key to account for derivatives
    const static unsigned max_num_qn = ((1 + (CGShell::max_qn+1)) * (2 + (CGShell::max_qn+1)) * (3 + (CGShell::max_qn+1)) /6);
    const static unsigned max_key = 2 * OriginDerivative::max_key * max_num_qn;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  private:
    /// key_l_offset[L] is the number of all possible CGFs with angular momentum less than L
    static unsigned key_l_offset[CGShell::max_key+2];
  };

  CGF operator+(const CGF& A, const CGF& B);
  CGF operator-(const CGF& A, const CGF& B);

  /**
     TrivialBFSet<T> defines static member result, which is true if T
     is a basis function set consisting of 1 function
  */
  template <class T>
    struct TrivialBFSet;
  template <>
    struct TrivialBFSet<CGShell> {
      static const bool result = false;
    };
  template <>
    struct TrivialBFSet<CGF> {
      static const bool result = true;
    };

};

#endif

