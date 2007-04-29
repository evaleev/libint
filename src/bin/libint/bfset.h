
#include <iostream>
#include <string>
#include <smart_ptr.h>
#include <polyconstr.h>
#include <hashable.h>

#ifndef _libint2_src_bin_libint_bfset_h_
#define _libint2_src_bin_libint_bfset_h_

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
  */
  class IncableBFSet : public BFSet {

  public:
    virtual ~IncableBFSet() {}

    /// Increment i-th quantum number. Do nothing if i is outside the allowed range.
    virtual void inc(unsigned int i) =0;
    /// Decrements i-th quantum number. Return false and do nothing is i is outside the allowed range, else return true.
    virtual bool dec(unsigned int i) =0;
    /// Returns the norm of the quantum numbers
    virtual unsigned int norm() const =0;
    /// norm() == 0
    bool zero() const { return norm() == 0; }

  protected:
    IncableBFSet() {}

  };

  class CGF;
  
  /// Cartesian Gaussian Shell
  class CGShell : public IncableBFSet, public Hashable<unsigned,ReferToKey> {

    unsigned int qn_[1];

  public:
    /// As far as SetIterator is concerned, CGShell is a set of one CGF
    typedef CGF iter_type;
    typedef IncableBFSet parent_type;

    /// Default constructor creates an s-type shell
    CGShell();
    CGShell(unsigned int qn);
    CGShell(unsigned int qn[1]);
    CGShell(const CGShell&);
    CGShell(const SafePtr<CGShell>&);
    CGShell(const SafePtr<parent_type>&);
    ~CGShell();
    CGShell& operator=(const CGShell&);

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set
    unsigned int num_bf() const { return (qn_[0]+1)*(qn_[0]+2)/2; };

    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz=0) const { return qn_[0]; }

    /// Comparison operator
    bool operator==(const CGShell&) const;

    /// Implements purely virtual IncableBFSet::dec
    bool dec(unsigned int i);
    /// Implements purely virtual IncableBFSet::inc
    void inc(unsigned int i);
    /// Implements IncableBFSet::norm()
    unsigned int norm() const;
    /// Implements Hashable<unsigned>::key()
    unsigned key() const { return qn_[0]; }
    /// The range of keys is [0,max_key]
    const static unsigned max_key = 19;
    //const static unsigned max_key = LIBINT_MAX_AM;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  };

  /// Cartesian Gaussian Function
  class CGF : public IncableBFSet, public Hashable<unsigned,ComputeKey> {

    unsigned int qn_[3];

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
    CGF(const SafePtr<CGF>&);
    CGF(const SafePtr<parent_type>&);
    CGF(const ConstructablePolymorphically&);
    CGF(const SafePtr<ConstructablePolymorphically>&);
    ~CGF();

    /// Return a compact label
    const std::string label() const;
    /// Returns the number of basis functions in the set (always 1)
    unsigned int num_bf() const { return 1; };

    /// Returns the angular momentum
    unsigned int qn(unsigned int xyz) const;

    /// Comparison operator
    bool operator==(const CGF&) const;
    
    /// Implements purely virtual IncableBFSet::dec
    bool dec(unsigned int i);
    /// Implements purely virtual IncableBFSet::inc
    void inc(unsigned int i);
    /// Implements IncableBFSet::norm()
    unsigned int norm() const;
    /// Implements Hashable<unsigned>::key()
    unsigned key() const {
      unsigned nxy = qn_[1] + qn_[2];
      unsigned l = nxy + qn_[0];
      unsigned key = nxy*(nxy+1)/2 + qn_[2];
      return key + key_l_offset[l];
    }
    /// The range of keys is [0,max_key). The formula is easily derived by summing (L+1)(L+2)/2 up to CGShell::max_key
    const static unsigned max_key = 1 + CGShell::max_key*CGShell::max_key + CGShell::max_key*(CGShell::max_key*CGShell::max_key + 11*CGShell::max_key)/6;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

  private:
    /// key_l_offset[L] is the number of all possible CGFs with angular momentum less than L
    static unsigned key_l_offset[CGShell::max_key+2];
  };

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

