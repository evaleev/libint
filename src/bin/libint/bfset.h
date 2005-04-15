
#include <iostream>
#include <string>
#include <smart_ptr.h>
#include <polyconstr.h>

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

    /// Increment i-th quantum number. Do nothing if i is outside the allowed range
    virtual void inc(unsigned int i) throw() =0;
    /// Decrements i-th quantum number. Do nothing is i is outside the allowed range
    virtual void dec(unsigned int i) =0;
    /// Returns true if all quanta are 0
    virtual bool zero() const =0;
    
  protected:
    IncableBFSet() {}

  };

  /// Cartesian Gaussian Function
  class CGF : public IncableBFSet {

    unsigned int qn_[3];

  public:
    /// As far as SetIterator is concerned, CGF is a set of one CGF
    typedef CGF iter_type;
    typedef IncableBFSet parent_type;

    /// Default constructor makes an s-type Gaussian
    CGF();
    CGF(unsigned int qn[3]);
    CGF(const CGF&);
    CGF(const SafePtr<CGF>&);
    CGF(const SafePtr<parent_type>&);
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
    
    /// Implements purely virtual IncableBFSet::dec, may throw InvalidDecrement
    void dec(unsigned int i);
    /// Implements purely virtual IncableBFSet::inc
    void inc(unsigned int i) throw();
    /// Implements IncableBFSet::zero()
    bool zero() const;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;
    
  };

  /// Cartesian Gaussian Shell
  class CGShell : public IncableBFSet {

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
    CGShell(const SafePtr<ConstructablePolymorphically>&);
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

    /// Implements purely virtual IncableBFSet::dec, may throw InvalidDecrement
    void dec(unsigned int i);
    /// Implements purely virtual IncableBFSet::inc
    void inc(unsigned int i) throw();
    /// Implements IncableBFSet::zero()
    bool zero() const;

    /// Print out the content
    void print(std::ostream& os = std::cout) const;

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

