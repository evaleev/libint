
#include <vector>
#include <typelist.h>
#include <smart_ptr.h>

#ifndef _libint2_src_bin_libint_quanta_h_
#define _libint2_src_bin_libint_quanta_h_

using namespace std;


namespace libint2 {

  /** QuantumSet is the base class for all (sets of) quantum numbers.
      QuantumSet's must be constructable using
      SafePtr<QuantumSet> or SafePtr<ConstructablePolymorphically>.
  */
  class QuantumSet : public ConstructablePolymorphically {
  public:
    virtual ~QuantumSet() {}
    
    virtual const std::string label() const =0;

    /// Number of quantum numbers in the set
    virtual const unsigned int num_quanta() const =0;
    /// Increment i-th quantum number
    virtual void inc(unsigned int i) =0;
    /// Decrement i-th quantum number
    virtual void dec(unsigned int i) =0;
  };

  /** QuantumNumbers<T,N> is a set of N quantum numbers of type T.
  */
  template<typename T, unsigned int N> class QuantumNumbers : public QuantumSet {

    vector<T> qn_;

  public:
    typedef QuantumSet parent_type;
    /// QuantumSet is a set of one QuantumSet
    typedef QuantumNumbers iter_type;

    QuantumNumbers(const vector<T>& qn);
    QuantumNumbers(const SafePtr<QuantumNumbers>&);
    QuantumNumbers(const SafePtr<QuantumSet>&);
    QuantumNumbers(const SafePtr<ConstructablePolymorphically>&);
    ~QuantumNumbers();
    
    bool operator==(const QuantumNumbers&) const;
    const std::string label() const;

    /// Increment quantum number i
    void inc(unsigned int i) { ++qn_.at(i); }
    /// Decrement quantum number i
    void dec(unsigned int i) {
      if (qn_.at(i) == T(0))
        throw std::runtime_error("QuantumNumber::dec -- quantum number already zero");
      --qn_.at(i);
    }
    
    /// Return i-th quantum number
    const T elem(unsigned int i) const {
      return qn_.at(i);
    }

    /// Return i-th quantum number
    const unsigned int num_quanta() const {
      return qn_.size();
    }

  };

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::QuantumNumbers(const vector<T>& qn) :
    qn_(qn)
    {
    }

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::QuantumNumbers(const SafePtr<QuantumNumbers>& sptr) :
    qn_(sptr->qn_)
    {
    }
  
  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::QuantumNumbers(const SafePtr<QuantumSet>& sptr)
    {
      const SafePtr< QuantumNumbers<T,N> > sptr_cast = dynamic_pointer_cast<QuantumNumbers,QuantumSet>(sptr);
      if (sptr_cast == 0)
        throw std::runtime_error("QuantumNumbers<T,N>::QuantumNumbers(const SafePtr<QuantumSet>& sptr) -- type of sptr is incompatible with QuantumNumbers");

      qn_ = sptr_cast->qn_;
    }

  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::QuantumNumbers(const SafePtr<ConstructablePolymorphically>& sptr)
    {
      const SafePtr< QuantumNumbers<T,N> > sptr_cast = dynamic_pointer_cast<QuantumNumbers,ConstructablePolymorphically>(sptr);
      if (sptr_cast == 0)
        throw std::runtime_error("QuantumNumbers<T,N>::QuantumNumbers(const SafePtr<ConstructablePolymorphically>& sptr) -- type of sptr is incompatible with QuantumNumbers");

      qn_ = sptr_cast->qn_;
    }
    
  template<typename T, unsigned int N>
    QuantumNumbers<T,N>::~QuantumNumbers()
    {
    }

  template<typename T, unsigned int N>
  bool
    QuantumNumbers<T,N>::operator==(const QuantumNumbers& a) const
    {
      return qn_ == a.qn_;
    }

  template<typename T, unsigned int N>
    const std::string
    QuantumNumbers<T,N>::label() const
    {
      std::string result = " ";
      for(int i=0; i<qn_.size(); i++)
        result += qn_[i];
      return result;
    }
  
  typedef QuantumNumbers<unsigned int,0> NullQuantumSet;
  
};

#endif
