
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
  class QuantumSet : public ConstructablePolymorphically,
                     public Hashable<LIBINT2_UINT_LEAST64,ComputeKey> {
  public:
    /// Quantum numbers lie in range [0,max_quantum_number)
    static const LIBINT2_UINT_LEAST64 max_quantum_number = 100;
    
    virtual ~QuantumSet() {}
    virtual const std::string label() const =0;

    /// Number of quantum numbers in the set
    virtual const unsigned int num_quanta() const =0;
    /// Increment i-th quantum number
    virtual void inc(unsigned int i) =0;
    /// Decrement i-th quantum number
    virtual void dec(unsigned int i) =0;
  };

  /** QuantumNumbers<T,N> is a set of N quantum numbers of type T implemented in terms of std::vector.
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
    QuantumNumbers(const ConstructablePolymorphically&);
    ~QuantumNumbers();
    
    bool operator==(const QuantumNumbers&) const;
    const std::string label() const;

    /// Increment quantum number i
    void inc(unsigned int i) { ++qn_.at(i); }
    /// Decrement quantum number i
    void dec(unsigned int i) {
      if (qn_.at(i) == T(0))
        throw std::runtime_error("QuantumNumbers::dec -- quantum number already zero");
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

    /// Implements Hashable::key()
    LIBINT2_UINT_LEAST64 key() const {
      LIBINT2_UINT_LEAST64 key = 0;
      LIBINT2_UINT_LEAST64 pfac = 1;
      const int maxi = ((int)num_quanta()) - 1;
      for(int i=maxi; i>=0; i--) {
        key += pfac*qn_[i];
        pfac *= QuantumSet::max_quantum_number;
      }
      return key;
    }

    /// key is in range [0,max_key())
    LIBINT2_UINT_LEAST64 max_key() const {
      LIBINT2_UINT_LEAST64 max_key = 1;
      const int maxi = ((int)num_quanta()) - 1;
      for(int i=maxi; i>=0; i--) {
        max_key *= QuantumSet::max_quantum_number;
      }
      return max_key;
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
    QuantumNumbers<T,N>::QuantumNumbers(const ConstructablePolymorphically& sptr)
    {
      const QuantumNumbers<T,N>& sptr_cast = dynamic_cast<const QuantumNumbers&>(sptr);
      qn_ = sptr_cast.qn_;
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
      std::ostringstream oss;
      oss << "{";
      if (qn_.size() > 0)
	oss << qn_[0];
      for(int i=1; i<qn_.size(); i++)
        oss << "," << qn_[i];
      oss << "}";
      return oss.str();
    }
  
  typedef QuantumNumbers<unsigned int,0> NullQuantumSet;
  
  
  /**
     QuantumNumbersA<T,N> is a set of N quantum numbers of type T implemented in terms of a C-style array.
     QuantumNumbersA is faster than QuantumNumbers but is not as safe!
  */
  template<typename T, unsigned int N> class QuantumNumbersA : public QuantumSet {

    T qn_[N];

  public:
    typedef QuantumSet parent_type;
    /// QuantumSet is a set of one QuantumSet
    typedef QuantumNumbersA iter_type;

    QuantumNumbersA(const T* qn);
    QuantumNumbersA(const vector<T>& qn);
    QuantumNumbersA(const SafePtr<QuantumNumbersA>&);
    QuantumNumbersA(const SafePtr<QuantumSet>&);
    QuantumNumbersA(const SafePtr<ConstructablePolymorphically>&);
    ~QuantumNumbersA();
    
    bool operator==(const QuantumNumbersA&) const;
    const std::string label() const;

    /// Increment quantum number i
    void inc(unsigned int i) { ++qn_[i]; }
    /// Decrement quantum number i
    void dec(unsigned int i) {
      if (qn_[i] == T(0))
        throw std::runtime_error("QuantumNumbersA::dec -- quantum number already zero");
      --qn_[i];
    }
    
    /// Return i-th quantum number
    const T elem(unsigned int i) const {
      return qn_[i];
    }

    /// Implementation of QuantumSet::num_quanta()
    const unsigned int num_quanta() const {
      return N;
    }

    /// Implements Hashable::key()
    LIBINT2_UINT_LEAST64 key() const {
      LIBINT2_UINT_LEAST64 key = 0;
      LIBINT2_UINT_LEAST64 pfac = 1;
      const int maxi = ((int)num_quanta()) - 1;
      for(int i=maxi; i>=0; i--) {
        key += pfac*qn_[i];
        pfac *= QuantumSet::max_quantum_number;
      }
      return key;
    }

    /// key is in range [0,max_key())
    LIBINT2_UINT_LEAST64 max_key() const {
      LIBINT2_UINT_LEAST64 max_key = 1;
      const int maxi = ((int)num_quanta()) - 1;
      for(int i=maxi; i>=0; i--) {
        max_key *= QuantumSet::max_quantum_number;
      }
      return max_key;
    }
    
  };

  template<typename T, unsigned int N>
    QuantumNumbersA<T,N>::QuantumNumbersA(const T* qn)
    {
      for(int i=0; i<N; i++)
        qn_[i] = qn[i];
    }

  template<typename T, unsigned int N>
    QuantumNumbersA<T,N>::QuantumNumbersA(const vector<T>& qn)
    {
      for(int i=0; i<N; i++)
        qn_[i] = qn[i];
    }

  template<typename T, unsigned int N>
    QuantumNumbersA<T,N>::QuantumNumbersA(const SafePtr<QuantumNumbersA>& sptr)
    {
      T* qn = sptr->qn_;
      for(int i=0; i<N; i++)
        qn_[i] = qn[i];
    }
  
  template<typename T, unsigned int N>
    QuantumNumbersA<T,N>::QuantumNumbersA(const SafePtr<QuantumSet>& sptr)
    {
      const SafePtr< QuantumNumbersA<T,N> > sptr_cast = dynamic_pointer_cast<QuantumNumbersA,QuantumSet>(sptr);
      if (sptr_cast == 0)
        throw std::runtime_error("QuantumNumbersA<T,N>::QuantumNumbersA(const SafePtr<QuantumSet>& sptr) -- type of sptr is incompatible with QuantumNumbersA");

      T* qn = sptr_cast->qn_;
      for(int i=0; i<N; i++)
        qn_[i] = qn[i];
    }

  template<typename T, unsigned int N>
    QuantumNumbersA<T,N>::QuantumNumbersA(const SafePtr<ConstructablePolymorphically>& sptr)
    {
      const SafePtr< QuantumNumbersA<T,N> > sptr_cast = dynamic_pointer_cast<QuantumNumbersA,ConstructablePolymorphically>(sptr);
      if (sptr_cast == 0)
        throw std::runtime_error("QuantumNumbersA<T,N>::QuantumNumbersA(const SafePtr<ConstructablePolymorphically>& sptr) -- type of sptr is incompatible with QuantumNumbersA");

      T* qn = sptr_cast->qn_;
      for(int i=0; i<N; i++)
        qn_[i] = qn[i];
    }
    
  template<typename T, unsigned int N>
    QuantumNumbersA<T,N>::~QuantumNumbersA()
    {
    }

  template<typename T, unsigned int N>
  bool
    QuantumNumbersA<T,N>::operator==(const QuantumNumbersA& a) const
    {
      const T* qn0 = qn_;
      const T* qn1 = a.qn_;
      for(int i=0; i<N; i++, ++qn0, ++qn1)
        if (*qn0 != *qn1)
          return false;

      return true;
    }

  template<typename T, unsigned int N>
    const std::string
    QuantumNumbersA<T,N>::label() const
    {
      std::ostringstream oss;
      oss << "{";
      if (N > 0)
	oss << qn_[0];
      for(int i=1; i<N; i++)
        oss << "," << qn_[i];
      oss << "}";
      return oss.str();
    }
  
  /** partial specialization of QuantumNumbersA for the case N=0 */
    /**
     QuantumNumbersA<T,N> is a set of N quantum numbers of type T implemented in terms of a C-style array.
     QuantumNumbersA is faster than QuantumNumbers but is not as safe!
  */
  template<typename T> class QuantumNumbersA<T,0> : public QuantumSet {

  public:
    typedef QuantumSet parent_type;
    /// QuantumSet is a set of one QuantumSet
    typedef QuantumNumbersA iter_type;

    QuantumNumbersA(const T* qn) {}
    QuantumNumbersA(const vector<T>& qn) {}
    QuantumNumbersA(const SafePtr<QuantumNumbersA>&) {}
    QuantumNumbersA(const SafePtr<QuantumSet>&) {}
    QuantumNumbersA(const SafePtr<ConstructablePolymorphically>&) {}
    ~QuantumNumbersA() {}
    
    bool operator==(const QuantumNumbersA&) const { return true; }
    const std::string label() const { return "{}"; }

    /// Increment quantum number i
    void inc(unsigned int i) { throw std::runtime_error("QuantumNumbersA<T,0>::inc -- no quantum numbers to increment"); }
    /// Decrement quantum number i
    void dec(unsigned int i) {
      throw std::runtime_error("QuantumNumbersA<T,0>::inc -- no quantum numbers to decrement");
    }
    /// Return i-th quantum number
    const T elem(unsigned int i) const { throw std::runtime_error("QuantumNumbersA<T,0>::inc -- no quantum numbers to return"); }
    /// Implementation of QuantumSet::num_quanta()
    const unsigned int num_quanta() const { return 0; }

    /// Implements Hashable::key()
    LIBINT2_UINT_LEAST64 key() const { return 0; }

    /// key is in range [0,max_key())
    LIBINT2_UINT_LEAST64 max_key() const { return 1; }
    
  };

};

#endif
