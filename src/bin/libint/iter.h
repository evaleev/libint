
#include <vector>
#include <smart_ptr.h>
#include <policy.h>
#include <polyconstr.h>

#ifndef _libint2_src_bin_libint_iter_h_
#define _libint2_src_bin_libint_iter_h_

// gcc 3.4 doesn't seem to allow
//#define ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES


using namespace std;

namespace libint2 {

  /** SubIterator provides a base class for all subiterator classes. Subiterators iterate through
      certain types of objects as if they were sets of some other data. For example, SubIterator can
      be implemented for iterating Gaussian functions within shells, or over integrals within shell
      sets of integrals.
  */
  class SubIterator {
  
  public:
    /// Returns a number of iterations (number of elements in a set over which to iterate).
    virtual const unsigned int num_iter() const =0;
    /// Initializes the iterator.
    virtual void init() =0;
    /// Iterates to the next element. Only prefix form is provided.
    virtual SubIterator& operator++() =0;
    /// This is used to check whether next element exists. Returns 1 if it does.
    virtual operator int() const =0;
    /** Returns pointer to current element via a pointer to base. This function can only be
      be implemented if elements are derived rom ConstructablePolymorphically. */
    virtual const SafePtr<ConstructablePolymorphically> pelem() const =0;

  protected:
    SubIterator();
    virtual ~SubIterator();

  private:
    SubIterator operator++(int);
  };

  /** SubIteratorBase<T> provides a base class for a sub-iterator class for T. It iterates through
      T as if it were a set of some data of type T::iter_type. Traits of class T (ordering of
      T::iter_type, etc.) are provided by Tr.
  */
  template <class T, class Tr = Policy<T> > class SubIteratorBase : public SubIterator {

  public:
    typedef typename T::iter_type iter_type;
    /// the only allowed constructor
    SubIteratorBase(const SafePtr<T>&);
    virtual ~SubIteratorBase();
    
    /// Returns current element
    const SafePtr<iter_type> elem() const;
    /// Returns current element. Implements SubIterator's pelem().
    const SafePtr<ConstructablePolymorphically> pelem() const;

    /// Returns a number of iterations (number of elements in a set over which to iterate).
    const unsigned int num_iter() const;
    /// Initializes the iterator.
    void init();
    /// Iterates to the next element. Only prefix form is provided.
    SubIterator& operator++();
    /// This is used to check whether current element exists. Returns 1 if it does.
    operator int() const;
    
  protected:
    const SafePtr<T> obj_;
    vector< SafePtr<iter_type> > subobj_;

  private:
    /// the iteration counter (starts at 0)
    unsigned int iter_;

    // These templates are used as a trick to make possible "partial specialization
    // of a template with multiple template params". Default implementations are not provided
    // so user must provide specialization for the case X=T
    template <class X> void init_subobj();
    template <class X> void delete_subobj();
    void init_subobj();
    void delete_subobj();

  };

  template <class T, class P>
    SubIteratorBase<T,P>::SubIteratorBase(const SafePtr<T>& obj) :
    obj_(obj), subobj_(0), iter_(0)
    {
#ifdef ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES
      init_subobj<T>();
#else
      init_subobj();
#endif
    }
  
  template <class T, class P>
    SubIteratorBase<T,P>::~SubIteratorBase()
    {
#ifdef ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES
      delete_subobj<T>();
#else
      delete_subobj();
#endif
    }

  template <class T, class P>
    const unsigned int
    SubIteratorBase<T,P>::num_iter() const
    {
      return subobj_.size();
    }
  
  template <class T, class P>
    const SafePtr<typename SubIteratorBase<T,P>::iter_type>
    SubIteratorBase<T,P>::elem() const
  {
      return subobj_.at(iter_);
  }

  template <class T, class P>
    const SafePtr<ConstructablePolymorphically>
    SubIteratorBase<T,P>::pelem() const
    {
      return dynamic_pointer_cast<ConstructablePolymorphically,iter_type>(elem());
    }
  
  template <class T, class P>
    void
    SubIteratorBase<T,P>::init()
  {
      iter_ = 0;
  }

  template <class T, class P>
    SubIterator&
    SubIteratorBase<T,P>::operator++()
  {
      ++iter_;
      return *this;
  }

  template <class T, class P>
    SubIteratorBase<T,P>::operator int() const
  {
    return (iter_ < num_iter()) ? 1 : 0;
  }

  template <class T, class P>
    void
    SubIteratorBase<T,P>::init_subobj()
  {
      P::init_subobj(obj_,subobj_);
  }

  template <class T, class P>
    void
    SubIteratorBase<T,P>::delete_subobj()
  {
      P::dealloc_subobj(subobj_);
  }
  
};

#endif
