
#include <vector>
#include <smart_ptr.h>
#include <policy.h>
#include <polyconstr.h>
#include <exception.h>

#ifndef _libint2_src_bin_libint_iter_h_
#define _libint2_src_bin_libint_iter_h_

// gcc 3.4 doesn't seem to allow
//#define ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES


using namespace std;

namespace libint2 {

  struct DummyIterator;
  
  /** Iterator provides a base class for all object iterator classes. It iterates over
      certain objects as if they were sets of some other data. For example, Iterator can
      be implemented for iterating over Gaussian functions within shells, or over integrals within shell
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
    /** Return current element via base class. These functions can only be
      be implemented if elements are derived from ConstructablePolymorphically.
      Default implementation throws, thus must be overridden by SubIteratorBase<T>. */
    virtual const ConstructablePolymorphically& pelem() const;
    virtual ~SubIterator();

  protected:
    SubIterator();

  private:
    //SubIterator operator++(int);
  };

  /** SubIteratorBase<T> provides a base class for a sub-iterator class for T. It iterates through
      T as if it were a set of some data of type T::iter_type. Traits of class T (ordering of
      T::iter_type, etc.) are provided by Tr<T>.
  */
  template <class T, template <class> class Tr = Policy > class SubIteratorBase : public SubIterator {

  public:
    typedef typename T::iter_type iter_type;
    typedef Tr<T> TPolicy;
    typedef typename TPolicy::obj_stype tref;
    typedef typename TPolicy::subobj_stype iref;
    /// Return reference to ConstructablePolymorphically as object of this type
    typedef const ConstructablePolymorphically& cp_rettype;
    SubIteratorBase(const tref&);
    virtual ~SubIteratorBase();
    
    /// Returns current element
    const iref& elem() const;
    /// Returns current element. Implements SubIterator's pelem().
    cp_rettype pelem() const;

    /// Returns a number of iterations (number of elements in a set over which to iterate).
    const unsigned int num_iter() const;
    /// Initializes the iterator.
    void init();
    /// Iterates to the next element. Only prefix form is provided.
    SubIterator& operator++();
    /// This is used to check whether current element exists. Returns 1 if it does.
    operator int() const;
    
  protected:
    const tref obj_;
    vector<iref> subobj_;

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

    // implementation of pelem()
    template <typename X, bool return_smart_ptr>
    struct PElemImpl {
    };
    template <typename X>
    struct PElemImpl<X,true> {
      static cp_rettype pelem(const iref& elem) {
        SafePtr<ConstructablePolymorphically> elem_cast = dynamic_pointer_cast<ConstructablePolymorphically,X>(elem);
        return *(elem_cast.get());
      }
    };
    template <typename X>
    struct PElemImpl<X,false> {
      static cp_rettype pelem(const iref& elem) {
        return elem;
      }
    };
    
  };

  template <class T, template <class> class P>
    SubIteratorBase<T,P>::SubIteratorBase(const tref& obj) :
    obj_(obj), subobj_(0), iter_(0)
    {
#ifdef ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES
      init_subobj<T>();
#else
      init_subobj();
#endif
    }
  
  template <class T, template <class> class P>
    SubIteratorBase<T,P>::~SubIteratorBase()
    {
#ifdef ALLOW_PARTIALLY_SPECIALIZED_NESTED_TEMPLATES
      delete_subobj<T>();
#else
      delete_subobj();
#endif
    }

  template <class T, template <class> class P>
    const unsigned int
    SubIteratorBase<T,P>::num_iter() const
    {
      return subobj_.size();
    }
  
  template <class T, template <class> class P>
    const typename SubIteratorBase<T,P>::iref&
    SubIteratorBase<T,P>::elem() const
    {
      return subobj_.at(iter_);
    }

  template <class T, template <class> class P>
    typename SubIteratorBase<T,P>::cp_rettype
    SubIteratorBase<T,P>::pelem() const
    {
      return PElemImpl<iter_type,IsSafePtr<iref>::result>::pelem(elem());
    }

#if 0
  template <class T, template <class> class P>
    const SafePtr<ConstructablePolymorphically>
    SubIteratorBase<T,P>::pelem() const
    {
      return dynamic_pointer_cast<ConstructablePolymorphically,iter_type>(elem());
    }
#endif
  
  template <class T, template <class> class P>
    void
    SubIteratorBase<T,P>::init()
  {
      iter_ = 0;
  }

  template <class T, template <class> class P>
    SubIterator&
    SubIteratorBase<T,P>::operator++()
  {
      ++iter_;
      return *this;
  }

  template <class T, template <class> class P>
    SubIteratorBase<T,P>::operator int() const
  {
    return (iter_ < num_iter()) ? 1 : 0;
  }

  template <class T, template <class> class P>
    void
    SubIteratorBase<T,P>::init_subobj()
  {
      P<T>::init_subobj(obj_,subobj_);
  }

  template <class T, template <class> class P>
    void
    SubIteratorBase<T,P>::delete_subobj()
  {
      P<T>::dealloc_subobj(subobj_);
  }
  
};

#endif
