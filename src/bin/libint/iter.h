
#ifndef _libint2_src_bin_libint_iter_h_
#define _libint2_src_bin_libint_iter_h_

#include <rr.h>
#include <integral.h>

using namespace std;


namespace libint2 {

  /** SetIterator<T> provides an iterator class for T. It iterates through
      T as if it were a set of some data of type T::iter_type.
      For example, it can be used to access individual elements in a shell quertet.
  */
  template <class T> class SetIterator {

  public:
    SetIterator(const T*);
    ~SetIterator();

    typedef typename T::iter_type iter_type;

    /// Returns a number of iterations (number of elements in a set over which to iterate).
    const unsigned int num_iter() const;
    /// Returns the first element.
    const iter_type* first();
    /// Returns next element.
    const iter_type* next();

  private:
    const T* obj_;
    vector<const iter_type*> data_;

  };

  template <class T>
    SetIterator<T>::SetIterator(const T* obj) :
    obj_(obj), data_(num_elem())
    {
    }
  
  template <class T>
    SetIterator<T>::~SetIterator()
    {
    }

  template <class T>
    const unsigned int
    SetIterator<T>::num_iter() const
    {
      throw std::runtime_error("SetIterator<T>::num_elem() -- no specialization is provided");
    }

  template <class T>
    const typename SetIterator<T>::iter_type*
    SetIterator<T>::first()
    {
      throw std::runtime_error("SetIterator<T>::first() -- no specialization is provided");
    }

  template <class T>
    const typename SetIterator<T>::iter_type*
    SetIterator<T>::next()
    {
      throw std::runtime_error("SetIterator<T>::next() -- no specialization is provided");
    }

  /** SetIterator_traits<I> is a collection of traits for SetIterator I.
      This is what is used to obtain information about I.
      Consider this a laundry list of what I must provide in
      order to be useable.
   */
  template <class SetIter> struct SetIterator_traits {
    typedef typename SetIter::iter_type iter_type;
  };


};

#endif
