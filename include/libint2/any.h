// taken from http://codereview.stackexchange.com/questions/20058/c11-any-class
// and modified slightly.
// this code is in public domain

#ifndef _libint2_include_libint2_any_h_
#define _libint2_include_libint2_any_h_

#include <type_traits>
#include <utility>
#include <typeinfo>
#include <string>
#include <cassert>

namespace libint2 {

  template<class T>
  using StorageType = typename std::decay<T>::type;

  /// emulates boost::any
  struct Any
  {
      bool is_null() const { return !ptr; }
      bool not_null() const { return ptr; }

      template<typename U> Any(U&& value)
          : ptr(new Derived<StorageType<U>>(std::forward<U>(value)))
      {

      }

      template<class U> bool is() const
      {
          typedef StorageType<U> T;

          auto derived = dynamic_cast<Derived<T>*> (ptr);

          return derived;
      }

      /// \note if NDEBUG is not defined, will throw if U is not the stored type
      template<class U>
      StorageType<U>& as()
      {
          typedef StorageType<U> T;

#if not defined(NDEBUG)
          auto derived = dynamic_cast<Derived<T>*> (ptr);
          if (!derived)
              throw std::bad_cast();
#else // NDEBUG
          auto derived = static_cast<Derived<T>*> (ptr);
#endif

          return derived->value;
      }

      /// \note if NDEBUG is not defined, will throw if U is not the stored type
      template<class U>
      const StorageType<U>& as() const
      {
          typedef StorageType<U> T;

#if not defined(NDEBUG)
          auto derived = dynamic_cast<Derived<T>*> (ptr);
          if (!derived)
              throw std::bad_cast();
#else // NDEBUG
          auto derived = static_cast<Derived<T>*> (ptr);
#endif

          return derived->value;
      }

      template<class U>
      operator U()
      {
          return as<StorageType<U>>();
      }

      Any()
          : ptr(nullptr)
      {

      }

      Any(Any& that)
          : ptr(that.clone())
      {

      }

      Any(Any&& that)
          : ptr(that.ptr)
      {
          that.ptr = nullptr;
      }

      Any(const Any& that)
          : ptr(that.clone())
      {

      }

      Any(const Any&& that)
          : ptr(that.clone())
      {

      }

      Any& operator=(const Any& a)
      {
          if (ptr == a.ptr)
              return *this;

          auto old_ptr = ptr;

          ptr = a.clone();

          if (old_ptr)
              delete old_ptr;

          return *this;
      }

      Any& operator=(Any&& a)
      {
          if (ptr == a.ptr)
              return *this;

          std::swap(ptr, a.ptr);

          return *this;
      }

      ~Any()
      {
          if (ptr)
              delete ptr;
      }

  private:
      struct Base
      {
          virtual ~Base() {}

          virtual Base* clone() const = 0;
      };

      template<typename T>
      struct Derived : Base
      {
          template<typename U> Derived(U&& value) : value(std::forward<U>(value)) { }

          T value;

          Base* clone() const { return new Derived<T>(value); }
      };

      Base* clone() const
      {
          if (ptr)
              return ptr->clone();
          else
              return nullptr;
      }

      Base* ptr;
  };

} // namespace libint2

#endif // header guard

