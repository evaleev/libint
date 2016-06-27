// taken from http://codereview.stackexchange.com/questions/20058/c11-any-class
// and modified slightly.
// this code is in public domain

#ifndef _libint2_include_libint2_util_any_h_
#define _libint2_include_libint2_util_any_h_

#include <type_traits>
#include <utility>
#include <typeinfo>
#include <string>
#include <cassert>

namespace libint2 {

  /// emulates boost::any
  struct any
  {
    public:

      template<class T>
      using StorageType = typename std::decay<T>::type;

      any() = default;
      any(const any& that) : handle_(that.clone()) {}
      any(any&& that) : handle_(std::move(that.handle_)) { }
      template <typename U, typename = typename std::enable_if<
                                not std::is_same<any, U>::value>::type>
      any(U&& value)
          : handle_(new handle<StorageType<U>>(std::forward<U>(value))) {}

      any& operator=(const any& a)
      {
        any tmp(a);
        std::swap(*this, tmp);
        return *this;
      }
      template <typename U, typename = typename std::enable_if<
                                not std::is_same<any, U>::value>::type>
      any& operator=(U a) {
        any tmp(std::forward<U>(a));
        std::swap(*this, tmp);
        return *this;
      }
      any& operator=(any&& a)
      {
        std::swap(handle_, a.handle_);
        return *this;
      }

      operator bool() const { return bool(handle_); }

      template<class U> bool is() const
      {
          typedef StorageType<U> T;
          auto derived = dynamic_cast<handle<T>*> (handle_.get());
          return derived;
      }

      /// \note if NDEBUG is not defined, will throw \c std::bad_cast if U is not the stored type
      template<class U>
      StorageType<U>& as()
      {
          typedef StorageType<U> T;

#if not defined(NDEBUG)
          auto derived = dynamic_cast<handle<T>*> (handle_.get());
          if (!derived)
              throw std::bad_cast();
#else // NDEBUG
          auto derived = static_cast<handle<T>*> (handle_.get());
#endif

          return derived->value;
      }

      /// \note if NDEBUG is not defined, will throw \c std::bad_cast if U is not the stored type
      template<class U>
      const StorageType<U>& as() const
      {
          typedef StorageType<U> T;

#if not defined(NDEBUG)
          auto derived = dynamic_cast<handle<T>*> (handle_.get());
          if (!derived)
              throw std::bad_cast();
#else // NDEBUG
          auto derived = static_cast<handle<T>*> (handle_.get());
#endif

          return derived->value;
      }

      template<class U>
      operator U()
      {
          return as<StorageType<U>>();
      }


    private:
      struct handle_base
      {
          virtual ~handle_base() {}

          virtual handle_base* clone() const = 0;
      };

      template<typename T>
      struct handle : handle_base
      {
          template<typename U> handle(U&& value) : value(std::forward<U>(value)) { }

          T value;

          handle_base* clone() const { return new handle<T>(value); }
      };

      handle_base* clone() const
      {
          if (handle_)
              return handle_->clone();
          else
              return nullptr;
      }

      std::unique_ptr<handle_base> handle_;
  };

} // namespace libint2

#endif // header guard

