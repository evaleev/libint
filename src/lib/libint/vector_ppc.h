
#ifndef _libint2_src_lib_libint_vectorppc_h_
#define _libint2_src_lib_libint_vectorppc_h_

#ifdef __VECTOR4DOUBLE__

namespace libint2 {

  struct VectorQPXDouble {

      typedef double T;
      vector4double d;

      VectorQPXDouble() {}

      VectorQPXDouble(T a) {
        d = vec_splats(a, a, a, a);
      }

      VectorQPXDouble& operator=(T a) {
        d = vec_splats(a, a, a, a);
        return *this;
      }

      VectorQPXDouble& operator+=(VectorQPXDouble a) {
        d = vec_add(d, a.d);
        return *this;
      }

      VectorQPXDouble& operator-=(VectorQPXDouble a) {
        d = vec_sub(d, a.d);
        return *this;
      }

      operator double() const {
        double d0 = vec_extract(d, 0);
        return d0;
      }

  };

  //@{ arithmetic operators
  inline VectorQPXDouble operator*(double a, VectorQPXDouble b) {
    VectorQPXDouble c;
    VectorQPXDouble _a(a);
    c.d = vec_mul(_a.d, b.d);
    return c;
  }

  inline VectorQPXDouble operator*(VectorQPXDouble a, double b) {
    VectorQPXDouble c;
    VectorQPXDouble _b(b);
    c.d = vec_mul(a.d, _b.d);
    return c;
  }

  inline VectorQPXDouble operator*(int a, VectorQPXDouble b) {
    if (a == 1)
      return b;
    else {
      VectorQPXDouble c;
      VectorQPXDouble _a((double)a);
      c.d = vec_mul(_a.d, b.d);
      return c;
    }
  }

  inline VectorQPXDouble operator*(VectorQPXDouble a, int b) {
    if (b == 1)
      return a;
    else {
      VectorQPXDouble c;
      VectorQPXDouble _b((double)b);
      c.d = vec_mul(a.d, _b.d);
      return c;
    }
  }

  inline VectorQPXDouble operator*(VectorQPXDouble a, VectorQPXDouble b) {
    VectorQPXDouble c;
    c.d = vec_mul(a.d, b.d);
    return c;
  }

  inline VectorQPXDouble operator+(VectorQPXDouble a, VectorQPXDouble b) {
    VectorQPXDouble c;
    c.d = vec_add(a.d, b.d);
    return c;
  }

  inline VectorQPXDouble operator-(VectorQPXDouble a, VectorQPXDouble b) {
    VectorQPXDouble c;
    c.d = vec_sub(a.d, b.d);
    return c;
  }

  inline VectorQPXDouble operator/(VectorQPXDouble a, VectorQPXDouble b) {
    VectorQPXDouble c;
    c.d = vec_swdiv(a.d, b.d);
    return c;
  }

  //@}

};

#endif // QPX-only

#endif // header guard

