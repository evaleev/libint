
#ifndef _libint2_src_lib_libint_vectorppc_h_
#define _libint2_src_lib_libint_vectorppc_h_

// clang on BG/Q uses this header file for intrinsics
#if defined(__clang__) && defined(__bgq__)
# include <qpxintrin.h>
#endif

#ifdef __VECTOR4DOUBLE__

namespace libint2 {

  struct VectorQPXDouble {

      typedef double T;
      vector4double d;

      VectorQPXDouble() {}

      VectorQPXDouble(T a) {
        d = vec_splats(a);
      }

      VectorQPXDouble& operator=(T a) {
        d = vec_splats(a);
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

  inline VectorQPXDouble fma_plus(VectorQPXDouble a, VectorQPXDouble b, VectorQPXDouble c) {
    VectorQPXDouble d;
    d.d = vec_madd(a.d, b.d, c.d);
    return d;
  }
  inline VectorQPXDouble fma_minus(VectorQPXDouble a, VectorQPXDouble b, VectorQPXDouble c) {
    VectorQPXDouble d;
    d.d = vec_msub(a.d, b.d, c.d);
    return d;
  }

  //@}

};

#endif // QPX-only

// only xlC on BG/L and BG/P supports FP2 (Double Hummer)instructions, not sure how to check if they are enabled
#if (defined(__xlC__) || defined(__clang__)) && (defined(__bgp__) || defined(__blrts__))

#if defined(__xlC__)
# include <builtins.h>
#endif
#if defined(__clang__)
# include <fp2intrin.h>
#endif

namespace libint2 {

  struct VectorFP2Double {

      typedef double T;
      double _Complex d; //< represents 2 doubles

      VectorFP2Double() {}

      VectorFP2Double(T a) {
        T a01[2]; a01[0] = a; a01[1] = a;
        d = __lfpd(&a01[0]);
      }

      VectorFP2Double& operator=(T a) {
        T a01[2]; a01[0] = a; a01[1] = a;
        d = __lfpd(&a01[0]);
        return *this;
      }

      VectorFP2Double& operator+=(VectorFP2Double a) {
        d = __fpadd(d, a.d);
        return *this;
      }

      VectorFP2Double& operator-=(VectorFP2Double a) {
        d = __fpsub(d, a.d);
        return *this;
      }

      operator double() const {
        double d0 = __creal(d);
        return d0;
      }

  };

  //@{ arithmetic operators
  inline VectorFP2Double operator*(double a, VectorFP2Double b) {
    VectorFP2Double c;
    c.d = __fxpmul(b.d, a);
    return c;
  }

  inline VectorFP2Double operator*(VectorFP2Double a, double b) {
    VectorFP2Double c;
    c.d = __fxpmul(a.d, b);
    return c;
  }

  inline VectorFP2Double operator*(int a, VectorFP2Double b) {
    if (a == 1)
      return b;
    else {
      VectorFP2Double c;
      c.d = __fxpmul(b.d, (double)a);
      return c;
    }
  }

  inline VectorFP2Double operator*(VectorFP2Double a, int b) {
    if (b == 1)
      return a;
    else {
      VectorFP2Double c;
      c.d = __fxpmul(a.d, (double)b);
      return c;
    }
  }

  inline VectorFP2Double operator*(VectorFP2Double a, VectorFP2Double b) {
    VectorFP2Double c;
    c.d = __fpmul(a.d, b.d);
    return c;
  }

  inline VectorFP2Double operator+(VectorFP2Double a, VectorFP2Double b) {
    VectorFP2Double c;
    c.d = __fpadd(a.d, b.d);
    return c;
  }

  inline VectorFP2Double operator-(VectorFP2Double a, VectorFP2Double b) {
    VectorFP2Double c;
    c.d = __fpsub(a.d, b.d);
    return c;
  }

  /* there's no division DH instruction that I can see
  inline VectorFP2Double operator/(VectorFP2Double a, VectorFP2Double b) {
    VectorFP2Double c;
  }
  */

  inline VectorFP2Double fma_plus(VectorFP2Double a, VectorFP2Double b, VectorFP2Double c) {
    VectorFP2Double d;
    d.d = __fpmadd(a.d, b.d, c.d);
    return d;
  }
  inline VectorFP2Double fma_minus(VectorFP2Double a, VectorFP2Double b, VectorFP2Double c) {
    VectorFP2Double d;
    d.d = __fpmsub(a.d, b.d, c.d);
    return d;
  }

  //@}

};

#endif // FP2-only

#endif // header guard

