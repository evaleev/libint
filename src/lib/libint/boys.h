// prototype for the Boys function engines (Boys function = Fm(T))
// the original Chebyshev extrapolation code is from ORCA, due to Frank Neese

#include <cmath>
#include <vector.h>
typedef libint2::simd::VectorSSEDouble REALTYPE; // for now REALTYPE will be SSE2 type, eventually this will be defined elsewhere and the Interpolate will become
                                                 // a template (or likely a macro since OpenCL does not support templates as of spec 1.2)


/* --------------------------------------------------------------------
   Define factorials
   -------------------------------------------------------------------- */
#define FAK0        1.0
#define FAK1        1.0
#define FAK2        2.0
#define FAK3        6.0
#define FAK4       24.0
#define FAK5      120.0
#define FAK6      720.0
#define FAK7     5040.0
#define FAK8    40320.0
#define FAK9   362880.0
#define FAK10 3628800.0
#define FAK11 11.0*FAK10
#define FAK12 12.0*FAK11
#define FAK13 13.0*FAK12
#define FAK14 14.0*FAK13
#define FAK15 15.0*FAK14
#define FAK16 16.0*FAK15
#define FAK17 17.0*FAK16
#define FAK18 18.0*FAK17
#define FAK19 19.0*FAK18
#define FAK20 20.0*FAK19
#define FAK21 21.0*FAK20
#define FAK22 22.0*FAK21
#define FAK23 23.0*FAK22
#define FAK24 24.0*FAK23
#define FAK25 25.0*FAK24
#define FAK26 26.0*FAK25
#define FAK27 27.0*FAK26
#define FAK28 28.0*FAK27
#define FAK29 29.0*FAK28
#define FAK30 30.0*FAK29
#define FAK31 31.0*FAK30
#define FAK32 32.0*FAK31

static const double
          Fak[33] =
          {       FAK0,  /*  0 */
                  FAK1,  /*  1 */
                  FAK2,  /*  2 */
                  FAK3,  /*  3 */
                  FAK4,  /*  4 */
                  FAK5,  /*  5 */
                  FAK6,  /*  6 */
                  FAK7,  /*  7 */
                  FAK8,  /*  8 */
                  FAK9,  /*  9 */
                  FAK10, /* 10 */
                  FAK11, /* 11 */
                  FAK12, /* 12 */
                  FAK13, /* 13 */
                  FAK14, /* 14 */
                  FAK15, /* 15 */
                  FAK16, /* 16 */
                  FAK17, /* 17 */
                  FAK18, /* 18 */
                  FAK19, /* 19 */
                  FAK20, /* 20 */
                  FAK21, /* 21 */
                  FAK22, /* 22 */
                  FAK23, /* 23 */
                  FAK24, /* 24 */
                  FAK25, /* 25 */
                  FAK26, /* 26 */
                  FAK27, /* 27 */
                  FAK28, /* 28 */
                  FAK29, /* 29 */
                  FAK30, /* 30 */
                  FAK31, /* 31 */
                  FAK32  /* 32 */
          };

#define _local_min_macro(a,b) ((a) > (b) ? (a) : (b))

/// Slow for the sake of precision -- only use for reference purposes
template <typename Real>
struct FmEval_Reference {

    /// computes a single value of Fm(T) using MacLaurin series.
    static Real eval(Real T, size_t m, Real absolute_precision) {
      Real denom = (m+0.5);
      Real term = 0.5 * exp(-T) / denom;
      Real sum = term;
      Real rel_error;
      Real epsilon;
      const Real relative_zero = 1e-15;
      const Real absolute_precision_o_10 = absolute_precision * 0.1;
      do {
        denom += 1.0;
        term *= T/denom;
        sum += term;
        rel_error = term/sum;
        // stop if adding a term smaller or equal to absolute_precision/10 and smaller than relative_zero * sum
        // When sum is small in absolute value, the second threshold is more important
        epsilon = _local_min_macro(absolute_precision_o_10, sum*relative_zero);
      } while (term > epsilon);

      return sum;
    }

    /// fills up Fm from the top
    static void eval(Real* Fm, Real T, size_t mmax, Real absolute_precision) {

      // evaluate for mmax using MacLaurin series
      // it converges fastest for the largest m -> use it to compute Fmmax(T)
      //  see JPC 94, 5564 (1990).
      Fm[mmax] = eval(T, mmax, absolute_precision);
      /* And then do downward recursion */
      if (mmax > 0) {
        const Real T2 = 2.0 * T;
        const Real exp_T = exp(-T);
        for (int m = mmax - 1; m >= 0; m--)
          Fm[m] = (Fm[m + 1] * T2 + exp_T) / (2*m + 1);
      }
    }

};

/// Computes the Boys function, Fm(T), using Chebyshev interpolation from a precomputed table of values
class FmEval_Chebyshev3 {

    const int FM_N = 2048;
    const int ORDER = 4;
    const double FM_MAX;
    const double FM_DELTA;

    int mmax; /* the maximum m that is tabulated */
    double **c; /* the Chebyshev coefficients for a given m */
    double twoi1[64];/* needed for downward recursion */
    double ihalf[64];/* needed for upward recursion */

  public:
    FmEval_Chebyshev3(int m_max) :
        mmax(m_max), FM_MAX(25.0), FM_DELTA(FM_MAX/(FM_N-1)) {
      for (int i = 0; i < 64; i++) {
        twoi1[i] = 1.0 / (2.0 * i + 1.0);
        // in the time critical step
        ihalf[i] = double(i) - 0.5;
      }

      init();
    }
    ~FmEval_Chebyshev3() {
      delete[] c[0];
      delete c;
    }

// -----------------------------------------------------
// the original Boys function from ORCA
// computes F_m(x) for all m = 0 ... mmax
// F_mmax(x) is evaluated by extrapolation, the rest
// by downward recursion
// -----------------------------------------------------
    inline void eval(double* Fm, double x, int mmax) {

      // ---------------------------------------------
      // large arguments. The total cost is:
      //  1       SQRT
      // +1       FLOP
      // +2*(m-1) FLOPS
      // ---------------------------------------------
      if (x>FM_MAX){
        Fm[0]= 0.88622692545275801365/sqrt(x);
        if (mmax==0) return;
        //for (i=1;i<=m;i++) Fm[i] = Fm[i-1]*(double(i)-0.5)/x;
        for (int i=1;i<=mmax;i++) Fm[i] = Fm[i-1]*ihalf[i]/x;
        return;
      }

      // ---------------------------------------------
      // small and intermediate arguments. The total
      // cost is 6       FLOPS
      //        +3*(m-1) FLOPS
      //        +2       FLOPS
      //        +1       EXP
      // ---------------------------------------------
      const double *d = c[mmax]; // a pointer to the correct m-vector
      // about which point on the grid to interpolate?
      const double xstep = FM_DELTA; // the interpolation interval -- hardwired in ORCA, should be determined by the target precision and the interpolation order
      const double xd = x / xstep;
      const int iv = int(xd); // the interval
      const int ofs = iv * 4; // the offset in the interpolation table
      // for the largest m evaluate by interpolation (6 FLOPS)
      Fm[mmax] = d[ofs] + xd * (d[ofs + 1] + xd * (d[ofs + 2] + xd * d[ofs + 3]));

      // check against the reference value
      {
        double refvalue = FmEval_Reference<double>::eval(x, mmax, 1e-15);
        if (abs(refvalue - Fm[mmax]) > 1e-10) {
          std::cout << "T = " << x << " m = " << mmax << " cheb = " << Fm[mmax] << " ref = " << refvalue << std::endl;
        }
      }

      // use downward recursion to make the other members (3 FLOPS/member)
      if (mmax > 0) {
        const double x2 = 2.0 * x;
        const double exp_x = exp(-x);
        for (int m = mmax - 1; m >= 0; m--)
          Fm[m] = (Fm[m + 1] * x2 + exp_x) * twoi1[m];
      }
    }

// -----------------------------------------------------
// the vectorized Boys function
// computes F_m(x) for all m = 0 ... mmax
// x is a scalar or a SIMD-type vector, hence the output is an array of scalars/SIMD vectors
// F_mmax(x) is evaluated by extrapolation, the rest
// by downward recursion
// -----------------------------------------------------
    inline void eval(REALTYPE *Fm, REALTYPE x, int mmax) {

      abort(); // the rest is to be implemented

#if 0
      // ---------------------------------------------
      // large arguments. The total cost is:
      //  1       SQRT
      // +1       FLOP
      // +2*(m-1) FLOPS
      // ---------------------------------------------
      if (x>FM_MAX){
        Fm[0]= 0.88622692545275801365/sqrt(x);
        if (mmax==0) return;
        //for (i=1;i<=m;i++) Fm[i] = Fm[i-1]*(double(i)-0.5)/x;
        for (int i=1;i<=mmax;i++) Fm[i] = Fm[i-1]*ihalf[i]/x;
        return;
      }
#endif

      // ---------------------------------------------
      // small and intermediate arguments. The total
      // cost is 6       FLOPS
      //        +3*(m-1) FLOPS
      //        +2       FLOPS
      //        +1       EXP
      // ---------------------------------------------
      const double *d = c[mmax]; // a pointer to the correct m-vector

      REALTYPE xstep(FM_DELTA,
                     FM_DELTA);
      REALTYPE xd = x / xstep; // the interpolation variable

#if 0
      // convert xd into iv and ofs
      // iv and ofs must be the integer scalar/vector types corresponding to REALTYPE, hence need type traits

      //
      // find value of the largest member by interpolation (6 FLOPS)
      //

      // these will be held in vector registers
      // note that these are initialized with values from noncontiguous memory locations:
      // d0.0 = d[ofs.0] and d0.1 = d[ofs.1], hence for any ofs.0 != ofs.1 the loads will not be packed (coalesced in gpu lingo)
      REALTYPE d0,// d[ofs]
      d1,// d[ofs+1]
      d2,// d[ofs+2]
      d3;// d[ofs+3]
      Fm[mmax] = d0 + xd * (d1 + xd * (d2 + xd * d3));

      if (mmax > 0) {
        // use downward recursion to make the other members (3 FLOPS/member)
        REALTYPE x2 = 2.0 * x;
        REALTYPE exp_x = exp(-x);
        for (i = m - 1; i >= 0; i--)
        Fm[i] = (Fm[i + 1] * x2 + exp_x) * twoi1_realtype[i];
      }
#endif
    }

  private:

    /* ----------------------------------------------------------------------------
       This function here creates the expansion coefficients for a single interval

               ON INPUT  a,b  : the interval boundaries
                         cc   : a pointer to the appropriate place in
                                the coefficient table
                         m    : the F[m] to generate
               ON OUTPUT cc   : cc[0]-cc[3] hold the coefficients
       ---------------------------------------------------------------------------- */
    void MakeCoeffs(double a,double b,double *cc,int m)
    { int k,j;
      double f[128],ac[128],Fm[128];
      double sum;
      // characterize the interval
      double TwoDelta = b-a;
      double Delta    = 0.5*TwoDelta;
      double HalfDelta= 0.5*Delta;
      double XXX      = a+Delta;

      const double absolute_precision = 1e-20;  // compute as precisely as possible
      FmEval_Reference<double>::eval(Fm,XXX,m+ORDER+20,absolute_precision);

      for (k=0;k<=ORDER+20;k++){
        if ((k%2)== 0) f[k]= Fm[k+m];
        else           f[k]=-Fm[k+m];
      }
      // calculate the coefficients a
      double fac;
      for (j=0;j<ORDER;j++){
        if (j==0) fac=1.0;
        else      fac=2.0*pow(HalfDelta,(double)j);
        sum  =0.0;
        for (k=0;k<10;k++)
          sum += f[j+2*k]*pow(HalfDelta,(double)(2*k))/Fak[k]/Fak[k+j];
        ac[j] = fac*sum;
      }
      // calculate the coefficients c that are Gill's f's
      double arg =-XXX/Delta;
      double arg2= arg*arg;
      double arg3= arg2*arg;
      cc[0]= (ac[0]-ac[2])
            +(ac[1]-3.0*ac[3])*arg
            +2.0*ac[2]*arg2+4.0*ac[3]*arg3;
      cc[1]= (2.0*ac[1]-6.0*ac[3])+8.0*ac[2]*arg+24.0*ac[3]*arg2;
      cc[2]= 8.0*ac[2]+48.0*ac[3]*arg;
      cc[3]= 32.0*ac[3];
    };

    /* ----------------------------------------------------------------------------
       This function here make the expansion coefficients for all intervals


               ON INPUT  m    : the highest F[m] to generate

               ON OUTPUT  c   : the coefficients c[m][i] are generated
       ---------------------------------------------------------------------------- */
     void init() {
       int iv,im;

       // get memory
       c= new double*[mmax+1];
       c[0] = new double[(mmax+1) * FM_N*ORDER];

       // loop over all m values and make the coefficients
       for (im=0;im<=mmax;im++){
         if (im > 0) c[im]= c[im - 1] + FM_N*ORDER;

         // make expansion coefficients for this particular m
         for (iv=0;iv<FM_N;iv++){
           const double a   =iv*FM_DELTA;
           const double b   =a+FM_DELTA;
           MakeCoeffs(a,b,&(c[im][iv*ORDER]),im);
         }
       }
    }

};

