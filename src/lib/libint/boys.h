// prototype for the Boys function engines (Boys function = Fm(T))
// the original Chebyshev extrapolation code is from ORCA, due to Frank Neese

#ifndef _libint2_src_lib_libint_boys_h_
#define _libint2_src_lib_libint_boys_h_

#if defined(__cplusplus)

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector.h>
#include <cassert>
#include <vector>

namespace libint2 {

  /// holds tables of expensive quantities
  template<typename Real>
  class ExpensiveNumbers {
    public:
      ExpensiveNumbers(int ifac = -1, int idf = -1, int ibc = -1) {
        if (ifac >= 0) {
          fac.resize(ifac + 1);
          fac[0] = 1.0;
          for (int i = 1; i <= ifac; i++) {
            fac[i] = i * fac[i - 1];
          }
        }

        if (idf >= 0) {
          df.resize(idf + 1);
          /* df[i] gives (i-1)!!, so that (-1)!! is defined... */
          df[0] = 1.0;
          if (idf >= 1)
            df[1] = 1.0;
          if (idf >= 2)
            df[2] = 1.0;
          for (int i = 3; i <= idf; i++) {
            df[i] = (i - 1) * df[i - 2];
          }
        }

        if (ibc >= 0) {
          bc_.resize((ibc+1)*(ibc+1));
          std::fill(bc_.begin(), bc_.end(), Real(0));
          bc.resize(ibc+1);
          bc[0] = &bc_[0];
          for(int i=1; i<=ibc; ++i)
            bc[i] = bc[i-1] + (ibc+1);

          for(int i=0; i<=ibc; i++)
            bc[i][0] = 1.0;
          for(int i=0; i<=ibc; i++)
            for(int j=1; j<=i; ++j)
              bc[i][j] = bc[i][j-1] * Real(i-j+1) / Real(j);
        }

        for (int i = 0; i < 64; i++) {
          twoi1[i] = 1.0 / (Real(2.0) * i + Real(1.0));
          ihalf[i] = Real(i) - Real(0.5);
        }

      }

      ~ExpensiveNumbers() {
      }

      std::vector<Real> fac;
      std::vector<Real> df;
      std::vector<Real*> bc;

      // these quantitites are needed with indices <= mmax
      // 64 is sufficient to handle up to 4 center integrals with up to L=15 basis functions
      Real twoi1[64]; /* 1/(2 i + 1); needed for downward recursion */
      Real ihalf[64]; /* i - 0.5, needed for upward recursion */

    private:
      std::vector<Real> bc_;
  };

#define _local_min_macro(a,b) ((a) > (b) ? (a) : (b))

  /// Slow for the sake of precision -- only use for reference purposes
  template<typename Real>
  struct FmEval_Reference {

      /// computes a single value of Fm(T) using MacLaurin series.
      static Real eval(Real T, size_t m, Real absolute_precision) {
        Real denom = (m + 0.5);
        Real term = 0.5 * exp(-T) / denom;
        Real sum = term;
        Real rel_error;
        Real epsilon;
        const Real relative_zero = 1e-15;
        const Real absolute_precision_o_10 = absolute_precision * 0.1;
        do {
          denom += 1.0;
          term *= T / denom;
          sum += term;
          rel_error = term / sum;
          // stop if adding a term smaller or equal to absolute_precision/10 and smaller than relative_zero * sum
          // When sum is small in absolute value, the second threshold is more important
          epsilon = _local_min_macro(absolute_precision_o_10, sum*relative_zero);
        } while (term > epsilon);

        return sum;
      }

      /// fills up Fm from the top using downward recursion
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
            Fm[m] = (Fm[m + 1] * T2 + exp_T) / (2 * m + 1);
        }
      }

  };

  /** Computes the Boys function, Fm(T), using Chebyshev interpolation from a precomputed table of values
    * based on the code from ORCA by Dr. Frank Neese.
    */
  class FmEval_Chebyshev3 {

      static const int FM_N = 2048;
      static const int ORDER = 4;
      const double FM_MAX;
      const double FM_DELTA;
      const double FM_one_over_DELTA;

      int mmax; /* the maximum m that is tabulated */
      double **c; /* the Chebyshev coefficients for a given m */
      ExpensiveNumbers<double> numbers_;

    public:
      FmEval_Chebyshev3(int m_max) :
          FM_MAX(30.0), // this translates in appr. 1e-15  error in upward recursion, see the note below
          FM_DELTA(FM_MAX / (FM_N - 1)),
          FM_one_over_DELTA(1.0 / FM_DELTA),
          mmax(m_max), numbers_(14, 0) {
        assert(mmax <= 63);
        init();
      }
      ~FmEval_Chebyshev3() {
        delete[] c[0];
        delete c;
      }

      inline void eval(double* Fm, double x, int mmax) const {

        // large T => use upward recursion
        // cost = 1 div + 1 sqrt + (1 + 2*(m-1)) muls
        if (x > FM_MAX) {
          const double one_over_x = 1.0/x;
          Fm[0] = 0.88622692545275801365 * sqrt(one_over_x); // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
          if (mmax == 0)
            return;
          // this upward recursion formula omits - e^(-x)/(2x), which for x>FM_MAX is <1e-15
          for (int i = 1; i <= mmax; i++)
            Fm[i] = Fm[i - 1] * numbers_.ihalf[i] * one_over_x; // see Eq. (9.8.13)
          return;
        }

        // ---------------------------------------------
        // small and intermediate arguments => interpolate Fm and downward recursion
        // cost is 6       FLOPS
        //        +3*(m-1) FLOPS
        //        +2       FLOPS
        //        +1       EXP
        // ---------------------------------------------
        const double *d = c[mmax]; // a pointer to the correct m-vector
        // about which point on the grid to interpolate?
        const double xd = x * FM_one_over_DELTA;
        const int iv = int(xd); // the interval
        const int ofs = iv * 4; // the offset in the interpolation table
        // for the largest m evaluate by interpolation (6 FLOPS)
        Fm[mmax] = d[ofs]
            + xd * (d[ofs + 1] + xd * (d[ofs + 2] + xd * d[ofs + 3]));

//        // check against the reference value
//        if (false) {
//          double refvalue = FmEval_Reference<double>::eval(x, mmax, 1e-15); // compute F(T) with m=mmax
//          if (abs(refvalue - Fm[mmax]) > 1e-10) {
//            std::cout << "T = " << x << " m = " << mmax << " cheb = "
//                << Fm[mmax] << " ref = " << refvalue << std::endl;
//          }
//        }

        // use downward recursion (Eq. (9.8.14))
        // cost = 1 exp + (1 + 3*m mul/add)
        if (mmax > 0) {
          const double x2 = 2.0 * x;
          const double exp_x = exp(-x);
          for (int m = mmax - 1; m >= 0; m--)
            Fm[m] = (Fm[m + 1] * x2 + exp_x) * numbers_.twoi1[m];
        }
      }

// -----------------------------------------------------
// the vectorized Boys function
// computes F_m(x) for all m = 0 ... mmax
// x is a scalar or a SIMD-type vector, hence the output is an array of scalars/SIMD vectors
// F_mmax(x) is evaluated by extrapolation, the rest
// by downward recursion
// -----------------------------------------------------
#define __COMMENT_OUT_VECTORIZED_EVAL 1
#if not(__COMMENT_OUT_VECTORIZED_EVAL)
      typedef libint2::simd::VectorSSEDouble REALTYPE; // for now REALTYPE will be SSE2 type, eventually this will be defined elsewhere and the Interpolate will become
                                                       // a template (or likely a macro since OpenCL does not support templates as of spec 1.2)
      inline void eval(REALTYPE *Fm, REALTYPE x, int mmax) const {

        abort(); // the rest is to be implemented

#if 0
        // ---------------------------------------------
        // large arguments. The total cost is:
        //  1       SQRT
        // +1       FLOP
        // +2*(m-1) FLOPS
        // ---------------------------------------------
        if (x>FM_MAX) {
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

        REALTYPE xstep(FM_DELTA, FM_DELTA);
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
#endif

    private:

      /* ----------------------------------------------------------------------------
       This function here creates the expansion coefficients for a single interval

       ON INPUT  a,b  : the interval boundaries
       cc   : a pointer to the appropriate place in
       the coefficient table
       m    : the F[m] to generate
       ON OUTPUT cc   : cc[0]-cc[3] hold the coefficients
       ---------------------------------------------------------------------------- */
      void MakeCoeffs(double a, double b, double *cc, int m) {
        int k, j;
        double f[128], ac[128], Fm[128];
        double sum;
        // characterize the interval
        double TwoDelta = b - a;
        double Delta = 0.5 * TwoDelta;
        double HalfDelta = 0.5 * Delta;
        double XXX = a + Delta;

        const double absolute_precision = 1e-20; // compute as precisely as possible
        FmEval_Reference<double>::eval(Fm, XXX, m + ORDER + 20,
                                       absolute_precision);

        for (k = 0; k <= ORDER + 20; k++) {
          if ((k % 2) == 0)
            f[k] = Fm[k + m];
          else
            f[k] = -Fm[k + m];
        }
        // calculate the coefficients a
        double fac;
        for (j = 0; j < ORDER; j++) {
          if (j == 0)
            fac = 1.0;
          else
            fac = 2.0 * pow(HalfDelta, (double) j);
          sum = 0.0;
          for (k = 0; k < 10; k++)
            sum += f[j + 2 * k] * pow(HalfDelta, (double) (2 * k)) / numbers_.fac[k]
                / numbers_.fac[k + j];
          ac[j] = fac * sum;
        }
        // calculate the coefficients c that are Gill's f's
        double arg = -XXX / Delta;
        double arg2 = arg * arg;
        double arg3 = arg2 * arg;
        cc[0] = (ac[0] - ac[2]) + (ac[1] - 3.0 * ac[3]) * arg
            + 2.0 * ac[2] * arg2 + 4.0 * ac[3] * arg3;
        cc[1] = (2.0 * ac[1] - 6.0 * ac[3]) + 8.0 * ac[2] * arg
            + 24.0 * ac[3] * arg2;
        cc[2] = 8.0 * ac[2] + 48.0 * ac[3] * arg;
        cc[3] = 32.0 * ac[3];
      }

      /* ----------------------------------------------------------------------------
       This function makes the expansion coefficients for all intervals


       ON INPUT  m    : the highest F[m] to generate

       ON OUTPUT  c   : the coefficients c[m][i] are generated
       ---------------------------------------------------------------------------- */
      void init() {
        int iv, im;

        // get memory
        c = new double*[mmax + 1];
        c[0] = new double[(mmax + 1) * FM_N * ORDER];

        // loop over all m values and make the coefficients
        for (im = 0; im <= mmax; im++) {
          if (im > 0)
            c[im] = c[im - 1] + FM_N * ORDER;

          // make expansion coefficients for this particular m
          for (iv = 0; iv < FM_N; iv++) {
            const double a = iv * FM_DELTA;
            const double b = a + FM_DELTA;
            MakeCoeffs(a, b, &(c[im][iv * ORDER]), im);
          }
        }
      }

  };

#ifndef STATIC_OON
#define STATIC_OON
  namespace {
    const double oon[] = {0.0, 1.0, 1.0/2.0, 1.0/3.0, 1.0/4.0, 1.0/5.0, 1.0/6.0, 1.0/7.0, 1.0/8.0, 1.0/9.0, 1.0/10.0, 1.0/11.0};
  }
#endif

  /// Uses Taylor interpolation of up to 8-th order to compute the Boys function
  template<typename Real, int INTERPOLATION_ORDER = 6>
  class FmEval_Taylor {
    public:
      static const int max_interp_order = 8;
      static const int INTERPOLATION_AND_RECURSION = 1; // compute F_lmax(T) and then iterate down to F_0(T)? Else use interpolation only
      const Real relative_zero_;
      const Real soft_zero_;

      FmEval_Taylor(unsigned int mmax, Real precision) :
          relative_zero_(1e-15), soft_zero_(1e-6), cutoff_(precision), numbers_(
              INTERPOLATION_ORDER + 1, 2 * (mmax + INTERPOLATION_ORDER - 1)) {

        assert(mmax <= 63);

        const Real sqrt_pi = std::sqrt(M_PI);

        /*---------------------------------------
         We are doing Taylor interpolation with
         n=TAYLOR_ORDER terms here:
         error <= delT^n/(n+1)!
         ---------------------------------------*/
        delT_ = 2.0
            * std::pow(cutoff_ * numbers_.fac[INTERPOLATION_ORDER + 1],
                       1.0 / INTERPOLATION_ORDER);
        oodelT_ = 1.0 / delT_;
        max_m_ = mmax + INTERPOLATION_ORDER - 1;

        T_crit_ = new Real[max_m_ + 1]; /*--- m=0 is included! ---*/
        max_T_ = 0;
        /*--- Figure out T_crit for each m and put into the T_crit ---*/
        for (int m = max_m_; m >= 0; --m) {
          /*------------------------------------------
           Damped Newton-Raphson method to solve
           T^{m-0.5}*exp(-T) = epsilon*Gamma(m+0.5)
           The solution is the max T for which to do
           the interpolation
           ------------------------------------------*/
          Real T = -log(cutoff_);
          const Real egamma = cutoff_ * sqrt_pi * numbers_.df[2 * m]
              / std::pow(2.0, m);
          Real T_new = T;
          Real func;
          do {
            const Real damping_factor = 0.2;
            T = T_new;
            /* f(T) = the difference between LHS and RHS of the equation above */
            func = std::pow(T, m - 0.5) * std::exp(-T) - egamma;
            const Real dfuncdT = ((m - 0.5) * std::pow(T, m - 1.5)
                - std::pow(T, m - 0.5)) * std::exp(-T);
            /* f(T) has 2 roots and has a maximum in between. If f'(T) > 0 we are to the left of the hump. Make a big step to the right. */
            if (dfuncdT > 0.0) {
              T_new *= 2.0;
            } else {
              /* damp the step */
              Real deltaT = -func / dfuncdT;
              const Real sign_deltaT = (deltaT > 0.0) ? 1.0 : -1.0;
              const Real max_deltaT = damping_factor * T;
              if (std::fabs(deltaT) > max_deltaT)
                deltaT = sign_deltaT * max_deltaT;
              T_new = T + deltaT;
            }
            if (T_new <= 0.0) {
              T_new = T / 2.0;
            }
          } while (std::fabs(func / egamma) >= soft_zero_);
          T_crit_[m] = T_new;
          const int T_idx = (int) std::floor(T_new / delT_);
          max_T_ = std::max(max_T_, T_idx);
        }

        // allocate the grid (see the comments below)
        {
          const int nrow = max_T_ + 1;
          const int ncol = max_m_ + 1;
          grid_ = new Real*[nrow];
          grid_[0] = new Real[nrow * ncol];
          for (int r = 1; r < nrow; ++r)
            grid_[r] = grid_[r - 1] + ncol;
        }

        /*-------------------------------------------------------
         Tabulate the gamma function from t=delT to T_crit[m]:
         1) include T=0 though the table is empty for T=0 since
         Fm(0) is simple to compute
         -------------------------------------------------------*/
        /*--- do the mmax first ---*/
        const double cutoff_o_10 = 0.1 * cutoff_;
        for (int T_idx = max_T_; T_idx >= 0; --T_idx) {
          const double T = T_idx * delT_;
          libint2::FmEval_Reference<Real>::eval(grid_[T_idx], T, max_m_, cutoff_);
        }
      }

      ~FmEval_Taylor() {
        delete[] T_crit_;
        delete[] grid_[0];
        delete[] grid_;
      }

      void eval(Real* Fm, Real T, int mmax) const {
        const double sqrt_pio2 = 1.2533141373155002512;
        const double two_T = 2.0 * T;

        // stop recursion at mmin
        const int mmin = INTERPOLATION_AND_RECURSION ? mmax : 0;
        /*-------------------------------------
         Compute Fm(T) from mmax down to mmin
         -------------------------------------*/
        const bool use_upward_recursion = true;
        if (use_upward_recursion) {
//          if (T > 30.0) {
          if (T > T_crit_[0]) {
            const double one_over_x = 1.0/T;
            Fm[0] = 0.88622692545275801365 * sqrt(one_over_x); // see Eq. (9.8.9) in Helgaker-Jorgensen-Olsen
            if (mmax == 0)
              return;
            // this upward recursion formula omits - e^(-x)/(2x), which for x>FM_MAX is <1e-15
            for (int i = 1; i <= mmax; i++)
              Fm[i] = Fm[i - 1] * numbers_.ihalf[i] * one_over_x; // see Eq. (9.8.13)
            return;
          }
        }

        // since Tcrit grows with mmax, this condition only needs to be determined once
        if (T > T_crit_[mmax]) {
          double pow_two_T_to_minusjp05 = std::pow(two_T, -mmax - 0.5);
          for (int m = mmax; m >= mmin; --m) {
            /*--- Asymptotic formula ---*/
            Fm[m] = numbers_.df[2 * m] * sqrt_pio2 * pow_two_T_to_minusjp05;
            pow_two_T_to_minusjp05 *= two_T;
          }
        }
        else
        {
          const int T_ind = (int) (0.5 + T * oodelT_);
          const double h = T_ind * delT_ - T;
          const double* F_row = grid_[T_ind] + mmax;

          for (int m = mmax; m >= mmin; --m, --F_row) {

            /*--- Taylor interpolation ---*/
            Real total = 0.0;
            for(int i=INTERPOLATION_ORDER; i>=1; --i)
              total = oon[i]*h*(F_row[i] + total);

            Fm[m] = F_row[0] + total;

          } // interpolation for F_m(T), mmin<=m<=mmax

          // check against the reference value
//          if (false) {
//            double refvalue = FmEval_Reference<double>::eval(T, mmax, 1e-15); // compute F(T) with m=mmax
//            if (abs(refvalue - Fm[mmax]) > 1e-14) {
//              std::cout << "T = " << T << " m = " << mmax << " cheb = "
//                  << Fm[mmax] << " ref = " << refvalue << std::endl;
//            }
//          }

        } // if T < T_crit

        /*------------------------------------
         And then do downward recursion in j
         ------------------------------------*/
        if (INTERPOLATION_AND_RECURSION) {
          if (mmin > 0) {
            double F_mp1 = Fm[mmin];
            const double exp_mT = std::exp(-T);
            for (int m = mmin - 1; m >= 0; --m) {
              const double F_m = (exp_mT + two_T * F_mp1) * numbers_.twoi1[m];
              Fm[m] = F_m;
              F_mp1 = F_m;
            }
          }
        }
      }

    private:
      Real **grid_; /* Table of "exact" Fm(T) values. Row index corresponds to
       values of T (max_T+1 rows), column index to values
       of m (max_m+1 columns) */
      Real delT_; /* The step size for T, depends on cutoff */
      Real oodelT_; /* 1.0 / delT_, see above */
      Real cutoff_; /* Tolerance cutoff used in all computations of Fm(T) */
      int max_m_; /* Maximum value of m in the table, depends on cutoff
       and the number of terms in Taylor interpolation */
      int max_T_; /* Maximum index of T in the table, depends on cutoff
       and m */
      Real *T_crit_; /* Maximum T for each row, depends on cutoff;
       for a given m and T_idx <= max_T_idx[m] use Taylor interpolation,
       for a given m and T_idx > max_T_idx[m] use the asymptotic formula */

      ExpensiveNumbers<Real> numbers_;

      /**
       * Power series estimate of the error introduced by replacing
       * \f$ F_m(T) = \int_0^1 \exp(-T t^2) t^{2 m} \, \mathrm{d} t \f$ with analytically
       * integrable \f$ \int_0^\infinity \exp(-T t^2) t^{2 m} \, \mathrm{d} t = \frac{(2m-1)!!}{2^{m+1}} \sqrt{\frac{\pi}{T^{2m+1}}} \f$
       * @param m
       * @param T
       * @return the error estimate
       */
      static double truncation_error(unsigned int m, double T) {
        const double m2= m * m;
        const double m3= m2 * m;
        const double m4= m2 * m2;
        const double T2= T * T;
        const double T3= T2 * T;
        const double T4= T2 * T2;
        const double T5= T2 * T3;

        const double result = exp(-T) * (105 + 16*m4 + 16*m3*(T - 8) - 30*T + 12*T2
            - 8*T3 + 16*T4 + 8*m2*(43 - 9*T + 2*T2) +
            4*m*(-88 + 23*T - 8*T2 + 4*T3))/(32*T5);
        return result;
      }
      /**
       * Leading-order estimate of the error introduced by replacing
       * \f$ F_m(T) = \int_0^1 \exp(-T t^2) t^{2 m} \, \mathrm{d} t \f$ with analytically
       * integrable \f$ \int_0^\infinity \exp(-T t^2) t^{2 m} \, \mathrm{d} t = \frac{(2m-1)!!}{2^{m+1}} \sqrt{\frac{\pi}{T^{2m+1}}} \f$
       * @param m
       * @param T
       * @return the error estimate
       */
      static double truncation_error(double T) {
        const double result = exp(-T) /(2*T);
        return result;
      }
  };


  //////////////////////////////////////////////////////////
  /// core integral for Yukawa and exponential interactions
  //////////////////////////////////////////////////////////

  /**
   * Evaluates core integral for the Yukawa potential \f$ \exp(- \zeta r) / r \f$
   *
   * P.S. Slow for the sake of precision -- only use for reference purposes
   */
  template<typename Real>
  struct YukawaGmEval_Reference {

      /// computes a single value of G(U,T) using MacLaurin series.
      static Real eval(Real T, size_t m, Real absolute_precision) {
        Real denom = (m + 0.5);
        Real term = 0.5 * exp(-T) / denom;
        Real sum = term;
        Real rel_error;
        Real epsilon;
        const Real relative_zero = 1e-15;
        const Real absolute_precision_o_10 = absolute_precision * 0.1;
        do {
          denom += 1.0;
          term *= T / denom;
          sum += term;
          rel_error = term / sum;
          // stop if adding a term smaller or equal to absolute_precision/10 and smaller than relative_zero * sum
          // When sum is small in absolute value, the second threshold is more important
          epsilon = _local_min_macro(absolute_precision_o_10, sum*relative_zero);
        } while (term > epsilon);

        return sum;
      }

      /// fills up Fm from the top using downward recursion
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
            Fm[m] = (Fm[m + 1] * T2 + exp_T) / (2 * m + 1);
        }
      }

  };

  //////////////////////////////////////////////////////////
  /// core integrals r12^k \sum_i \exp(- a_i r_12^2)
  //////////////////////////////////////////////////////////

  /**
   * Evaluates core integral Gm(rho, T) over a general contracted
   * Gaussian geminal \f$ r_{12}^k \sum_i c_i \exp(- a_i r_{12}^2), \quad k = -1, 0, 2 \f$ .
   * The integrals are needed in R12/F12 methods with STG-nG correlation factors.
   * Specifically, for a correlation factor \f$ f(r_{12}) = \sum_i c_i \exp(- a_i r_{12}^2) \f$
   * integrals with the following kernels are needed:
   * <ul>
   *   <li> \f$ f(r_{12}) \f$  (k=0) </li>
   *   <li> \f$ f(r_{12}) / r_{12} \f$  (k=-1) </li>
   *   <li> \f$ f(r_{12})^2 \f$ (k=0, @sa GaussianGmEval::eval ) </li>
   *   <li> \f$ [f(r_{12}), [\hat{T}_1, f(r_{12})]] \f$ (k=2, @sa GaussianGmEval::eval ) </li>
   * </ul>
   *
   * N.B. ``Asymmetric'' kernels, \f$ f(r_{12}) g(r_{12}) \f$ and
   *   \f$ [f(r_{12}), [\hat{T}_1, g(r_{12})]] \f$, where f and g are two different geminals,
   *   can also be handled straightforwardly.
   *
   */
  template<typename Real, int k>
  struct GaussianGmEval {

      /**
       * @param[in] mmax the evaluator will be used to compute Gm(T) for m <= mmax
       */
      GaussianGmEval(int mmax, Real precision) : mmax_(mmax),
          precision_(precision), fm_eval_(0), numbers_(-1,-1,mmax),
          g_i(mmax+1), r_i(mmax+1), oorhog_i(mmax+1) {
        assert(k == -1 || k == 0 || k == 2);
        // for k=-1 need to evaluate the Boys function
        if (k == -1) {
//          fm_eval_ = new FmEval_Taylor<Real>(mmax_, precision_);
          fm_eval_ = new FmEval_Chebyshev3(mmax_);
          Fm_.resize(mmax_ + 1);
        }
      }

      ~GaussianGmEval() {
        delete fm_eval_;
        fm_eval_ = 0;
      }

      /** fills up Gm(T) from the top using downward recursion.
       *
       * @param[out] Gm
       * @param[in] rho
       * @param[in] T
       * @param[in] mmax
       * @param[in] geminal
       */
      void eval(Real* Gm, Real rho, Real T, size_t mmax,
                const std::vector<std::pair<Real, Real> >& geminal) {

        std::fill(Gm, Gm+mmax+1, Real(0));

        const double sqrt_rho = sqrt(rho);

        typedef typename std::vector<std::pair<Real, Real> >::const_iterator citer;
        const citer gend = geminal.end();
        for(citer i=geminal.begin(); i!= gend; ++i) {

          const double gamma = i->first;
          const double gcoef = i->second;
          const double rhog = rho + gamma;
          const double oorhog = 1.0/rhog;

          const double gorg = gamma * oorhog;
          const double rorg = rho * oorhog;
          const double sqrt_rho_org = sqrt_rho * oorhog;
          const double sqrt_rhog = sqrt(rhog);
          const double sqrt_rorg = sqrt_rho_org * sqrt_rhog;

          /// (ss|g12|ss)
          const double SS_K0G12_SS = rorg * sqrt_rorg * exp(-gorg*T);

          if (k == -1) {
            const double rorgT = rorg * T;
            fm_eval_->eval(&Fm_[0], rorgT, mmax);

#if 0
            const Real const_2_SQRTPI(1.12837916709551257389615890312154517);   /* 2/sqrt(pi)     */
            double pfac = const_2_SQRTPI * sqrt_rhog * SS_K0G12_SS;
            g_i[0] = 1.0;
            r_i[0] = 1.0;
            oorhog_i[0] = 1.0;
            for(int i=1; i<=mmax; i++) {
              g_i[i] = g_i[i-1] * gamma;
              r_i[i] = r_i[i-1] * rho;
              oorhog_i[i] = oorhog_i[i-1] * oorhog;
            }
            for(int m=0; m<=mmax; m++) {
              double ssss = 0.0;
              for(int k=0; k<=m; k++) {
                ssss += numbers_.bc[m][k] * r_i[k] * g_i[m-k] * Fm_[k];
              }
              Gm[m] += gcoef * pfac * ssss * oorhog_i[m];
            }
#endif
            return;
          }

          if (k == 0) {

            double ss_oper_ss_m = SS_K0G12_SS * gcoef;
            Gm[0] += ss_oper_ss_m;
            for(int m=1; m<=mmax; ++m) {
              ss_oper_ss_m *= gorg;
              Gm[m] += ss_oper_ss_m;
            }
            return;
          }

          if (k == 2) {

            /// (ss|g12*r12^2|ss)
            const double rorgT = rorg * T;
            const double SS_K2G12_SS_0 = (1.5 + rorgT) * (SS_K0G12_SS * oorhog);
            const double SS_K2G12_SS_m1 = rorg * (SS_K0G12_SS * oorhog);

            double SS_K2G12_SS_gorg_m = SS_K2G12_SS_0 * gcoef;
            double SS_K2G12_SS_gorg_m1 = SS_K2G12_SS_m1 * gcoef;
            Gm[0] += SS_K2G12_SS_gorg_m;
            for(int m=1; m<=mmax; ++m) {
              SS_K2G12_SS_gorg_m *= gorg;
              Gm[m] += SS_K2G12_SS_gorg_m - m * SS_K2G12_SS_gorg_m1;
              SS_K2G12_SS_gorg_m1 *= gorg;
            }
            return;
          }

        }

      }

    private:
      int mmax_;
      Real precision_; //< absolute precision
//      FmEval_Taylor<Real>* fm_eval_;
      FmEval_Chebyshev3* fm_eval_;
      std::vector<Real> Fm_;

      std::vector<double> g_i;
      std::vector<double> r_i;
      std::vector<double> oorhog_i;


      ExpensiveNumbers<Real> numbers_;
  };


} // end of namespace libint2

#endif // C++ only
#endif // header guard
