/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Library General Public License, version 2,
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU Library General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

#ifndef _libint2_src_lib_libint_solidharmonics_h_
#define _libint2_src_lib_libint_solidharmonics_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <array>
#include <vector>

#include <libint2/shell.h>
#include <libint2/cgshell_ordering.h>
#include <Eigen/Dense>

namespace {
  template <typename Int>
  signed char parity(Int i) {
    return i%2 ? -1 : 1;
  }
}

namespace libint2 {

  namespace solidharmonics {

    // to avoid overhead of Eigen::SparseMatrix will roll our own

    /// Transformation coefficients from unnormalized Cartesian Gaussians (rows) to unit-normalized real Solid Harmonics Gaussians.
    /// \note Implemented as a simple fixed-size CSR sparse matrix
    template <typename Real>
    class SolidHarmonicsCoefficients {
      public:
        typedef ::libint2::real_t real_t;

        SolidHarmonicsCoefficients() : l_(-1) {
        }
        SolidHarmonicsCoefficients(unsigned char l) : l_(l) {
          assert(l <= std::numeric_limits<signed char>::max());
          init();
        }
        ~SolidHarmonicsCoefficients() {
          cleanup();
        }
        void init(unsigned char l) {
          assert(l <= std::numeric_limits<signed char>::max());
          l_ = l;
          init();
        }

        /// returns ptr to row values
        const Real* row_values(size_t r) const {
          return values_ + row_offset_[r];
        }
        /// returns ptr to row indices
        const unsigned char* row_idx(size_t r) const {
          return colidx_ + row_offset_[r];
        }
        /// number of nonzero elements in row \c r
        unsigned char nnz(size_t r) const {
          return row_offset_[r+1] - row_offset_[r];
        }

      private:
        Real* values_;  // elements
        unsigned short* row_offset_; // "pointer" to the beginning of each row
        unsigned char* colidx_;  // column indices
        signed char l_;        // the angular momentum quantum number

        void init() {
          const unsigned short npure = 2*l_ + 1;
          const unsigned short ncart = (l_ + 1) * (l_ + 2) / 2;
          std::vector<Real> full_coeff(npure * ncart);

          for(signed char m=-l_; m<=l_; ++m) {
            const signed char pure_idx = m + l_;
            signed char cart_idx = 0;
            signed char lx, ly, lz;
            FOR_CART(lx, ly, lz, l_)
              full_coeff[pure_idx * ncart + cart_idx] = coeff(l_, m, lx, ly, lz);
              ++cart_idx;
            END_FOR_CART
          }

          if (false) {
            std::cout << "SolidHarmonicsCoefficients: l = " << l_ << "\n";
            typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;
            Eigen::Map<Matrix> fcoeff(&full_coeff[0], npure, ncart);
            std::cout << fcoeff << "\n";
          }

          // compress rows
          // 1) count nonzeroes
          size_t nnz = 0;
          for(size_t i=0; i!=full_coeff.size(); ++i)
            nnz += full_coeff[i] == 0.0 ? 0 : 1;
          // 2) allocate
          values_ = new Real[nnz];
          colidx_ = new unsigned char[nnz];
          row_offset_ = new unsigned short[npure+1];
          // 3) copy
          {
            unsigned short pc = 0;
            unsigned short cnt = 0;
            for(unsigned short p=0; p!=npure; ++p) {
              row_offset_[p] = cnt;
              for(unsigned short c=0; c!=ncart; ++c, ++pc) {
                if (full_coeff[pc] != 0.0) {
                  values_[cnt] = full_coeff[pc];
                  colidx_[cnt] = c;
                  ++cnt;
                }
              }
            }
            row_offset_[npure] = cnt;
          }
          // done
        }

        void cleanup() {
          delete[] values_;
          delete[] row_offset_;
          delete[] colidx_;
        }

        /*!---------------------------------------------------------------------------------------------
          Computes coefficient of a cartesian Gaussian in a real solid harmonic Gaussian
          See IJQC 54, 83 (1995), eqn (15). If m is negative, imaginary part is computed, whereas
          a positive m indicates that the real part of spherical harmonic Ylm is requested.
         ---------------------------------------------------------------------------------------------*/
        static double coeff(int l, int m, int lx, int ly, int lz) {
          using libint2::math::fac;
          using libint2::math::df_Kminus1;
          using libint2::math::bc;

          auto abs_m = std::abs(m);
          if ((lx + ly - abs_m)%2)
            return 0.0;

          auto j = (lx + ly - abs_m)/2;
          if (j < 0)
            return 0.0;

          /*----------------------------------------------------------------------------------------
            Checking whether the cartesian polynomial contributes to the requested component of Ylm
           ----------------------------------------------------------------------------------------*/
          auto comp = (m >= 0) ? 1 : -1;
          /*  if (comp != ((abs_m-lx)%2 ? -1 : 1))*/
          auto i = abs_m-lx;
          if (comp != parity(abs(i)))
            return 0.0;

          assert(l <= 10); // libint2::math::fac[] is only defined up to 20
          Real pfac = sqrt(Real(fac[2*lx]*fac[2*ly]*fac[2*lz]*fac[l-abs_m])/(fac[2*l]*fac[l]*fac[lx]*fac[ly]*fac[lz]*fac[l+abs_m]));
          /*  pfac = sqrt(fac[l-abs_m]/(fac[l]*fac[l]*fac[l+abs_m]));*/
          pfac /= (1L << l);
          if (m < 0)
            pfac *= parity((i-1)/2);
          else
            pfac *= parity(i/2);

          auto i_min = j;
          auto i_max = (l-abs_m)/2;
          Real sum = 0;
          for(auto i=i_min;i<=i_max;i++) {
            Real pfac1 = bc(l,i)*bc(i,j);
            pfac1 *= (Real(parity(i)*fac[2*(l-i)])/fac[l-abs_m-2*i]);
            Real sum1 = 0.0;
            const int k_min = std::max((lx-abs_m)/2,0);
            const int k_max = std::min(j,lx/2);
            for(int k=k_min;k<=k_max;k++) {
              if (lx-2*k <= abs_m)
                sum1 += bc(j,k)*bc(abs_m,lx-2*k)*parity(k);
            }
            sum += pfac1*sum1;
          }
          sum *= sqrt(Real(df_Kminus1[2*l])/(df_Kminus1[2*lx]*df_Kminus1[2*ly]*df_Kminus1[2*lz]));

          if (m == 0)
            return pfac*sum;
          else
            return M_SQRT2*pfac*sum;
        }

    };

    using shg_coefs_type = SolidHarmonicsCoefficients<real_t>;
    static std::vector<shg_coefs_type> shg_coefs;

    inline void init() {
      shg_coefs.resize(LIBINT_MAX_AM+1);
      for(size_t i=0; i<shg_coefs.size(); ++i) {
        shg_coefs[i].init(static_cast<unsigned char>(i));
      }
    }
    inline void cleanup() {
      std::vector<shg_coefs_type> null;
      std::swap(shg_coefs, null);
    }

    // generic transforms

    template <typename Real>
    void transform_first(size_t l, size_t n2, const Real *src, Real *tgt)
    {
      const shg_coefs_type& coefs = shg_coefs[l];

      const auto n = 2*l+1;
      memset(tgt,0,n*n2*sizeof(Real));

      // loop over shg
      for(size_t s=0; s!=n; ++s) {
        const auto nc_s = coefs.nnz(s);      // # of cartesians contributing to shg s
        const auto* c_idxs = coefs.row_idx(s); // indices of cartesians contributing to shg s
        const auto* c_vals = coefs.row_values(s); // coefficients of cartesians contributing to shg s

        const auto tgt_blk_s_offset = s * n2;

        for(size_t ic=0; ic!=nc_s; ++ic) { // loop over contributing cartesians
          const auto c = c_idxs[ic];
          const auto s_c_coeff = c_vals[ic];

          auto src_blk_s = src + c * n2;
          auto tgt_blk_s = tgt + tgt_blk_s_offset;

          // loop over other dims
          for(size_t i2=0; i2!=n2; ++i2, ++src_blk_s, ++tgt_blk_s) {
            *tgt_blk_s += s_c_coeff * *src_blk_s;
          }
        }
      }

    }

    /// transforms two first dimensions of tensor from cartesian to real solid harmonic basis
    template <typename Real>
    void transform_first2(int l1, int l2, size_t inner_dim, const Real* source_blk, Real* target_blk) {
      const shg_coefs_type& coefs1 = shg_coefs[l1];
      const shg_coefs_type& coefs2 = shg_coefs[l2];

      const auto ncart1 = (l1+1)*(l1+2)/2;
      const auto ncart2 = (l2+1)*(l2+2)/2;
      const auto npure1 = 2*l1+1;
      const auto npure2 = 2*l2+1;
      const auto ncart2inner = ncart2 * inner_dim;
      const auto npure2inner = npure2 * inner_dim;
      memset(target_blk, 0, npure1*npure2inner*sizeof(real_t));

      // loop over blocks of inner dimension
      const size_t inner_blk_size = 8;
      const size_t nblks = (inner_dim+inner_blk_size-1)/inner_blk_size;
      for(size_t blk=0; blk!=nblks; ++blk) {
        const auto blk_begin = blk * inner_blk_size;
        const auto blk_end = std::min(blk_begin + inner_blk_size,inner_dim);
        const auto blk_size = blk_end - blk_begin;

        // loop over first shg
        for(size_t s1=0; s1!=npure1; ++s1) {
          const auto nc1 = coefs1.nnz(s1);      // # of cartesians contributing to shg s1
          const auto* c1_idxs = coefs1.row_idx(s1); // indices of cartesians contributing to shg s1
          const auto* c1_vals = coefs1.row_values(s1); // coefficients of cartesians contributing to shg s1

          auto target_blk_s1 = target_blk + s1 * npure2inner + blk_begin;

          // loop over second shg
          for(size_t s2=0; s2!=npure2; ++s2) {
            const auto nc2 = coefs2.nnz(s2);      // # of cartesians contributing to shg s2
            const auto* c2_idxs = coefs2.row_idx(s2); // indices of cartesians contributing to shg s2
            const auto* c2_vals = coefs2.row_values(s2); // coefficients of cartesians contributing to shg s2
            const auto s2inner = s2 * inner_dim;
            const auto target_blk_s1_blk_begin = target_blk_s1 + s2inner;

            for(size_t ic1=0; ic1!=nc1; ++ic1) { // loop over contributing cartesians
              auto c1 = c1_idxs[ic1];
              auto s1_c1_coeff = c1_vals[ic1];

              auto source_blk_c1 = source_blk + c1 * ncart2inner + blk_begin;

              for(size_t ic2=0; ic2!=nc2; ++ic2) { // loop over contributing cartesians
                auto c2 = c2_idxs[ic2];
                auto s2_c2_coeff = c2_vals[ic2];
                const auto c2inner = c2 * inner_dim;

                const auto coeff = s1_c1_coeff * s2_c2_coeff;
                const auto source_blk_c1_blk_begin = source_blk_c1 + c2inner;
                for(auto b=0; b<blk_size; ++b)
                  target_blk_s1_blk_begin[b] += source_blk_c1_blk_begin[b] * coeff;

              } // cart2

            } //cart1

          } // shg2

        } // shg1

      } // blk

    } // transform_first2()

    template <typename Real>
    void transform_inner(size_t n1, size_t l, size_t n2, const Real *src, Real *tgt)
    {
      const shg_coefs_type& coefs = shg_coefs[l];

      const auto nc = (l+1)*(l+2)/2;
      const auto n = 2*l+1;
      const auto nc_n2 = nc * n2;
      const auto n_n2 = n * n2;
      memset(tgt,0,n1*n_n2*sizeof(Real));

      // loop over shg
      for(size_t s=0; s!=n; ++s) {
        const auto nc_s = coefs.nnz(s);      // # of cartesians contributing to shg s
        const auto* c_idxs = coefs.row_idx(s); // indices of cartesians contributing to shg s
        const auto* c_vals = coefs.row_values(s); // coefficients of cartesians contributing to shg s

        const auto tgt_blk_s_offset = s * n2;

        for(size_t ic=0; ic!=nc_s; ++ic) { // loop over contributing cartesians
          const auto c = c_idxs[ic];
          const auto s_c_coeff = c_vals[ic];

          auto src_blk_s = src + c * n2;
          auto tgt_blk_s = tgt + tgt_blk_s_offset;

          // loop over other dims
          for(size_t i1=0; i1!=n1; ++i1, src_blk_s+=nc_n2, tgt_blk_s+=n_n2) {
            for(size_t i2=0; i2!=n2; ++i2) {
              tgt_blk_s[i2] += s_c_coeff * src_blk_s[i2];
            }
          }
        }
      }

    }

    /// transforms the last dimension of \c src from cartesian to solid harmonic Gaussians, stores result to \c tgt
    template <typename Real>
    void transform_last(size_t n1, size_t l, const Real *src, Real *tgt)
    {
      const shg_coefs_type& coefs = shg_coefs[l];

      const auto nc = (l+1)*(l+2)/2;
      const auto n = 2*l+1;
      memset(tgt,0,n1*n*sizeof(Real));

      // loop over shg
      for(size_t s=0; s!=n; ++s) {
        const auto nc_s = coefs.nnz(s);      // # of cartesians contributing to shg s
        const auto* c_idxs = coefs.row_idx(s); // indices of cartesians contributing to shg s
        const auto* c_vals = coefs.row_values(s); // coefficients of cartesians contributing to shg s

        const auto tgt_blk_s_offset = s;

        for(size_t ic=0; ic!=nc_s; ++ic) { // loop over contributing cartesians
          const auto c = c_idxs[ic];
          const auto s_c_coeff = c_vals[ic];

          auto src_blk_s = src + c;
          auto tgt_blk_s = tgt + tgt_blk_s_offset;

          // loop over other dims
          for(size_t i1=0; i1!=n1; ++i1, src_blk_s+=nc, tgt_blk_s+=n) {
              *tgt_blk_s += s_c_coeff * *src_blk_s;
          }
        }
      }

    }

    /// transforms the last two dimensions of \c src from cartesian to solid harmonic Gaussians, stores result to \c tgt
    template <typename Real>
    void tform_last2(size_t n1, int l_row, int l_col, const Real* source_blk, Real* target_blk) {
      const shg_coefs_type& coefs_row = shg_coefs[l_row];
      const shg_coefs_type& coefs_col = shg_coefs[l_col];

      const auto ncart_row = (l_row+1)*(l_row+2)/2;
      const auto ncart_col = (l_col+1)*(l_col+2)/2;
      const auto ncart = ncart_row * ncart_col;
      const auto npure_row = 2*l_row+1;
      const auto npure_col = 2*l_col+1;
      const auto npure = npure_row * npure_col;
      memset(target_blk, 0, n1*npure*sizeof(Real));

      for(size_t i1=0; i1!=n1; ++i1, source_blk+=ncart, target_blk+=npure) {
        // loop over row shg
        for(size_t s1=0; s1!=npure_row; ++s1) {
          const auto nc1 = coefs_row.nnz(s1);      // # of cartesians contributing to shg s1
          const auto* c1_idxs = coefs_row.row_idx(s1); // indices of cartesians contributing to shg s1
          const auto* c1_vals = coefs_row.row_values(s1); // coefficients of cartesians contributing to shg s1

          auto target_blk_s1 = target_blk + s1 * npure_col;

          // loop over col shg
          for(size_t s2=0; s2!=npure_col; ++s2) {
            const auto nc2 = coefs_col.nnz(s2);      // # of cartesians contributing to shg s2
            const auto* c2_idxs = coefs_col.row_idx(s2); // indices of cartesians contributing to shg s2
            const auto* c2_vals = coefs_col.row_values(s2); // coefficients of cartesians contributing to shg s2

            for(size_t ic1=0; ic1!=nc1; ++ic1) { // loop over contributing cartesians
              auto c1 = c1_idxs[ic1];
              auto s1_c1_coeff = c1_vals[ic1];

              auto source_blk_c1 = source_blk + c1 * ncart_col;

              for(size_t ic2=0; ic2!=nc2; ++ic2) { // loop over contributing cartesians
                auto c2 = c2_idxs[ic2];
                auto s2_c2_coeff = c2_vals[ic2];

                target_blk_s1[s2] += source_blk_c1[c2] * s1_c1_coeff * s2_c2_coeff;
              } // cart2

            } //cart1

          } // shg2

        } // shg1
      }
    } // tform()

    /// multiplies rows and columns of matrix \c source_blk, stores result to \c target_blk
    template <typename Real>
    void tform(int l_row, int l_col, const Real* source_blk, Real* target_blk) {
      const shg_coefs_type& coefs_row = shg_coefs[l_row];
      const shg_coefs_type& coefs_col = shg_coefs[l_col];

      const auto ncart_row = (l_row+1)*(l_row+2)/2;
      const auto ncart_col = (l_col+1)*(l_col+2)/2;
      const auto npure_row = 2*l_row+1;
      const auto npure_col = 2*l_col+1;
      memset(target_blk, 0, npure_row*npure_col*sizeof(real_t));

      // loop over row shg
      for(size_t s1=0; s1!=npure_row; ++s1) {
        const auto nc1 = coefs_row.nnz(s1);      // # of cartesians contributing to shg s1
        const auto* c1_idxs = coefs_row.row_idx(s1); // indices of cartesians contributing to shg s1
        const auto* c1_vals = coefs_row.row_values(s1); // coefficients of cartesians contributing to shg s1

        auto target_blk_s1 = target_blk + s1 * npure_col;

        // loop over col shg
        for(size_t s2=0; s2!=npure_col; ++s2) {
          const auto nc2 = coefs_col.nnz(s2);      // # of cartesians contributing to shg s2
          const auto* c2_idxs = coefs_col.row_idx(s2); // indices of cartesians contributing to shg s2
          const auto* c2_vals = coefs_col.row_values(s2); // coefficients of cartesians contributing to shg s2

          for(size_t ic1=0; ic1!=nc1; ++ic1) { // loop over contributing cartesians
            auto c1 = c1_idxs[ic1];
            auto s1_c1_coeff = c1_vals[ic1];

            auto source_blk_c1 = source_blk + c1 * ncart_col;

            for(size_t ic2=0; ic2!=nc2; ++ic2) { // loop over contributing cartesians
              auto c2 = c2_idxs[ic2];
              auto s2_c2_coeff = c2_vals[ic2];

              target_blk_s1[s2] += source_blk_c1[c2] * s1_c1_coeff * s2_c2_coeff;
            } // cart2

          } //cart1

        } // shg2

      } // shg1
    } // transform_last2()

    /// multiplies columns of matrix \c source_blk, stores result to \c target_blk
    template <typename Real>
    void tform_cols(size_t nrow, int l_col, const Real* source_blk, Real* target_blk) {
      return transform_last(nrow, l_col, source_blk, target_blk);
      const shg_coefs_type& coefs_col = shg_coefs[l_col];

      const auto ncart_col = (l_col+1)*(l_col+2)/2;
      const auto npure_col = 2*l_col+1;

      // loop over rows
      for(size_t r1=0; r1!=nrow; ++r1) {

        auto source_blk_r1 = source_blk + r1 * ncart_col;
        auto target_blk_r1 = target_blk + r1 * npure_col;

        // loop over col shg
        for(size_t s2=0; s2!=npure_col; ++s2) {
          const auto nc2 = coefs_col.nnz(s2);      // # of cartesians contributing to shg s2
          const auto* c2_idxs = coefs_col.row_idx(s2); // indices of cartesians contributing to shg s2
          const auto* c2_vals = coefs_col.row_values(s2); // coefficients of cartesians contributing to shg s2

          Real r1_s2_value = 0.0;

          for(size_t ic2=0; ic2!=nc2; ++ic2) { // loop over contributing cartesians
            auto c2 = c2_idxs[ic2];
            auto s2_c2_coeff = c2_vals[ic2];

            r1_s2_value += source_blk_r1[c2] * s2_c2_coeff;

          } // cart2

          target_blk_r1[s2] = r1_s2_value;

        } // shg1

      } // rows

    } // tform_cols()

    /// multiplies rows of matrix \c source_blk, stores result to \c target_blk
    template <typename Real>
    void tform_rows(int l_row, size_t ncol,  const Real* source_blk, Real* target_blk) {
      return transform_first(l_row, ncol, source_blk, target_blk);
      const shg_coefs_type& coefs_row = shg_coefs[l_row];

      const auto ncart_row = (l_row+1)*(l_row+2)/2;
      const auto npure_row = 2*l_row+1;

      // loop over row shg
      for(size_t s1=0; s1!=npure_row; ++s1) {
        const auto nc1 = coefs_row.nnz(s1);      // # of cartesians contributing to shg s1
        const auto* c1_idxs = coefs_row.row_idx(s1); // indices of cartesians contributing to shg s1
        const auto* c1_vals = coefs_row.row_values(s1); // coefficients of cartesians contributing to shg s1

        auto target_blk_s1 = target_blk + s1 * ncol;

        // loop over cols
        for(size_t c2=0; c2!=ncol; ++c2) {

          Real s1_c2_value = 0.0;
          auto source_blk_c2_offset = source_blk + c2;

          for(size_t ic1=0; ic1!=nc1; ++ic1) { // loop over contributing cartesians
            auto c1 = c1_idxs[ic1];
            auto s1_c1_coeff = c1_vals[ic1];

            s1_c2_value += source_blk_c2_offset[c1 * ncol] * s1_c1_coeff;

          } //cart1

          target_blk_s1[c2] = s1_c2_value;

        } // shg2

      } // shg1
    } // tform_rows();

    /// transforms matrix from cartesian to real solid harmonic basis
    template <typename Real, typename Shell> // Shell = libint2::Shell::Contraction
    void tform(const Shell& shell_row, const Shell& shell_col, const Real* source_blk, Real* target_blk) {
      const auto trow = shell_row.pure;
      const auto tcol = shell_col.pure;
      if (trow) {
        if (tcol) {
          //tform(shell_row.l, shell_col.l, source_blk, target_blk);
          Real localscratch[500];
          tform_cols(shell_row.cartesian_size(), shell_col.l, source_blk, &localscratch[0]);
          tform_rows(shell_row.l, shell_col.size(), &localscratch[0], target_blk);
        }
        else
          tform_rows(shell_row.l, shell_col.cartesian_size(), source_blk, target_blk);
      }
      else
        tform_cols(shell_row.cartesian_size(), shell_col.l, source_blk, target_blk);
    }

  } // namespace libint2::solidharmonics

  inline void libint2_init_shg() {
    libint2::solidharmonics::init();
  }
  inline void libint2_cleanup_shg() {
    libint2::solidharmonics::cleanup();
  }

} // namespace libint2

#endif /* _libint2_src_lib_libint_solidharmonics_h_ */
