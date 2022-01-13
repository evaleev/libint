/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_lib_libint_diis_h_
#define _libint2_src_lib_libint_diis_h_

#include <libint2/util/cxxstd.h>
#if LIBINT2_CPLUSPLUS_STD < 2011
# error "libint2/diis.h requires C++11 support"
#endif

#include <deque>

#pragma GCC diagnostic push
#pragma GCC system_header
#include <Eigen/Core>
#include <Eigen/QR>
#pragma GCC diagnostic pop

namespace libint2 {

  namespace diis {

    template <typename D>
    struct traits;

    /*
    template <typename D>
    typename traits<D>::element_type
    dot_product(const D& d1, const D& d2);

    template <typename D>
    void
    zero(D& d);

    template <typename D, typename Scalar>
    void
    axpy(const D& y, Scalar a, const D& x);
    */

    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    struct traits<Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >> {
        typedef Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols > D;
        typedef typename D::Scalar element_type;
    };

    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    typename traits<Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >>::element_type
    dot_product(const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >& d1,
                const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >& d2) {
      return d1.cwiseProduct(d2).sum();
    }

    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols>
    void
    zero(Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >& d) {
      d.setZero(d.rows(), d.cols());
    }

    template<typename _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols, typename Scalar>
    void
    axpy(Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >& y,
         Scalar a,
         const Eigen::Matrix< _Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols >& x) {
      y += a*x;
    }

  }

  /// DIIS (``direct inversion of iterative subspace'') extrapolation

  /// The DIIS class provides DIIS extrapolation to an iterative solver of
  /// (systems of) linear or nonlinear equations of the \f$ f(x) = 0 \f$ form,
  /// where \f$ f(x) \f$ is a (non-linear) function of \f$ x \f$ (in general,
  /// \f$ x \f$ is a set of numeric values). Such equations are usually solved
  /// iteratively as follows:
  /// \li given a current guess at the solution, \f$ x_i \f$, evaluate the error
  ///     (``residual'') \f$ e_i = f(x_i) \f$ (NOTE that the dimension of
  ///     \f$ x \f$ and \f$ e \f$ do not need to coincide);
  /// \li use the error to compute an updated guess \f$ x_{i+1} = x_i + g(e_i) \f$;
  /// \li proceed until a norm of the error is less than the target precision
  ///     \f$ \epsilon \f$. Another convergence criterion may include
  ///     \f$ ||x_{i+1} - x_i|| < \epsilon \f$ .
  /// \\
  /// For example, in the Hartree-Fock method in the density form, one could
  /// choose \f$ x \equiv \mathbf{P} \f$, the one-electron density matrix, and
  /// \f$ f(\mathbf{P}) \equiv [\mathbf{F}, \mathbf{P}] \f$ , where
  /// \f$ \mathbf{F} = \mathbf{F}(\mathbf{P}) \f$ is the Fock matrix, a linear
  /// function of the density. Because \f$ \mathbf{F} \f$ is a linear function
  /// of the density and DIIS uses a linear extrapolation, it is possible to
  /// just extrapolate the Fock matrix itself, i.e. \f$ x \equiv \mathbf{F} \f$
  /// and \f$ f(\mathbf{F}) \equiv [\mathbf{F}, \mathbf{P}] \f$ .
  /// \\
  /// Similarly, in the Hartree-Fock method in the molecular orbital
  /// representation, DIIS is used to extrapolate the Fock matrix, i.e.
  /// \f$ x \equiv \mathbf{F} \f$ and \f$ f(\mathbf{F}) \equiv \{ F_i^a \} \f$ ,
  /// where \f$ i \f$ and \f$ a \f$ are the occupied and unoccupied orbitals,
  /// respectively.
  /// \\
  /// Here's a short description of the DIIS method. Given a set of solution
  /// guess vectors \f$ \{ x_k \}, k=0..i \f$ and the corresponding error
  /// vectors \f$ \{ e_k \} \f$ DIIS tries to find a linear combination of
  /// \f$ x \f$ that would minimize the error by solving a simple linear system
  /// set up from the set of errors. The solution is a vector of coefficients
  /// \f$ \{ C_k \} \f$ that can be used to obtain an improved \f$ x \f$:
  /// \f$ x_{\mathrm{extrap},i+1} = \sum\limits_{k=0}^i C_{k,i} x_{k} \f$
  /// A more complicated version of DIIS introduces mixing:
  /// \f$ x_{\mathrm{extrap},i+1} = \sum\limits_{k=0}^i C_{k,i} ( (1-f) x_{k} + f x_{extrap,k} ) \f$
  /// Note that the mixing is not used in the first iteration.
  /// \\
  /// The original DIIS reference: P. Pulay, Chem. Phys. Lett. 73, 393 (1980).
  ///
  /// \tparam D type of \c x
  template <typename D>
  class DIIS {
    public:
      typedef typename diis::traits<D>::element_type value_type;

      /// Constructor

      /// \param strt The DIIS extrapolation will begin on the iteration given
      ///   by this integer (default = 1).
      /// \param ndi This integer maximum number of data sets to retain (default
      ///   = 5).
      /// \param dmp This nonnegative floating point number is used to dampen
      ///   the DIIS extrapolation (default = 0.0).
      /// \param ngr The number of iterations in a DIIS group. DIIS
      ///   extrapolation is only used for the first \c ngrdiis of these
      ///   iterations (default = 1). If \c ngr is 1 and \c ngrdiis is
      ///   greater than 0, then DIIS will be used on all iterations after and
      ///   including the start iteration.
      /// \param ngrdiis The number of DIIS extrapolations to do at the
      ///   beginning of an iteration group.  See the documentation for \c ngr
      ///   (default = 1).
      /// \param mf This real number in [0,1] is used to dampen the DIIS
      ///   extrapolation by mixing the input data with the output data for each
      ///   iteration (default = 0.0), which performs no mixing. The approach
      ///   described in Kerker, Phys. Rev. B, 23, p3082, 1981.
      DIIS(unsigned int strt=1,
           unsigned int ndi=5,
           value_type dmp =0,
           unsigned int ngr=1,
           unsigned int ngrdiis=1,
           value_type mf=0) :
             error_(0), errorset_(false),
             start(strt), ndiis(ndi),
             iter(0), ngroup(ngr),
             ngroupdiis(ngr),
             damping_factor(dmp),
             mixing_fraction(mf)
           {
            init();
           }
      ~DIIS() {
        x_.clear();
        errors_.clear();
        x_extrap_.clear();
      }

      /// \param[in,out] x On input, the most recent solution guess; on output,
      ///   the extrapolated guess
      /// \param[in,out] error On input, the most recent error; on output, the
      ///   if \c extrapolate_error \c == \c true will be the extrapolated
      ///   error, otherwise the value unchanged
      /// \param extrapolate_error whether to extrapolate the error (default =
      ///   false).
      void extrapolate(D& x,
                       D& error,
                       bool extrapolate_error = false)
      {
        using namespace ::libint2::diis;

        const value_type zero_determinant = std::numeric_limits<value_type>::epsilon();
        const value_type zero_norm = 1.0e-10;

        iter++;

        const bool do_mixing = (mixing_fraction != 0.0);
        const value_type scale = 1.0 + damping_factor;

        // if have ndiis vectors
        if (errors_.size() == ndiis) { // holding max # of vectors already? drop the least recent {x, error} pair
          x_.pop_front();
          errors_.pop_front();
          if (not x_extrap_.empty()) x_extrap_.pop_front();
          EigenMatrixX Bcrop = B_.bottomRightCorner(ndiis-1,ndiis-1);
          Bcrop.conservativeResize(ndiis,ndiis);
          B_ = Bcrop;
        }

        // push {x, error} to the set
        x_.push_back(x);
        errors_.push_back(error);
        const unsigned int nvec = errors_.size();
        assert(x_.size() == errors_.size());

        // and compute the most recent elements of B, B(i,j) = <ei|ej>
        for (unsigned int i=0; i < nvec-1; i++)
          B_(i,nvec-1) = B_(nvec-1,i) = dot_product(errors_[i], errors_[nvec-1]);
        B_(nvec-1,nvec-1) = dot_product(errors_[nvec-1], errors_[nvec-1]);

        if (iter == 1) { // the first iteration
          if (not x_extrap_.empty() && do_mixing) {
            zero(x);
            axpy(x, (1.0-mixing_fraction), x_[0]);
            axpy(x, mixing_fraction, x_extrap_[0]);
          }
        }
        else if (iter > start && (((iter - start) % ngroup) < ngroupdiis)) { // not the first iteration and need to extrapolate?

          EigenVectorX c;

          value_type absdetA;
          unsigned int nskip = 0; // how many oldest vectors to skip for the sake of conditioning?
                                        // try zero
          // skip oldest vectors until found a numerically stable system
          do {

            const unsigned int rank = nvec - nskip + 1; // size of matrix A

            // set up the DIIS linear system: A c = rhs
            EigenMatrixX A(rank, rank);
            A.col(0).setConstant(-1.0);
            A.row(0).setConstant(-1.0);
            A(0,0) = 0.0;
            EigenVectorX rhs = EigenVectorX::Zero(rank);
            rhs[0] = -1.0;

            value_type norm = 1.0;
            if (B_(nskip,nskip) > zero_norm)
              norm = 1.0/B_(nskip,nskip);

            A.block(1, 1, rank-1, rank-1) = B_.block(nskip, nskip, rank-1, rank-1) * norm;
            A.diagonal() *= scale;
            //for (unsigned int i=1; i < rank ; i++) {
            //  for (unsigned int j=1; j <= i ; j++) {
            //    A(i, j) = A(j, i) = B_(i+nskip-1, j+nskip-1) * norm;
            //    if (i==j) A(i, j) *= scale;
            //  }
            //}

#if 0
            std::cout << "DIIS: iter=" << iter << " nskip=" << nskip << " nvec=" << nvec << std::endl;
            std::cout << "DIIS: B=" << B_ << std::endl;
            std::cout << "DIIS: A=" << A << std::endl;
            std::cout << "DIIS: rhs=" << rhs << std::endl;
#endif

            // finally, solve the DIIS linear system
            Eigen::ColPivHouseholderQR<EigenMatrixX> A_QR = A.colPivHouseholderQr();
            c = A_QR.solve(rhs);
            absdetA = A_QR.absDeterminant();

            //std::cout << "DIIS: |A|=" << absdetA << " sol=" << c << std::endl;

            ++nskip;

          } while (absdetA < zero_determinant && nskip < nvec); // while (system is poorly conditioned)

          // failed?
          if (absdetA < zero_determinant) {
            std::ostringstream oss;
            oss << "DIIS::extrapolate: poorly-conditioned system, |A| = " << absdetA;
            throw std::domain_error(oss.str());
          }
          --nskip; // undo the last ++ :-(

          {

            zero(x);
            for (unsigned int k=nskip, kk=1; k < nvec; ++k, ++kk) {
              if (not do_mixing || x_extrap_.empty()) {
                //std::cout << "contrib " << k << " c=" << c[kk] << ":" << std::endl << x_[k] << std::endl;
                axpy(x, c[kk], x_[k]);
                if (extrapolate_error)
                  axpy(error, c[kk], errors_[k]);
              } else {
                axpy(x, c[kk] * (1.0 - mixing_fraction), x_[k]);
                axpy(x, c[kk] * mixing_fraction, x_extrap_[k]);
              }
            }
          }
        } // do DIIS

        // only need to keep extrapolated x if doing mixing
        if (do_mixing) x_extrap_.push_back(x);
      }

      /// calling this function forces the extrapolation to start upon next call
      /// to \c extrapolate() even if this object was initialized with start
      /// value greater than the current iteration index.
      void start_extrapolation() {
        if (start > iter) start = iter+1;
      }

      void reinitialize(const D* data = 0) {
        iter=0;
        if (data) {
          const bool do_mixing = (mixing_fraction != 0.0);
          if (do_mixing) x_extrap_.push_front(*data);
        }
      }

    private:
      value_type error_;
      bool errorset_;

      unsigned int start;
      unsigned int ndiis;
      unsigned int iter;
      unsigned int ngroup;
      unsigned int ngroupdiis;
      value_type damping_factor;
      value_type mixing_fraction;

      typedef Eigen::Matrix<value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixX;
      typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> EigenVectorX;

      EigenMatrixX B_; //!< B(i,j) = <ei|ej>

      std::deque<D> x_; //!< set of most recent x given as input (i.e. not exrapolated)
      std::deque<D> errors_; //!< set of most recent errors
      std::deque<D> x_extrap_; //!< set of most recent extrapolated x

      void set_error(value_type e) { error_ = e; errorset_ = true; }
      value_type error() { return error_; }

      void init() {
        iter = 0;

        B_ = EigenMatrixX::Zero(ndiis,ndiis);

        x_.clear();
        errors_.clear();
        x_extrap_.clear();
        //x_.resize(ndiis);
        //errors_.resize(ndiis);
        // x_extrap_ is bigger than the other because
        // it must hold data associated with the next iteration
        //x_extrap_.resize(diis+1);
      }

  }; // class DIIS

} // namespace libint2

#include <libint2/engine.h>

#endif /* _libint2_src_lib_libint_diis_h_ */
