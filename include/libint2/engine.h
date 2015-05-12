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

#ifndef _libint2_src_lib_libint_engine_h_
#define _libint2_src_lib_libint_engine_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>
#include <map>

#include <libint2.h>
#include <libint2/boys.h>
#include <libint2/shell.h>
#include <libint2/timer.h>
#include <libint2/solidharmonics.h>

#include <Eigen/Core>

// the engine will be profiled by default if library was configured with --enable-profile
#ifdef LIBINT2_PROFILE
#  define LIBINT2_ENGINE_TIMERS
// uncomment if want to profile each integral class
#  define LIBINT2_ENGINE_PROFILE_CLASS
#endif
// uncomment if want to profile the engine even if library was configured without --enable-profile
//#  define LIBINT2_ENGINE_TIMERS

namespace libint2 {

#ifdef LIBINT2_SUPPORT_ONEBODY
  /// OneBodyEngine computes integrals of 1-body operators, e.g. overlap, kinetic energy, dipole moment, etc.

  /**
   * OneBodyEngine computes integrals of types given by OneBodyEngine::integral_type
   */
  class OneBodyEngine {
    public:
      enum integral_type {
        overlap,
        kinetic,
        nuclear,
        _invalid
      };

      typedef libint2::FmEval_Taylor<real_t, 7> coulomb_core_eval_t;

      /// creates a default (unusable) OneBodyEngine; to be used as placeholder for copying a usable engine
      OneBodyEngine() : type_(_invalid), primdata_(), lmax_(-1) {}

      /// Constructs a (usable) OneBodyEngine

      /// \param t integral type, see OneBodyEngine::integral_type
      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_level if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_level
      /// \note if integral_type == nuclear, must specify charges using set_q()
      /// \warning currently only the following integral types are suported: \c overlap, \c kinetic, \c nuclear
      /// \warning currently derivative integrals are not supported
      /// \warning currently solid harmonics Gaussians are not supported
      OneBodyEngine(integral_type t, size_t max_nprim, int max_l, int deriv_order = 0) :
        type_(t), primdata_(max_nprim * max_nprim), lmax_(max_l), deriv_order_(deriv_order),
        fm_eval_(t == nuclear ? coulomb_core_eval_t::instance(2*max_l+deriv_order, 1e-25) : 0)
        //fm_eval_(0)
      {
        initialize();
      }

      /// move constructor
      OneBodyEngine(OneBodyEngine&& other) = default;

      /// (deep) copy constructor
      OneBodyEngine(const OneBodyEngine& other) :
        type_(other.type_),
        primdata_(other.primdata_.size()),
        lmax_(other.lmax_),
        deriv_order_(other.deriv_order_),
        q_(other.q_),
        fm_eval_(other.fm_eval_) {
        initialize();
      }

      ~OneBodyEngine() {
        finalize();
      }

      /// move assignment
      OneBodyEngine& operator=(OneBodyEngine&& other) = default;

      /// (deep) copy assignment
      OneBodyEngine& operator=(const OneBodyEngine& other) {
        type_ = other.type_;
        primdata_.resize(other.primdata_.size());
        lmax_ = other.lmax_;
        deriv_order_ = other.deriv_order_;
        q_ = other.q_;
        fm_eval_ = other.fm_eval_;
        initialize();
        return *this;
      }

      /// returns the integral type used by the engine
      integral_type type() const {return type_;}

      /// specifies the nuclear charges
      /// \param q vector of {charge,Cartesian coordinate} pairs
      void set_q(const std::vector<std::pair<double, std::array<double, 3>>>& q) {
        q_ = q;
      }

      /// computes shell set of integrals
      /// \note result is stored in row-major order
      const real_t* compute(const libint2::Shell& s1,
                            const libint2::Shell& s2) {

        // can only handle 1 contraction at a time
        assert(s1.ncontr() == 1 && s2.ncontr() == 1);
        // derivatives not supported for now
        assert(deriv_order_ == 0);

        const auto l1 = s1.contr[0].l;
        const auto l2 = s2.contr[0].l;

        // if want nuclear, make sure there is at least one nucleus .. otherwise the user likely forgot to call set_q
        if (type_ == nuclear and q_.size() == 0)
          throw std::runtime_error("libint2::OneBodyEngine(integral_type = nuclear), but no nuclei found; forgot to call set_q()?");

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // make sure bra.l >= ket.l
        const auto swap = (l1 < l2);
#else // make sure bra.l <= ket.l
        const auto swap = (l1 > l2);
#endif
        const auto& bra = !swap ? s1 : s2;
        const auto& ket = !swap ? s2 : s1;

        const auto n1 = s1.size();
        const auto n2 = s2.size();
        const auto ncart1 = s1.cartesian_size();
        const auto ncart2 = s2.cartesian_size();

        const bool use_scratch = (swap || type_ == nuclear);

        // assert # of primitive pairs
        const auto nprim_bra = bra.nprim();
        const auto nprim_ket = ket.nprim();
        const auto nprimpairs = nprim_bra * nprim_ket;
        assert(nprimpairs <= primdata_.size());

        // adjust max angular momentum, if needed
        const auto lmax = std::max(l1, l2);
        assert (lmax <= lmax_);
        if (lmax == 0) // (s|s) ints will be accumulated in the first element of stack
          primdata_[0].stack[0] = 0;
        else if (use_scratch)
          memset(static_cast<void*>(&scratch_[0]), 0, sizeof(real_t)*ncart1*ncart2);

        // loop over operator components
        const auto num_operset = type_ == nuclear ? q_.size() : 1u;
        for(auto oset=0u; oset!=num_operset; ++oset) {

            auto p12 = 0;
            for(auto pb=0; pb!=nprim_bra; ++pb) {
              for(auto pk=0; pk!=nprim_ket; ++pk, ++p12) {
                compute_primdata(primdata_[p12],bra,ket,pb,pk,oset);
              }
            }
            primdata_[0].contrdepth = p12;

            if (lmax == 0 && type_ != kinetic) { // (s|s) or (s|V|s)
              auto& result = primdata_[0].stack[0];
              switch (type_) {
                case overlap:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12]._0_Overlap_0_x[0]
                          * primdata_[p12]._0_Overlap_0_y[0]
                          * primdata_[p12]._0_Overlap_0_z[0];
                  break;
                case nuclear:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12].LIBINT_T_S_ELECPOT_S(0)[0];
                  break;
                default:
                  assert(false);
              }
              primdata_[0].targets[0] = primdata_[0].stack;
            }
            else {
              switch (type_) {
                case overlap:
                  LIBINT2_PREFIXED_NAME(libint2_build_overlap)[bra.contr[0].l][ket.contr[0].l](&primdata_[0]);
                  break;
                case kinetic:
                  LIBINT2_PREFIXED_NAME(libint2_build_kinetic)[bra.contr[0].l][ket.contr[0].l](&primdata_[0]);
                  break;
                case nuclear:
                  LIBINT2_PREFIXED_NAME(libint2_build_elecpot)[bra.contr[0].l][ket.contr[0].l](&primdata_[0]);
                  break;
                default:
                  assert(false);
              }
              if (use_scratch) {
                const auto ncart_bra = bra.cartesian_size();
                const auto ncart_ket = ket.cartesian_size();
                constexpr auto using_scalar_real = std::is_same<double,real_t>::value || std::is_same<float,real_t>::value;
                static_assert(using_scalar_real, "Libint2 C++11 API only supports fundamental real types");
                typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;
                Eigen::Map<Matrix> braket(primdata_[0].targets[0], ncart_bra, ncart_ket);
                Eigen::Map<Matrix> set12(&scratch_[0], ncart1, ncart2);
                if (swap)
                  set12 += braket.transpose();
                else
                  set12 += braket;
              }
            } // ltot != 0

        } // oset (operator components, artifact of nuclear)

        auto cartesian_ints = (use_scratch && lmax != 0) ? &scratch_[0] : primdata_[0].targets[0];

        auto result = cartesian_ints;

        if (s1.contr[0].pure || s2.contr[0].pure) {
          auto* spherical_ints = (cartesian_ints == &scratch_[0]) ? primdata_[0].targets[0] : &scratch_[0];
          if (s1.contr[0].pure && s2.contr[0].pure) {
            libint2::solidharmonics::tform(l1, l2, cartesian_ints, spherical_ints);
          }
          else {
            if (s1.contr[0].pure)
              libint2::solidharmonics::tform_rows(l1, n2, cartesian_ints, spherical_ints);
            else
              libint2::solidharmonics::tform_cols(n1, l2, cartesian_ints, spherical_ints);
          }

          result = spherical_ints;
        } // tform to solids

        return result;
      }

      void compute_primdata(Libint_t& primdata,
                            const Shell& s1, const Shell& s2,
                            size_t p1, size_t p2,
                            size_t oset) {

        const auto& A = s1.O;
        const auto& B = s2.O;

        const auto alpha1 = s1.alpha[p1];
        const auto alpha2 = s2.alpha[p2];

        const auto c1 = s1.contr[0].coeff[p1];
        const auto c2 = s2.contr[0].coeff[p2];

        const auto gammap = alpha1 + alpha2;
        const auto oogammap = 1.0 / gammap;
        const auto rhop = alpha1 * alpha2 * oogammap;
        const auto Px = (alpha1 * A[0] + alpha2 * B[0]) * oogammap;
        const auto Py = (alpha1 * A[1] + alpha2 * B[1]) * oogammap;
        const auto Pz = (alpha1 * A[2] + alpha2 * B[2]) * oogammap;
        const auto AB_x = A[0] - B[0];
        const auto AB_y = A[1] - B[1];
        const auto AB_z = A[2] - B[2];
        const auto AB2_x = AB_x * AB_x;
        const auto AB2_y = AB_y * AB_y;
        const auto AB2_z = AB_z * AB_z;
        const auto AB2 = AB2_x + AB2_y + AB2_z;

        assert (LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD);

        // overlap and kinetic energy ints don't use HRR, hence VRR on both centers
        // Coulomb potential do HRR on center 1 only
#if LIBINT2_DEFINED(eri,PA_x)
          primdata.PA_x[0] = Px - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
          primdata.PA_y[0] = Py - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
          primdata.PA_z[0] = Pz - A[2];
#endif

        if (type_ == overlap || type_ == kinetic ) {

#if LIBINT2_DEFINED(eri,PB_x)
          primdata.PB_x[0] = Px - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
          primdata.PB_y[0] = Py - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
          primdata.PB_z[0] = Pz - B[2];
#endif
        }

#if LIBINT2_DEFINED(eri,oo2z)
        primdata.oo2z[0] = 0.5*oogammap;
#endif

        if (type_ == nuclear) { // additional factor for electrostatic potential
          const auto& C = q_[oset].second;
#if LIBINT2_DEFINED(eri,PC_x)
          primdata.PC_x[0] = Px - C[0];
#endif
#if LIBINT2_DEFINED(eri,PC_y)
          primdata.PC_y[0] = Py - C[1];
#endif
#if LIBINT2_DEFINED(eri,PC_z)
          primdata.PC_z[0] = Pz - C[2];
#endif
          // elecpot uses HRR
#if LIBINT2_DEFINED(eri,AB_x)
          primdata.AB_x[0] = A[0] - B[0];
#endif
#if LIBINT2_DEFINED(eri,AB_y)
          primdata.AB_y[0] = A[1] - B[1];
#endif
#if LIBINT2_DEFINED(eri,AB_z)
          primdata.AB_z[0] = A[2] - B[2];
#endif

        }

        if (deriv_order_ > 0) {
          // prefactors for derivative overlap relations
          assert(false);
        }

        decltype(c1) sqrt_PI(1.77245385090551602729816748334);
        const auto xyz_pfac = sqrt_PI * sqrt(oogammap);
        const auto ovlp_ss_x = exp(- rhop * AB2_x) * xyz_pfac * c1 * c2;
        const auto ovlp_ss_y = exp(- rhop * AB2_y) * xyz_pfac;
        const auto ovlp_ss_z = exp(- rhop * AB2_z) * xyz_pfac;

        primdata._0_Overlap_0_x[0] = ovlp_ss_x;
        primdata._0_Overlap_0_y[0] = ovlp_ss_y;
        primdata._0_Overlap_0_z[0] = ovlp_ss_z;

        if (type_ == kinetic) {
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
          primdata.two_alpha0_bra[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
          primdata.two_alpha0_ket[0] = 2.0 * alpha2;
#endif
        }

        if (type_ == nuclear) {
#if LIBINT2_DEFINED(eri,PC_x) && LIBINT2_DEFINED(eri,PC_y) && LIBINT2_DEFINED(eri,PC_z)
          const auto PC2 = primdata.PC_x[0] * primdata.PC_x[0] +
                           primdata.PC_y[0] * primdata.PC_y[0] +
                           primdata.PC_z[0] * primdata.PC_z[0];
          const auto U = gammap * PC2;
          const auto ltot = s1.contr[0].l + s2.contr[0].l;
          auto* fm_ptr = &(primdata.LIBINT_T_S_ELECPOT_S(0)[0]);
          fm_eval_->eval(fm_ptr, U, ltot);

          decltype(U) two_o_sqrt_PI(1.12837916709551257389615890312);
          const auto pfac = - q_[oset].first * sqrt(gammap) * two_o_sqrt_PI * ovlp_ss_x * ovlp_ss_y * ovlp_ss_z;
          const auto ltot_p1 = ltot + 1;
          for(auto m=0; m!=ltot_p1; ++m) {
            fm_ptr[m] *= pfac;
          }
#endif
        }

      }

    private:
      integral_type type_;
      std::vector<Libint_t> primdata_;
      int lmax_;
      size_t deriv_order_;
      std::vector<std::pair<double, std::array<double,3>>> q_;

      std::shared_ptr<coulomb_core_eval_t> fm_eval_;

      std::vector<real_t> scratch_; // for transposes and/or transforming to solid harmonics

      void initialize() {
        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;

        switch(type_) {
          case overlap: assert(lmax_ <= LIBINT2_MAX_AM_overlap); break;
          case kinetic: assert(lmax_ <= LIBINT2_MAX_AM_kinetic); break;
          case nuclear: assert(lmax_ <= LIBINT2_MAX_AM_elecpot); break;
          default: assert(false);
        }
        assert(deriv_order_ <= LIBINT2_DERIV_ONEBODY_ORDER);

        if (type_ == overlap) {
          switch(deriv_order_) {

            case 0:
              libint2_init_overlap(&primdata_[0], lmax_, 0);
              scratch_.resize(ncart_max*ncart_max);
              break;
            case 1:
#if LIBINT2_DERIV_ONEBODY_ORDER > 0
              libint2_init_overlap1(&primdata_[0], lmax_, 0);
              scratch_.resize(3 * ncart_max*ncart_max);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
              libint2_init_overlap2(&primdata_[0], lmax_, 0);
              scratch_.resize(6 * ncart_max*ncart_max);
#endif
              break;
            default: assert(deriv_order_ < 3);
          }

          return;
        }

        if (type_ == kinetic) {
          switch(deriv_order_) {

            case 0:
              libint2_init_kinetic(&primdata_[0], lmax_, 0);
              scratch_.resize(ncart_max*ncart_max);
              break;
            case 1:
#if LIBINT2_DERIV_ONEBODY_ORDER > 0
              libint2_init_kinetic1(&primdata_[0], lmax_, 0);
              scratch_.resize(3 * ncart_max*ncart_max);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
              libint2_init_kinetic2(&primdata_[0], lmax_, 0);
              scratch_.resize(6 * ncart_max*ncart_max);
#endif
              break;
            default: assert(deriv_order_ < 3);
          }

          return;
        }

        if (type_ == nuclear) {

          switch(deriv_order_) {

            case 0:
              libint2_init_elecpot(&primdata_[0], lmax_, 0);
              scratch_.resize(ncart_max*ncart_max); // one more set to be able to accumulate
              break;
            case 1:
#if LIBINT2_DERIV_ONEBODY_ORDER > 0
              libint2_init_elecpot1(&primdata_[0], lmax_, 0);
              scratch_.resize(3 * ncart_max*ncart_max);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
              libint2_init_elecpot2(&primdata_[0], lmax_, 0);
              scratch_.resize(6 * ncart_max*ncart_max);
#endif
              break;
            default: assert(deriv_order_ < 3);
          }

          return;
        }

        assert(type_ == overlap || type_ == kinetic || type_ == nuclear);
      } // initialize()

      void finalize() {
        if (primdata_.size() != 0) {

          if (type_ == overlap) {
            switch(deriv_order_) {

              case 0:
              libint2_cleanup_overlap(&primdata_[0]);
              break;
              case 1:
  #if LIBINT2_DERIV_ONEBODY_ORDER > 0
              libint2_cleanup_overlap1(&primdata_[0]);
  #endif
              break;
              case 2:
  #if LIBINT2_DERIV_ONEBODY_ORDER > 1
              libint2_cleanup_overlap2(&primdata_[0]);
  #endif
              break;
            }

            return;
          }

          if (type_ == kinetic) {
            switch(deriv_order_) {

              case 0:
              libint2_cleanup_kinetic(&primdata_[0]);
              break;
              case 1:
  #if LIBINT2_DERIV_ONEBODY_ORDER > 0
              libint2_cleanup_kinetic1(&primdata_[0]);
  #endif
              break;
              case 2:
  #if LIBINT2_DERIV_ONEBODY_ORDER > 1
              libint2_cleanup_kinetic2(&primdata_[0]);
  #endif
              break;
            }

            return;
          }

          if (type_ == nuclear) {

            switch(deriv_order_) {

              case 0:
              libint2_cleanup_elecpot(&primdata_[0]);
              break;
              case 1:
  #if LIBINT2_DERIV_ONEBODY_ORDER > 0
              libint2_cleanup_elecpot1(&primdata_[0]);
  #endif
              break;
              case 2:
  #if LIBINT2_DERIV_ONEBODY_ORDER > 1
              libint2_cleanup_elecpot2(&primdata_[0]);
  #endif
              break;
            }

            return;
          }
        }

      } // finalize()


  }; // struct OneBodyEngine
#endif // LIBINT2_SUPPORT_ONEBODY

  /// types of multiplicative spherically-symmetric two-body kernels known by TwoBodyEngine
  enum MultiplicativeSphericalTwoBodyKernel {
    Coulomb,            //!< \f$ 1/r_{12} \f$
    cGTG,               //!< contracted Gaussian geminal = \f$ \sum_i c_i \exp(- \alpha r_{12}^2) \f$
    cGTG_times_Coulomb, //!< contracted Gaussian geminal times Coulomb
    DelcGTG_square      //!< (\f$ \nabla \f$ cGTG) \f$ \cdot \f$ (\f$ \nabla \f$ cGTG)
  };

  /// contracted Gaussian geminal = \f$ \sum_i c_i \exp(- \alpha r_{12}^2) \f$, represented as a vector of
  /// {\f$ \alpha_i \f$, \f$ c_i \f$ } pairs
  typedef std::vector<std::pair<double,double>> ContractedGaussianGeminal;

  namespace detail {
    template <int K> struct R12_K_G12_to_Kernel;
    template <> struct R12_K_G12_to_Kernel<-1> {
        static const MultiplicativeSphericalTwoBodyKernel value = cGTG_times_Coulomb;
    };
    template <> struct R12_K_G12_to_Kernel<0> {
        static const MultiplicativeSphericalTwoBodyKernel value = cGTG;
    };
    template <> struct R12_K_G12_to_Kernel<2> {
        static const MultiplicativeSphericalTwoBodyKernel value = DelcGTG_square;
    };

    template <MultiplicativeSphericalTwoBodyKernel Kernel> struct TwoBodyEngineDispatcher;

  } // namespace detail

  template <MultiplicativeSphericalTwoBodyKernel Kernel> struct TwoBodyEngineTraits;
  template <> struct TwoBodyEngineTraits<Coulomb> {
      typedef libint2::FmEval_Chebyshev3<double> core_eval_type;
      //typedef libint2::FmEval_Taylor<double, 7> core_eval_type;
      typedef struct {} oper_params_type;
  };
  template <> struct TwoBodyEngineTraits<cGTG> {
      typedef libint2::GaussianGmEval<real_t, 0> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
  };
  template <> struct TwoBodyEngineTraits<cGTG_times_Coulomb> {
      typedef libint2::GaussianGmEval<real_t, -1> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
  };
  template <> struct TwoBodyEngineTraits<DelcGTG_square> {
      typedef libint2::GaussianGmEval<real_t, 2> core_eval_type;
      typedef ContractedGaussianGeminal oper_params_type;
  };


#ifdef LIBINT2_SUPPORT_ERI
  /**
   * TwoBodyEngine computes (ab|O|cd) (i.e. <em>four-center</em>) integrals over
   * a two-body kernel of type MultiplicativeSphericalTwoBodyKernel using Obara-Saika-Ahlrichs relations.
   *
   * \tparam Kernel kernel type, the supported values are enumerated by MultiplicativeSphericalTwoBodyKernel
   */
  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  class TwoBodyEngine {
    public:
      typedef typename libint2::TwoBodyEngineTraits<Kernel>::oper_params_type oper_params_type;

      /// creates a default (unusable) TwoBodyEngine
      TwoBodyEngine() : primdata_(), lmax_(-1), core_eval_(0) {
        set_precision(std::numeric_limits<real_t>::epsilon());
      }

      /// Constructs a (usable) TwoBodyEngine

      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_level if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_level
      /// \param precision specifies the target precision with which the integrals will be computed; the default is the "epsilon"
      ///        of \c real_t type, given by \c std::numeric_limits<real_t>::epsilon(). The precision
      ///        control is somewhat empirical, hence be conservative. \sa set_precision()
      /// \param oper_params specifies the operator parameters. The type of \c oper_params depends on \c Kernel as follows:
      ///        <ol>
      ///        <li> \c Coulomb : empty type (does not need to be provided) </li>
      ///        <li> \c cGTG : ContractedGaussianGeminal </li>
      ///        <li> \c cGTG_times_Coulomb : ContractedGaussianGeminal </li>
      ///        <li> \c DelcGTG_square : ContractedGaussianGeminal </li>
      ///        </ol>
      /// \warning currently only the following kernel types are suported: \c Coulomb
      /// \warning currently derivative integrals are not supported
      /// \warning currently only one-contraction Shell objects are supported; i.e. generally-contracted Shells are not yet supported
      TwoBodyEngine(size_t max_nprim, int max_l,
                    int deriv_order = 0,
                    real_t precision = std::numeric_limits<real_t>::epsilon(),
                    const oper_params_type& oparams = oper_params_type()) :
        primdata_(max_nprim * max_nprim * max_nprim * max_nprim),
        spbra_(max_nprim), spket_(max_nprim),
        lmax_(max_l), deriv_order_(deriv_order),
        core_eval_(core_eval_type::instance(4*lmax_ + deriv_order, std::min(std::numeric_limits<real_t>::epsilon(),precision)))
      {
        set_precision(precision);
        initialize();
        init_core_ints_params(oparams);
      }

      /// move constructor
      TwoBodyEngine(TwoBodyEngine&& other) = default;

      /// (deep) copy constructor
      TwoBodyEngine(const TwoBodyEngine& other) :
        primdata_(other.primdata_),
        spbra_(other.spbra_), spket_(other.spket_),
        lmax_(other.lmax_), deriv_order_(other.deriv_order_),
        precision_(other.precision_), ln_precision_(other.ln_precision_),
        core_eval_(other.core_eval_), core_ints_params_(other.core_ints_params_)
      {
        initialize();
      }

      ~TwoBodyEngine() {
        finalize();
      }

      /// move assignment
      TwoBodyEngine& operator=(TwoBodyEngine&& other) = default;

      /// (deep) copy assignment
      TwoBodyEngine& operator=(const TwoBodyEngine& other)
      {
        primdata_ = other.primdata_;
        lmax_ = other.lmax_;
        deriv_order_ = other.deriv_order_;
        precision_ = other.precision_;
        ln_precision_ = other.ln_precision_;
        core_eval_ = other.core_eval_;
        core_ints_params_ = other.core_ints_params_;
        initialize();
        return *this;
      }

#ifdef LIBINT2_ENGINE_TIMERS
      Timers<3> timers; // timers[0] -> prereqs
                        // timers[1] -> build (only meaningful if LIBINT2_PROFILE is not defined
                        // timers[2] -> tform
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
      struct class_id {
          size_t l[4];
          template <typename I> class_id(I l0, I l1, I l2, I l3) {
            l[0] = l0;
            l[1] = l1;
            l[2] = l2;
            l[3] = l3;
          }
          bool operator<(const class_id& other) const {
            return ordinal(l) < ordinal(other.l);
          }
          static size_t ordinal(const size_t (&l)[4]) {
            return ((l[0]*LIBINT2_MAX_AM_ERI + l[1])*LIBINT2_MAX_AM_ERI + l[2])*LIBINT2_MAX_AM_ERI + l[3];
          }
          std::string to_string() const {
            std::ostringstream oss;
            oss << "(" << Shell::am_symbol(l[0]) << Shell::am_symbol(l[1])
                << "|" << Shell::am_symbol(l[2]) << Shell::am_symbol(l[3]) << ")";
            return oss.str();
          }
      };
      struct class_profile {
          double prereqs;
          double build_hrr;
          double build_vrr;
          double tform;
          size_t nshellset; // total number of shell sets
          size_t nprimset;  // total number of primitive sets
          class_profile() { clear(); }
          class_profile(const class_profile& other) = default;
          void clear() {
            prereqs = build_hrr = build_vrr = tform = 0.;
            nprimset = nshellset = 0;
          }
      };
      std::map<class_id, class_profile> class_profiles;
#  endif
#endif

      /// computes shell set of integrals
      /// \note result is stored in the "chemists" form, i.e. (tbra1 tbra2 |tket1 tket2), in row-major order
      const real_t* compute(const libint2::Shell& tbra1,
                            const libint2::Shell& tbra2,
                            const libint2::Shell& tket1,
                            const libint2::Shell& tket2) {

        //
        // i.e. bra and ket refer to chemists bra and ket
        //

        // can only handle 1 contraction at a time
        assert(tbra1.ncontr() == 1 && tbra2.ncontr() == 1 &&
               tket1.ncontr() == 1 && tket2.ncontr() == 1);

        // derivatives not supported for now
        assert(deriv_order_ == 0);

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // standard angular momentum ordering
        auto swap_bra = (tbra1.contr[0].l < tbra2.contr[0].l);
        auto swap_ket = (tket1.contr[0].l < tket2.contr[0].l);
        auto swap_braket = (tbra1.contr[0].l + tbra2.contr[0].l > tket1.contr[0].l + tket2.contr[0].l);
#else // orca angular momentum ordering
        auto swap_bra = (tbra1.contr[0].l > tbra2.contr[0].l);
        auto swap_ket = (tket1.contr[0].l > tket2.contr[0].l);
        auto swap_braket = (tbra1.contr[0].l + tbra2.contr[0].l < tket1.contr[0].l + tket2.contr[0].l);
#endif
        const auto& bra1 = swap_braket ? (swap_ket ? tket2 : tket1) : (swap_bra ? tbra2 : tbra1);
        const auto& bra2 = swap_braket ? (swap_ket ? tket1 : tket2) : (swap_bra ? tbra1 : tbra2);
        const auto& ket1 = swap_braket ? (swap_bra ? tbra2 : tbra1) : (swap_ket ? tket2 : tket1);
        const auto& ket2 = swap_braket ? (swap_bra ? tbra1 : tbra2) : (swap_ket ? tket1 : tket2);

        const bool tform = tbra1.contr[0].pure || tbra2.contr[0].pure || tket1.contr[0].pure || tket2.contr[0].pure;
        const bool use_scratch = (swap_braket || swap_bra || swap_ket || tform);

        // assert # of primitive pairs
        auto nprim_bra1 = bra1.nprim();
        auto nprim_bra2 = bra2.nprim();
        auto nprim_ket1 = ket1.nprim();
        auto nprim_ket2 = ket2.nprim();

        // adjust max angular momentum, if needed
        auto lmax = std::max(std::max(bra1.contr[0].l, bra2.contr[0].l), std::max(ket1.contr[0].l, ket2.contr[0].l));
        assert (lmax <= lmax_);
        if (lmax == 0) // (ss|ss) ints will be accumulated in the first element of stack
          primdata_[0].stack[0] = 0;

#ifdef LIBINT2_ENGINE_PROFILE_CLASS
        class_id id(bra1.contr[0].l, bra2.contr[0].l, ket1.contr[0].l, ket2.contr[0].l);
        if (class_profiles.find(id) == class_profiles.end()) {
          class_profile dummy;
          class_profiles[id] = dummy;
        }
#endif

        // compute primitive data
#ifdef LIBINT2_ENGINE_TIMERS
        timers.start(0);
#endif
        {
          auto p = 0;
          spbra_.init(bra1, bra2, ln_precision_);
          spket_.init(ket1, ket2, ln_precision_);
          const auto npbra = spbra_.primpairs.size();
          const auto npket = spket_.primpairs.size();
          for(auto pb=0; pb!=npbra; ++pb) {
            for(auto pk=0; pk!=npket; ++pk) {
              if (spbra_.primpairs[pb].scr + spket_.primpairs[pk].scr > ln_precision_) {
                if (compute_primdata(primdata_[p],bra1,bra2,ket1,ket2,spbra_,pb,spket_,pk)) {
                  ++p;
                }
              }
            }
          }
          primdata_[0].contrdepth = p;
        }

#ifdef LIBINT2_ENGINE_TIMERS
        const auto t0 = timers.stop(0);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
        class_profiles[id].prereqs += t0.count();
        if (primdata_[0].contrdepth != 0) {
          class_profiles[id].nshellset += 1;
          class_profiles[id].nprimset += primdata_[0].contrdepth;
        }
#  endif
#endif

        // all primitive combinations screened out? return zeroes
        if (primdata_[0].contrdepth == 0) {
          const size_t n = bra1.size() * bra2.size() * ket1.size() * ket2.size();
          memset(primdata_[0].stack, 0, sizeof(real_t)*n);
          return primdata_[0].stack;
        }

        real_t* result = nullptr;

        if (lmax == 0) { // (ss|ss)
#ifdef LIBINT2_ENGINE_TIMERS
          timers.start(1);
#endif
          auto& stack = primdata_[0].stack[0];
          for(auto p=0; p != primdata_[0].contrdepth; ++p)
            stack += primdata_[p].LIBINT_T_SS_EREP_SS(0)[0];
          primdata_[0].targets[0] = primdata_[0].stack;
#ifdef LIBINT2_ENGINE_TIMERS
          const auto t1 = timers.stop(1);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
          class_profiles[id].build_vrr += t1.count();
#  endif
#endif

          result = primdata_[0].targets[0];
        }
        else { // not (ss|ss)
#ifdef LIBINT2_ENGINE_TIMERS
#    ifdef LIBINT2_PROFILE
          const auto t1_hrr_start = primdata_[0].timers->read(0);
          const auto t1_vrr_start = primdata_[0].timers->read(1);
#    endif
          timers.start(1);
#endif
          LIBINT2_PREFIXED_NAME(libint2_build_eri)[bra1.contr[0].l][bra2.contr[0].l][ket1.contr[0].l][ket2.contr[0].l](&primdata_[0]);
#ifdef LIBINT2_ENGINE_TIMERS
          const auto t1 = timers.stop(1);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
#    ifndef LIBINT2_PROFILE
          class_profiles[id].build_vrr += t1.count();
#    else
          class_profiles[id].build_hrr += primdata_[0].timers->read(0) - t1_hrr_start;
          class_profiles[id].build_vrr += primdata_[0].timers->read(1) - t1_vrr_start;
#    endif
#  endif
#endif

          result = primdata_[0].targets[0];

#ifdef LIBINT2_ENGINE_TIMERS
          timers.start(2);
#endif

          // if needed, permute and transform
          if (use_scratch) {

            constexpr auto using_scalar_real = std::is_same<double,real_t>::value || std::is_same<float,real_t>::value;
            static_assert(using_scalar_real, "Libint2 C++11 API only supports fundamental real types");
            typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;

            // a 2-d view of the 4-d source tensor
            const auto nr1_cart = bra1.cartesian_size();
            const auto nr2_cart = bra2.cartesian_size();
            const auto nc1_cart = ket1.cartesian_size();
            const auto nc2_cart = ket2.cartesian_size();
            const auto ncol_cart = nc1_cart * nc2_cart;
            const auto nr1 = bra1.size();
            const auto nr2 = bra2.size();
            const auto nc1 = ket1.size();
            const auto nc2 = ket2.size();
            const auto nrow = nr1 * nr2;
            const auto ncol = nc1 * nc2;

            // a 2-d view of the 4-d target tensor
            const auto nr1_tgt = tbra1.size();
            const auto nr2_tgt = tbra2.size();
            const auto nc1_tgt = tket1.size();
            const auto nc2_tgt = tket2.size();
            const auto ncol_tgt = nc1_tgt * nc2_tgt;

            // transform to solid harmonics first, then unpermute, if necessary
            auto mainbuf = result;
            auto scratchbuf = &scratch_[0];
            if (bra1.contr[0].pure) {
              libint2::solidharmonics::transform_first(bra1.contr[0].l, nr2_cart*ncol_cart,
                                                       mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }
            if (bra2.contr[0].pure) {
              libint2::solidharmonics::transform_inner(bra1.size(), bra2.contr[0].l, ncol_cart,
                                                       mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }
            if (ket1.contr[0].pure) {
              libint2::solidharmonics::transform_inner(nrow, ket1.contr[0].l, nc2_cart,
                                                       mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }
            if (ket2.contr[0].pure) {
              libint2::solidharmonics::transform_last(bra1.size()*bra2.size()*ket1.size(), ket2.contr[0].l,
                                                      mainbuf, scratchbuf);
              std::swap(mainbuf, scratchbuf);
            }

            // loop over rows of the source matrix
            const auto* src_row_ptr = mainbuf;
            auto tgt_ptr = scratchbuf;
            for(auto r1=0; r1!=nr1; ++r1) {
              for(auto r2=0; r2!=nr2; ++r2, src_row_ptr+=ncol) {

                typedef Eigen::Map<const Matrix> ConstMap;
                typedef Eigen::Map<Matrix> Map;
                typedef Eigen::Map<Matrix, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic> > StridedMap;

                // represent this source row as a matrix
                ConstMap src_blk_mat(src_row_ptr, nc1, nc2);

                // and copy to the block of the target matrix
                if (swap_braket) {
                  // if swapped bra and ket, a row of source becomes a column of target
                  // source row {r1,r2} is mapped to target column {r1,r2} if !swap_ket, else to {r2,r1}
                  const auto tgt_col_idx = !swap_ket ? r1 * nr2 + r2 : r2 * nr1 + r1;
                  StridedMap tgt_blk_mat(tgt_ptr + tgt_col_idx, nr1_tgt, nr2_tgt,
                                         Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(nr2_tgt*ncol_tgt,ncol_tgt));
                  if (swap_bra)
                    tgt_blk_mat = src_blk_mat.transpose();
                  else
                    tgt_blk_mat = src_blk_mat;
                }
                else {
                  // source row {r1,r2} is mapped to target row {r1,r2} if !swap_bra, else to {r2,r1}
                  const auto tgt_row_idx = !swap_bra ? r1 * nr2 + r2 : r2 * nr1 + r1;
                  Map tgt_blk_mat(tgt_ptr + tgt_row_idx*ncol, nc1_tgt, nc2_tgt);
                  if (swap_ket)
                    tgt_blk_mat = src_blk_mat.transpose();
                  else
                    tgt_blk_mat = src_blk_mat;
                }

              } // end of loop
            }   // over rows of source

            result = scratchbuf;

          } // if need_scratch => needed to transpose

#ifdef LIBINT2_ENGINE_TIMERS
          const auto t2 = timers.stop(2);
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
          class_profiles[id].tform += t2.count();
#  endif
#endif

        } // not (ss|ss)

        return result;
      }

      /** this specifies target precision for computing the integrals.
       *  target precision \f$ \epsilon \f$ is used in 3 ways:
       *  (1) to screen out primitive pairs in ShellPair object for which
       *      \f$ {\rm scr}_{12} = \max|c_1| \max|c_2| \exp(-\rho_{12} |AB|^2)/\gamma_{12} < \epsilon \f$ ;
       *  (2) to screen out primitive quartets outside compute_primdata() for which \f$ {\rm scr}_{12} {\rm scr}_{34} <  \epsilon \f$;
       *  (3) to screen out primitive quartets inside compute_primdata() for which the prefactor of \f$ F_m(\rho, T) \f$ is smaller
       *      than \f$ \epsilon \f$ .
       */
      void set_precision(real_t prec) {
        if (prec <= 0.) {
          precision_ = 0.;
          ln_precision_ = std::numeric_limits<real_t>::lowest();
        }
        else {
          precision_ = prec;
          ln_precision_ = std::log(precision_);
        }
      }
      /// @return the target precision for computing the integrals
      /// @sa set_precision(real_t)
      real_t precision() const {
        return precision_;
      }

      void print_timers() {
#ifdef LIBINT2_ENGINE_TIMERS
        std::cout << "timers: prereq = " << timers.read(0);
#  ifndef LIBINT2_PROFILE // if libint's profiling was on, engine's build timer will include its overhead
                          // do not report it, detailed profiling from libint will be printed below
        std::cout << " build = " << timers.read(1);
#  endif
        std::cout << " tform = " << timers.read(2) << std::endl;
#endif
#ifdef LIBINT2_PROFILE
        std::cout << "build timers: hrr = " << primdata_[0].timers->read(0)
                << " vrr = " << primdata_[0].timers->read(1) << std::endl;
#endif
#ifdef LIBINT2_ENGINE_TIMERS
#  ifdef LIBINT2_ENGINE_PROFILE_CLASS
        for(const auto& p: class_profiles) {
          printf("{\"%s\", %10.5lf, %10.5lf, %10.5lf, %10.5lf, %ld, %ld},\n",
                 p.first.to_string().c_str(),
                 p.second.prereqs,
                 p.second.build_vrr,
                 p.second.build_hrr,
                 p.second.tform,
                 p.second.nshellset,
                 p.second.nprimset);
        }
#  endif
#endif
      }

    private:

      inline bool compute_primdata(Libint_t& primdata,
                                   const Shell& sbra1,
                                   const Shell& sbra2,
                                   const Shell& sket1,
                                   const Shell& sket2,
                                   const ShellPair& spbra, size_t pbra,
                                   const ShellPair& spket, size_t pket);

      std::vector<Libint_t> primdata_;
      ShellPair spbra_, spket_;
      int lmax_;
      size_t deriv_order_;
      real_t precision_;
      real_t ln_precision_;

      typedef typename libint2::TwoBodyEngineTraits<Kernel>::core_eval_type core_eval_type;
      std::shared_ptr<core_eval_type> core_eval_;

      typedef oper_params_type core_ints_params_type; // currently core ints params are always same type as operator params
      core_ints_params_type core_ints_params_;
      /// converts operator parameters to core ints params
      void init_core_ints_params(const oper_params_type& oparams);

      std::vector<real_t> scratch_; // for transposes and/or transforming to solid harmonics

      friend struct detail::TwoBodyEngineDispatcher<Kernel>;

      void initialize() {
        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;
        const auto max_shellpair_size = ncart_max * ncart_max;
        const auto max_shellset_size = max_shellpair_size * max_shellpair_size;

        assert(lmax_ <= LIBINT2_MAX_AM_ERI);
        assert(deriv_order_ <= LIBINT2_DERIV_ONEBODY_ORDER);

        switch(deriv_order_) {

            case 0:
              libint2_init_eri(&primdata_[0], lmax_, 0);
              scratch_.resize(max_shellset_size);
              break;
            case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0
              libint2_init_eri1(&primdata_[0], lmax_, 0);
              scratch_.resize(9 * max_shellset_size);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1
              libint2_init_eri2(&primdata_[0], lmax_, 0);
              scratch_.resize(45 * max_shellset_size);
#endif
              break;
            default: assert(deriv_order_ < 3);
        }

#ifdef LIBINT2_ENGINE_TIMERS
        timers.set_now_overhead(25);
#endif
#ifdef LIBINT2_PROFILE
        primdata_[0].timers->set_now_overhead(25);
#endif
      }

      void finalize() {
        if (primdata_.size() != 0) {
          switch(deriv_order_) {

            case 0:
              libint2_cleanup_eri(&primdata_[0]);
              break;
            case 1:
#if LIBINT2_DERIV_ERI_ORDER > 0
              libint2_cleanup_eri1(&primdata_[0]);
#endif
              break;
            case 2:
#if LIBINT2_DERIV_ERI_ORDER > 1
              libint2_cleanup_eri2(&primdata_[0]);
#endif
              break;
          }
        }
      }

  }; // struct TwoBodyEngine

  namespace detail {
    template <> struct TwoBodyEngineDispatcher<Coulomb> {
        static void core_eval(TwoBodyEngine<Coulomb>* engine,
                              real_t* Fm,
                              int mmax,
                              real_t T,
                              real_t) {
          engine->core_eval_->eval(Fm, T, mmax);
        }
    };

    template <>
    struct TwoBodyEngineDispatcher<cGTG_times_Coulomb> {
        static void core_eval(TwoBodyEngine<cGTG_times_Coulomb>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_->eval(Gm, rho, T, mmax, engine->core_ints_params_);
        }
    };
    template <>
    struct TwoBodyEngineDispatcher<cGTG> {
        static void core_eval(TwoBodyEngine<cGTG>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_->eval(Gm, rho, T, mmax, engine->core_ints_params_);
        }
    };
    template <>
    struct TwoBodyEngineDispatcher<DelcGTG_square> {
        static void core_eval(TwoBodyEngine<DelcGTG_square>* engine,
                              real_t* Gm,
                              int mmax,
                              real_t T,
                              real_t rho) {
          engine->core_eval_->eval(Gm, rho, T, mmax, engine->core_ints_params_);
        }
    };
  }

  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  inline bool TwoBodyEngine<Kernel>::compute_primdata(Libint_t& primdata,
                                                      const Shell& sbra1,
                                                      const Shell& sbra2,
                                                      const Shell& sket1,
                                                      const Shell& sket2,
                                                      const ShellPair& spbra, size_t pbra,
                                                      const ShellPair& spket, size_t pket) {

    const auto& A = sbra1.O;
    const auto& B = sbra2.O;
    const auto& C = sket1.O;
    const auto& D = sket2.O;
    const auto& AB = spbra.AB;
    const auto& CD = spket.AB;

    const auto& spbrapp = spbra.primpairs[pbra];
    const auto& spketpp = spket.primpairs[pket];
    const auto& pbra1 = spbrapp.p1;
    const auto& pbra2 = spbrapp.p2;
    const auto& pket1 = spketpp.p1;
    const auto& pket2 = spketpp.p2;

    const auto alpha0 = sbra1.alpha[pbra1];
    const auto alpha1 = sbra2.alpha[pbra2];
    const auto alpha2 = sket1.alpha[pket1];
    const auto alpha3 = sket2.alpha[pket2];

    const auto c0 = sbra1.contr[0].coeff[pbra1];
    const auto c1 = sbra2.contr[0].coeff[pbra2];
    const auto c2 = sket1.contr[0].coeff[pket1];
    const auto c3 = sket2.contr[0].coeff[pket2];

    const auto amtot = sbra1.contr[0].l + sket1.contr[0].l +
                       sbra2.contr[0].l + sket2.contr[0].l;

    const auto gammap = alpha0 + alpha1;
    const auto oogammap = spbrapp.one_over_gamma;
    const auto rhop = alpha0 * alpha1 * oogammap;

    const auto gammaq = alpha2 + alpha3;
    const auto oogammaq = spketpp.one_over_gamma;
    const auto rhoq = alpha2 * alpha3 * oogammaq;

    const auto& P = spbrapp.P;
    const auto& Q = spketpp.P;
    const auto PQx = P[0] - Q[0];
    const auto PQy = P[1] - Q[1];
    const auto PQz = P[2] - Q[2];
    const auto PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

    const auto K12 = spbrapp.K * spketpp.K;
    decltype(K12) two_times_M_PI_to_25(34.986836655249725693); // (2 \pi)^{5/2}
    const auto gammapq = gammap + gammaq;
    const auto sqrt_gammapq = sqrt(gammapq);
    const auto oogammapq = 1.0 / (gammapq);
    auto pfac = two_times_M_PI_to_25 * K12 * sqrt_gammapq * oogammapq;
    pfac *= c0 * c1 * c2 * c3;

    if (std::abs(pfac) < precision_)
      return false;

    const auto rho = gammap * gammaq * oogammapq;
    const auto T = PQ2*rho;
    auto* fm_ptr = &(primdata.LIBINT_T_SS_EREP_SS(0)[0]);
    const auto mmax = amtot + deriv_order_;

    detail::TwoBodyEngineDispatcher<Kernel>::core_eval(this, fm_ptr, mmax, T, rho);

    for(auto m=0; m!=mmax+1; ++m) {
      fm_ptr[m] *= pfac;
    }

    if (mmax == 0)
      return true;

#if LIBINT2_DEFINED(eri,PA_x)
    primdata.PA_x[0] = P[0] - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
    primdata.PA_y[0] = P[1] - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
    primdata.PA_z[0] = P[2] - A[2];
#endif
#if LIBINT2_DEFINED(eri,PB_x)
    primdata.PB_x[0] = P[0] - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
    primdata.PB_y[0] = P[1] - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
    primdata.PB_z[0] = P[2] - B[2];
#endif

#if LIBINT2_DEFINED(eri,QC_x)
    primdata.QC_x[0] = Q[0] - C[0];
#endif
#if LIBINT2_DEFINED(eri,QC_y)
    primdata.QC_y[0] = Q[1] - C[1];
#endif
#if LIBINT2_DEFINED(eri,QC_z)
    primdata.QC_z[0] = Q[2] - C[2];
#endif
#if LIBINT2_DEFINED(eri,QD_x)
    primdata.QD_x[0] = Q[0] - D[0];
#endif
#if LIBINT2_DEFINED(eri,QD_y)
    primdata.QD_y[0] = Q[1] - D[1];
#endif
#if LIBINT2_DEFINED(eri,QD_z)
    primdata.QD_z[0] = Q[2] - D[2];
#endif

#if LIBINT2_DEFINED(eri,AB_x)
    primdata.AB_x[0] = AB[0];
#endif
#if LIBINT2_DEFINED(eri,AB_y)
    primdata.AB_y[0] = AB[1];
#endif
#if LIBINT2_DEFINED(eri,AB_z)
    primdata.AB_z[0] = AB[2];
#endif
#if LIBINT2_DEFINED(eri,BA_x)
    primdata.BA_x[0] = -AB[0];
#endif
#if LIBINT2_DEFINED(eri,BA_y)
    primdata.BA_y[0] = -AB[1];
#endif
#if LIBINT2_DEFINED(eri,BA_z)
    primdata.BA_z[0] = -AB[2];
#endif

#if LIBINT2_DEFINED(eri,CD_x)
    primdata.CD_x[0] = CD[0];
#endif
#if LIBINT2_DEFINED(eri,CD_y)
    primdata.CD_y[0] = CD[1];
#endif
#if LIBINT2_DEFINED(eri,CD_z)
    primdata.CD_z[0] = CD[2];
#endif
#if LIBINT2_DEFINED(eri,DC_x)
    primdata.DC_x[0] = -CD[0];
#endif
#if LIBINT2_DEFINED(eri,DC_y)
    primdata.DC_y[0] = -CD[1];
#endif
#if LIBINT2_DEFINED(eri,DC_z)
    primdata.DC_z[0] = -CD[2];
#endif

    const auto gammap_o_gammapgammaq = oogammapq * gammap;
    const auto gammaq_o_gammapgammaq = oogammapq * gammaq;

    const auto Wx = (gammap_o_gammapgammaq * P[0] + gammaq_o_gammapgammaq * Q[0]);
    const auto Wy = (gammap_o_gammapgammaq * P[1] + gammaq_o_gammapgammaq * Q[1]);
    const auto Wz = (gammap_o_gammapgammaq * P[2] + gammaq_o_gammapgammaq * Q[2]);

#if LIBINT2_DEFINED(eri,WP_x)
    primdata.WP_x[0] = Wx - P[0];
#endif
#if LIBINT2_DEFINED(eri,WP_y)
    primdata.WP_y[0] = Wy - P[1];
#endif
#if LIBINT2_DEFINED(eri,WP_z)
    primdata.WP_z[0] = Wz - P[2];
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
    primdata.WQ_x[0] = Wx - Q[0];
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
    primdata.WQ_y[0] = Wy - Q[1];
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
    primdata.WQ_z[0] = Wz - Q[2];
#endif
#if LIBINT2_DEFINED(eri,oo2z)
    primdata.oo2z[0] = 0.5*oogammap;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
    primdata.oo2e[0] = 0.5*oogammaq;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
    primdata.oo2ze[0] = 0.5*oogammapq;
#endif
#if LIBINT2_DEFINED(eri,roz)
    primdata.roz[0] = rho*oogammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
    primdata.roe[0] = rho*oogammaq;
#endif

    // using ITR?
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_x)
    primdata.TwoPRepITR_pfac0_0_0_x[0] = - (alpha1*AB[0] + alpha3*CD[0]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_y)
    primdata.TwoPRepITR_pfac0_0_0_y[0] = - (alpha1*AB[1] + alpha3*CD[1]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_0_z)
    primdata.TwoPRepITR_pfac0_0_0_z[0] = - (alpha1*AB[2] + alpha3*CD[2]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_x)
    primdata.TwoPRepITR_pfac0_1_0_x[0] = - (alpha1*AB[0] + alpha3*CD[0]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_y)
    primdata.TwoPRepITR_pfac0_1_0_y[0] = - (alpha1*AB[1] + alpha3*CD[1]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_0_z)
    primdata.TwoPRepITR_pfac0_1_0_z[0] = - (alpha1*AB[2] + alpha3*CD[2]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_x)
    primdata.TwoPRepITR_pfac0_0_1_x[0] = (alpha0*AB[0] + alpha2*CD[0]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_y)
    primdata.TwoPRepITR_pfac0_0_1_y[0] = (alpha0*AB[1] + alpha2*CD[1]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_0_1_z)
    primdata.TwoPRepITR_pfac0_0_1_z[0] = (alpha0*AB[2] + alpha2*CD[2]) * oogammap;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_x)
    primdata.TwoPRepITR_pfac0_1_1_x[0] = (alpha0*AB[0] + alpha2*CD[0]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_y)
    primdata.TwoPRepITR_pfac0_1_1_y[0] = (alpha0*AB[1] + alpha2*CD[1]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,TwoPRepITR_pfac0_1_1_z)
    primdata.TwoPRepITR_pfac0_1_1_z[0] = (alpha0*AB[2] + alpha2*CD[2]) * oogammaq;
#endif
#if LIBINT2_DEFINED(eri,eoz)
    primdata.eoz[0] = gammaq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,zoe)
    primdata.zoe[0] = gammap*oogammaq;
#endif

    // prefactors for derivative ERI relations
    if (deriv_order_ > 0) {
#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
      primdata.alpha1_rho_over_zeta2[0] = alpha0 * rho / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
      primdata.alpha2_rho_over_zeta2[0] = alpha1 * rho / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
      primdata.alpha3_rho_over_eta2[0] = alpha2 * rho / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
      primdata.alpha4_rho_over_eta2[0] = alpha3 * rho / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha1_over_zetapluseta)
      primdata.alpha1_over_zetapluseta[0] = alpha0 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha2_over_zetapluseta)
      primdata.alpha2_over_zetapluseta[0] = alpha1 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha3_over_zetapluseta)
      primdata.alpha3_over_zetapluseta[0] = alpha2 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_over_zetapluseta)
      primdata.alpha4_over_zetapluseta[0] = alpha3 / (gammap + gammaq);
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
      primdata.rho12_over_alpha1[0] = rhop / alpha0;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
      primdata.rho12_over_alpha2[0] = rhop / alpha1;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha3)
      primdata.rho34_over_alpha3[0] = rhoq / alpha2;
#endif
#if LIBINT2_DEFINED(eri,rho34_over_alpha4)
      primdata.rho34_over_alpha4[0] = rhoq / alpha3;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_bra)
      primdata.two_alpha0_bra[0] = 2.0 * alpha0;
#endif
#if LIBINT2_DEFINED(eri,two_alpha0_ket)
      primdata.two_alpha0_ket[0] = 2.0 * alpha1;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_bra)
      primdata.two_alpha1_bra[0] = 2.0 * alpha2;
#endif
#if LIBINT2_DEFINED(eri,two_alpha1_ket)
      primdata.two_alpha1_ket[0] = 2.0 * alpha3;
#endif
    }

    return true;
  }

  template <>
  inline void TwoBodyEngine<DelcGTG_square>::init_core_ints_params(
      const oper_params_type& oparams) {
    // [g12,[- \Del^2, g12] = 2 (\Del g12) \cdot (\Del g12)
    // (\Del exp(-a r_12^2) \cdot (\Del exp(-b r_12^2) = 4 a b (r_{12}^2 exp(- (a+b) r_{12}^2) )
    // i.e. need to scale each coefficient by 4 a b
    const auto ng = oparams.size();
    core_ints_params_.reserve(ng*(ng+1)/2);
    for(size_t b=0; b<ng; ++b)
      for(size_t k=0; k<=b; ++k) {
        const auto gexp = oparams[b].first + oparams[k].first;
        const auto gcoeff = oparams[b].second * oparams[k].second * (b == k ? 1 : 2); // if a != b include ab and ba
        const auto gcoeff_rescaled = 4 * oparams[b].first * oparams[k].first * gcoeff;
        core_ints_params_.push_back(std::make_pair(gexp, gcoeff_rescaled));
      }
  }

  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  inline void TwoBodyEngine<Kernel>::init_core_ints_params(
      const oper_params_type& oparams) {
    core_ints_params_ = oparams;
  }


#endif // LIBINT2_SUPPORT_ERI

} // namespace libint2

#endif /* _libint2_src_lib_libint_engine_h_ */
