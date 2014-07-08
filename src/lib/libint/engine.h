
#ifndef _libint2_src_lib_libint_compute_h_
#define _libint2_src_lib_libint_compute_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>

#include <libint2.h>
#include <libint2/boys.h>
#include <libint2/shell.h>

#include <Eigen/Core>

namespace libint2 {

#ifdef LIBINT2_SUPPORT_ONEBODY
  /// OneBodyEngine computes integrals of 1-body operators, e.g. overlap, kinetic energy, dipole moment, etc.

  /**
   * OneBodyEngine computes integrals of types given by OneBodyEngine::type
   */
  class OneBodyEngine {
    public:
      enum type {
        overlap,
        kinetic,
        nuclear,
        _invalid
      };

      /// creates a default (unusable) OneBodyEngine; to be used as placeholder for copying a usable engine
      OneBodyEngine() : type_(_invalid), primdata_(), lmax_(-1), fm_eval_(-1) {}

      OneBodyEngine(const OneBodyEngine& other) = default;
      OneBodyEngine(OneBodyEngine&& other) = default;

      /// Constructs a (usable) OneBodyEngine

      /// \param t integral type, see OneBodyEngine::type
      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_level if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_level
      /// \note if type == nuclear, must specify charges using set_q()
      /// \warning currently only the following types are suported: \c overlap, \c kinetic, \c nuclear
      /// \warning currently derivative integrals are not supported
      /// \warning currently solid harmonics Gaussians are not supported
      OneBodyEngine(type t, size_t max_nprim, int max_l, int deriv_order = 0) :
        type_(t), primdata_(max_nprim * max_nprim), lmax_(max_l), deriv_order_(deriv_order),
        fm_eval_(t == nuclear ? static_cast<int>(2*lmax_ + deriv_order) : -1) {

        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;

        switch(type_) {
          case overlap: assert(max_l <= LIBINT2_MAX_AM_overlap); break;
          case kinetic: assert(max_l <= LIBINT2_MAX_AM_kinetic); break;
          case nuclear: assert(max_l <= LIBINT2_MAX_AM_elecpot); break;
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
      }

      /// specifies the nuclear charges
      /// \param q vector of {charge,Cartesian coordinate} pairs
      void set_q(const std::vector<std::pair<double, std::array<double, 3>>>& q) {
        q_ = q;
      }

      /// computes shell set of integrals
      /// \note result is stored in row-major order
      LIBINT2_REALTYPE* compute(const libint2::Shell& s1,
                                const libint2::Shell& s2) {

        // can only handle 1 contraction at a time
        assert(s1.ncontr() == 1 && s2.ncontr() == 1);

        // can only handle cartesian shells for now
        assert(s1.contr[0].pure == false && s2.contr[0].pure == false);

        // derivatives not supported for now
        assert(deriv_order_ == 0);

        // if want nuclear, make sure there is at least one nucleus .. otherwise the user likely forgot to call set_q
        if (type_ == nuclear and q_.size() == 0)
          throw std::runtime_error("libint2::OneBodyEngine(type = nuclear), but no nuclei found; forgot to call set_q()?");

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // make sure bra.l >= ket.l
        auto swap = (s1.contr[0].l < s2.contr[0].l);
#else // make sure bra.l <= ket.l
        auto swap = (s1.contr[0].l > s2.contr[0].l);
#endif
        const auto& bra = !swap ? s1 : s2;
        const auto& ket = !swap ? s2 : s1;

        const auto n1 = s1.size();
        const auto n2 = s2.size();

        const bool use_scratch = (swap || type_ == nuclear);

        // assert # of primitive pairs
        auto nprim_bra = bra.nprim();
        auto nprim_ket = ket.nprim();
        auto nprimpairs = nprim_bra * nprim_ket;
        assert(nprimpairs <= primdata_.size());

        // adjust max angular momentum, if needed
        auto lmax = std::max(bra.contr[0].l, ket.contr[0].l);
        assert (lmax <= lmax_);
        if (lmax == 0) // (s|s) ints will be accumulated in the first element of stack
          primdata_[0].stack[0] = 0;
        else if (use_scratch)
          memset(static_cast<void*>(&scratch_[0]), 0, sizeof(LIBINT2_REALTYPE)*s1.size()*s2.size());

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

            if (lmax == 0) { // (s|s)
              auto& result = primdata_[0].stack[0];
              switch (type_) {
                case overlap:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12].LIBINT_T_S_OVERLAP_S[0];
                  break;
                case kinetic:
                for(auto p12=0; p12 != primdata_[0].contrdepth; ++p12)
                  result += primdata_[p12].LIBINT_T_S_KINETIC_S[0];
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
                const auto nbra = bra.size();
                const auto nket = ket.size();
                constexpr auto using_scalar_real = std::is_same<double,LIBINT2_REALTYPE>::value || std::is_same<float,LIBINT2_REALTYPE>::value;
                static_assert(using_scalar_real, "Libint2 C++11 API only supports fundamental real types");
                typedef Eigen::Matrix<LIBINT2_REALTYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;
                Eigen::Map<Matrix> braket(primdata_[0].targets[0], nbra, nket);
                Eigen::Map<Matrix> set12(&scratch_[0], n1, n2);
                if (swap)
                  set12 += braket.transpose();
                else
                  set12 += braket;
              }
            } // ltot != 0

        } // oset (operator components, artifact of nuclear)

        if (use_scratch && lmax != 0)
          return &scratch_[0];
        else
          return primdata_[0].targets[0];
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
        const auto AB2 = AB_x*AB_x + AB_y*AB_y + AB_z*AB_z;

        if (LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // always VRR on bra, and HRR to bra (overlap, coulomb)
            || type_ == kinetic // kinetic energy ints don't use HRR, hence VRR on both centers
           ) {

#if LIBINT2_DEFINED(eri,PA_x)
        primdata.PA_x[0] = Px - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
        primdata.PA_y[0] = Py - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
        primdata.PA_z[0] = Pz - A[2];
#endif
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
        else { // always VRR on ket, HRR to ket (overlap, coulomb), or no HRR at all (kinetic energy ints)

#if LIBINT2_DEFINED(eri,PB_x)
        primdata.PB_x[0] = Px - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
        primdata.PB_y[0] = Py - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
        primdata.PB_z[0] = Pz - B[2];
#endif
#if LIBINT2_DEFINED(eri,BA_x)
        primdata.BA_x[0] = B[0] - A[0];
#endif
#if LIBINT2_DEFINED(eri,BA_y)
        primdata.BA_y[0] = B[1] - A[1];
#endif
#if LIBINT2_DEFINED(eri,BA_z)
        primdata.BA_z[0] = B[2] - A[2];
#endif
        }

#if LIBINT2_DEFINED(eri,oo2z)
        primdata.oo2z[0] = 0.5*oogammap;
#endif

        if (type_ == kinetic) { // additional factors for kinetic energy
#if LIBINT2_DEFINED(eri,rho12_over_alpha1)
        primdata.rho12_over_alpha1[0] = alpha2 * oogammap;
#endif
#if LIBINT2_DEFINED(eri,rho12_over_alpha2)
        primdata.rho12_over_alpha2[0] = alpha1 * oogammap;
#endif
#if LIBINT2_DEFINED(eri,two_rho12)
        primdata.two_rho12[0] = 2. * rhop;
#endif
        }

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
        }

        if (deriv_order_ > 0) {
          // prefactors for derivative overlap relations
          assert(false);
        }

        const auto K1 = exp(- rhop * AB2) * oogammap;
        decltype(K1) sqrt_PI_cubed(5.56832799683170784528481798212);
        const auto ovlp_ss = sqrt_PI_cubed * sqrt(oogammap) * K1 * c1 * c2;

        primdata.LIBINT_T_S_OVERLAP_S[0] = ovlp_ss;

        if (type_ == kinetic) {
          primdata.LIBINT_T_S_KINETIC_S[0] = rhop * (3. - 2.*rhop*AB2) * ovlp_ss;
        }

        if (type_ == nuclear) {
#if LIBINT2_DEFINED(eri,PC_x) && LIBINT2_DEFINED(eri,PC_y) && LIBINT2_DEFINED(eri,PC_z)
          const auto PC2 = primdata.PC_x[0] * primdata.PC_x[0] +
                           primdata.PC_y[0] * primdata.PC_y[0] +
                           primdata.PC_z[0] * primdata.PC_z[0];
          const auto U = gammap * PC2;
          const auto ltot = s1.contr[0].l + s2.contr[0].l;
          double* fm_ptr = &(primdata.LIBINT_T_S_ELECPOT_S(0)[0]);
          fm_eval_.eval(fm_ptr, U, ltot);

          decltype(U) two_o_sqrt_PI(1.12837916709551257389615890312);
          const auto pfac = - q_[oset].first * sqrt(gammap) * two_o_sqrt_PI * ovlp_ss;
          const auto ltot_p1 = ltot + 1;
          for(auto m=0; m!=ltot_p1; ++m) {
            fm_ptr[m] *= pfac;
          }
#endif
        }

      }

    private:
      type type_;
      std::vector<Libint_t> primdata_;
      int lmax_;
      size_t deriv_order_;
      std::vector<std::pair<double, std::array<double,3>>> q_;

      libint2::FmEval_Chebyshev3 fm_eval_;

      std::vector<LIBINT2_REALTYPE> scratch_; // for transposes and/or transforming to solid harmonics

  }; // struct OneBodyEngine
#endif // LIBINT2_SUPPORT_ONEBODY

  /// types of multiplicative spherically-symmetric two-body kernels known by TwoBodyEngine
  enum MultiplicativeSphericalTwoBodyKernel {
    Coulomb,            //!< \f$ 1/r_{12} \f$
    cGTG,               //!< contracted Gaussian geminal
    cGTG_times_Coulomb, //!< contracted Gaussian geminal times Coulomb
    DelcGTG_square      //!< \f$ \Nabla \f$ cGTG squared
  };

  namespace detail {
    template <MultiplicativeSphericalTwoBodyKernel Kernel> struct twobodyengine_traits;
    template <> struct twobodyengine_traits<Coulomb> {
        typedef libint2::FmEval_Chebyshev3 core_eval_type;
    };
    template <> struct twobodyengine_traits<cGTG> {
        typedef libint2::GaussianGmEval<LIBINT2_REALTYPE, 0> core_eval_type;
    };
    template <> struct twobodyengine_traits<cGTG_times_Coulomb> {
        typedef libint2::GaussianGmEval<LIBINT2_REALTYPE, -1> core_eval_type;
    };
    template <> struct twobodyengine_traits<DelcGTG_square> {
        typedef libint2::GaussianGmEval<LIBINT2_REALTYPE, 2> core_eval_type;
    };
  }

#ifdef LIBINT2_SUPPORT_ERI
  /// TwoBodyEngine computes (ab|O|cd) (i.e. <em>four-center</em>) integrals over
  /// a two-body kernel of type MultiplicativeSphericalTwoBodyKernel using Obara-Saika-Ahlrichs relations

  /**
   * \tparam Kernel kernel type
   */
  template <MultiplicativeSphericalTwoBodyKernel Kernel>
  class TwoBodyEngine {
    public:
      /// creates a default (unusable) TwoBodyEngine
      TwoBodyEngine() : primdata_(), lmax_(-1), core_eval_(-1) {}

      TwoBodyEngine(const TwoBodyEngine& other) = default;
      TwoBodyEngine(TwoBodyEngine&& other) = default;

      /// Initializes a TwoBodyEngine

      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_level if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_level
      /// \warning currently only the following kernel types are suported: \c Coulomb
      /// \warning currently derivative integrals are not supported
      /// \warning currently solid harmonics Gaussians are not supported
      TwoBodyEngine(size_t max_nprim, int max_l, int deriv_order = 0) :
        primdata_(max_nprim * max_nprim * max_nprim * max_nprim), lmax_(max_l), deriv_order_(deriv_order),
        core_eval_(4*lmax_ + deriv_order) {

        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;
        const auto max_shellset_size = ncart_max * ncart_max * ncart_max * ncart_max;

        assert(max_l <= LIBINT2_MAX_AM_ERI);
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

      }

      /// computes shell set of integrals
      /// \note result is stored in the "chemists" form, i.e. (tbra1 tbra2 |tket1 tket2), in row-major order
      LIBINT2_REALTYPE* compute(const libint2::Shell& tbra1,
                                const libint2::Shell& tbra2,
                                const libint2::Shell& tket1,
                                const libint2::Shell& tket2) {

        //
        // i.e. bra and ket refer to chemists bra and ket
        //

        // can only handle 1 contraction at a time
        assert(tbra1.ncontr() == 1 && tbra2.ncontr() == 1 &&
               tket1.ncontr() == 1 && tket2.ncontr() == 1);

        // can only handle cartesian shells for now
        assert(tbra1.contr[0].pure == false && tbra2.contr[0].pure == false &&
               tket1.contr[0].pure == false && tket2.contr[0].pure == false);

        // derivatives not supported for now
        assert(deriv_order_ == 0);

        // assert that user provided kernel params if needed ... for now just support Coulomb kernel
        assert(Kernel == Coulomb);

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

        const bool use_scratch = (swap_braket || swap_bra || swap_ket);

        // assert # of primitive pairs
        auto nprim_bra1 = bra1.nprim();
        auto nprim_bra2 = bra2.nprim();
        auto nprim_ket1 = ket1.nprim();
        auto nprim_ket2 = ket2.nprim();
        auto nprimsets = nprim_bra1 * nprim_bra2 * nprim_ket1 * nprim_bra2;
        assert(nprimsets <= primdata_.size());

        // adjust max angular momentum, if needed
        auto lmax = std::max(std::max(bra1.contr[0].l, bra2.contr[0].l), std::max(ket1.contr[0].l, ket2.contr[0].l));
        assert (lmax <= lmax_);
        if (lmax == 0) // (ss|ss) ints will be accumulated in the first element of stack
          primdata_[0].stack[0] = 0;

        // compute primitive data
        {
          auto p = 0;
          for(auto pb1=0; pb1!=nprim_bra1; ++pb1) {
            for(auto pb2=0; pb2!=nprim_bra2; ++pb2) {
              for(auto pk1=0; pk1!=nprim_ket1; ++pk1) {
                for(auto pk2=0; pk2!=nprim_ket2; ++pk2, ++p) {
                  compute_primdata(primdata_[p],bra1,bra2,ket1,ket2,pb1,pb2,pk1,pk2);
                }
              }
            }
          }
          primdata_[0].contrdepth = p;
        }

        if (lmax == 0) { // (ss|ss)
          auto& result = primdata_[0].stack[0];
          for(auto p=0; p != primdata_[0].contrdepth; ++p)
            result += primdata_[p].LIBINT_T_SS_EREP_SS(0)[0];
          primdata_[0].targets[0] = primdata_[0].stack;
        }
        else { // not (ss|ss)
          LIBINT2_PREFIXED_NAME(libint2_build_eri)[bra1.contr[0].l][bra2.contr[0].l][ket1.contr[0].l][ket2.contr[0].l](&primdata_[0]);

          // if needed, permute (and transform ... soon :)
          if (use_scratch) {

            constexpr auto using_scalar_real = std::is_same<double,LIBINT2_REALTYPE>::value || std::is_same<float,LIBINT2_REALTYPE>::value;
            static_assert(using_scalar_real, "Libint2 C++11 API only supports fundamental real types");
            typedef Eigen::Matrix<LIBINT2_REALTYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > Matrix;

            // a 2-d view of the 4-d source tensor
            const auto nr1 = bra1.size();
            const auto nr2 = bra2.size();
            const auto nc1 = ket1.size();
            const auto nc2 = ket2.size();
            const auto ncol = nc1 * nc2;

            // a 2-d view of the 4-d target tensor
            const auto nr1_tgt = tbra1.size();
            const auto nr2_tgt = tbra2.size();
            const auto nc1_tgt = tket1.size();
            const auto nc2_tgt = tket2.size();
            const auto ncol_tgt = nc1_tgt * nc2_tgt;

            auto tgt_ptr = &scratch_[0];

            // loop over rows of the source matrix
            const auto* src_row_ptr = primdata_[0].targets[0];
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
                  // source row {r1,r2} is mapped to target column {r1,r2} if !swap_bra, else to {r2,r1}
                  const auto tgt_col_idx = !swap_bra ? r1 * nr2 + r2 : r2 * nr1 + r1;
                  StridedMap tgt_blk_mat(tgt_ptr + tgt_col_idx, nr1_tgt, nr2_tgt, Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>(nr2_tgt*ncol_tgt,ncol_tgt));
                  if (swap_bra)
                    tgt_blk_mat = src_blk_mat.transpose();
                  else
                    tgt_blk_mat = src_blk_mat;
                }
                else {
                  const auto r12 = r1 * nr2 + r2;
                  Map tgt_blk_mat(tgt_ptr + r12*ncol, nc1_tgt, nc2_tgt);
                  if (swap_ket)
                    tgt_blk_mat = src_blk_mat.transpose();
                  else
                    tgt_blk_mat = src_blk_mat;
                }

              } // end of loop
            }   // over rows of source

          } // if need_scratch => needed to transpose

        } // not (ss|ss)

        if (use_scratch)
          return &scratch_[0];
        else
          return primdata_[0].targets[0];
      }

      void compute_primdata(Libint_t& primdata,
                            const Shell& sbra1,
                            const Shell& sbra2,
                            const Shell& sket1,
                            const Shell& sket2,
                            size_t pbra1,
                            size_t pbra2,
                            size_t pket1,
                            size_t pket2) {

        const auto& A = sbra1.O;
        const auto& B = sbra2.O;
        const auto& C = sket1.O;
        const auto& D = sket2.O;

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
        const auto oogammap = 1.0 / gammap;
        const auto rhop = alpha0 * alpha1 * oogammap;
        const auto Px = (alpha0 * A[0] + alpha1 * B[0]) * oogammap;
        const auto Py = (alpha0 * A[1] + alpha1 * B[1]) * oogammap;
        const auto Pz = (alpha0 * A[2] + alpha1 * B[2]) * oogammap;
        const auto AB_x = A[0] - B[0];
        const auto AB_y = A[1] - B[1];
        const auto AB_z = A[2] - B[2];
        const auto AB2 = AB_x*AB_x + AB_y*AB_y + AB_z*AB_z;

#if LIBINT2_DEFINED(eri,PA_x)
            primdata.PA_x[0] = Px - A[0];
#endif
#if LIBINT2_DEFINED(eri,PA_y)
            primdata.PA_y[0] = Py - A[1];
#endif
#if LIBINT2_DEFINED(eri,PA_z)
            primdata.PA_z[0] = Pz - A[2];
#endif
#if LIBINT2_DEFINED(eri,PB_x)
            primdata.PB_x[0] = Px - B[0];
#endif
#if LIBINT2_DEFINED(eri,PB_y)
            primdata.PB_y[0] = Py - B[1];
#endif
#if LIBINT2_DEFINED(eri,PB_z)
            primdata.PB_z[0] = Pz - B[2];
#endif

#if LIBINT2_DEFINED(eri,AB_x)
            primdata.AB_x[0] = AB_x;
#endif
#if LIBINT2_DEFINED(eri,AB_y)
            primdata.AB_y[0] = AB_y;
#endif
#if LIBINT2_DEFINED(eri,AB_z)
            primdata.AB_z[0] = AB_z;
#endif
#if LIBINT2_DEFINED(eri,BA_x)
            primdata.BA_x[0] = -AB_x;
#endif
#if LIBINT2_DEFINED(eri,BA_y)
            primdata.BA_y[0] = -AB_y;
#endif
#if LIBINT2_DEFINED(eri,BA_z)
            primdata.BA_z[0] = -AB_z;
#endif
#if LIBINT2_DEFINED(eri,oo2z)
            primdata.oo2z[0] = 0.5*oogammap;
#endif

            const auto gammaq = alpha2 + alpha3;
            const auto oogammaq = 1.0 / gammaq;
            const auto rhoq = alpha2 * alpha3 * oogammaq;
            const auto gammapq = gammap * gammaq / (gammap + gammaq);
            const auto gammap_o_gammapgammaq = gammapq * oogammaq;
            const auto gammaq_o_gammapgammaq = gammapq * oogammap;
            const auto Qx = (alpha2 * C[0] + alpha3 * D[0]) * oogammaq;
            const auto Qy = (alpha2 * C[1] + alpha3 * D[1]) * oogammaq;
            const auto Qz = (alpha2 * C[2] + alpha3 * D[2]) * oogammaq;
            const auto CD_x = C[0] - D[0];
            const auto CD_y = C[1] - D[1];
            const auto CD_z = C[2] - D[2];
            const auto CD2 = CD_x * CD_x + CD_y * CD_y + CD_z * CD_z;

#if LIBINT2_DEFINED(eri,QC_x)
            primdata.QC_x[0] = Qx - C[0];
#endif
#if LIBINT2_DEFINED(eri,QC_y)
            primdata.QC_y[0] = Qy - C[1];
#endif
#if LIBINT2_DEFINED(eri,QC_z)
            primdata.QC_z[0] = Qz - C[2];
#endif
#if LIBINT2_DEFINED(eri,QD_x)
            primdata.QD_x[0] = Qx - D[0];
#endif
#if LIBINT2_DEFINED(eri,QD_y)
            primdata.QD_y[0] = Qy - D[1];
#endif
#if LIBINT2_DEFINED(eri,QD_z)
            primdata.QD_z[0] = Qz - D[2];
#endif

#if LIBINT2_DEFINED(eri,CD_x)
            primdata.CD_x[0] = CD_x;
#endif
#if LIBINT2_DEFINED(eri,CD_y)
            primdata.CD_y[0] = CD_y;
#endif
#if LIBINT2_DEFINED(eri,CD_z)
            primdata.CD_z[0] = CD_z;
#endif
#if LIBINT2_DEFINED(eri,DC_x)
            primdata.DC_x[0] = -CD_x;
#endif
#if LIBINT2_DEFINED(eri,DC_y)
            primdata.DC_y[0] = -CD_y;
#endif
#if LIBINT2_DEFINED(eri,DC_z)
            primdata.DC_z[0] = -CD_z;
#endif
#if LIBINT2_DEFINED(eri,oo2e)
            primdata.oo2e[0] = 0.5*oogammaq;
#endif

            const auto PQx = Px - Qx;
            const auto PQy = Py - Qy;
            const auto PQz = Pz - Qz;
            const auto PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;
            const auto Wx = (gammap_o_gammapgammaq * Px + gammaq_o_gammapgammaq * Qx);
            const auto Wy = (gammap_o_gammapgammaq * Py + gammaq_o_gammapgammaq * Qy);
            const auto Wz = (gammap_o_gammapgammaq * Pz + gammaq_o_gammapgammaq * Qz);

#if LIBINT2_DEFINED(eri,WP_x)
            primdata.WP_x[0] = Wx - Px;
#endif
#if LIBINT2_DEFINED(eri,WP_y)
            primdata.WP_y[0] = Wy - Py;
#endif
#if LIBINT2_DEFINED(eri,WP_z)
            primdata.WP_z[0] = Wz - Pz;
#endif
#if LIBINT2_DEFINED(eri,WQ_x)
            primdata.WQ_x[0] = Wx - Qx;
#endif
#if LIBINT2_DEFINED(eri,WQ_y)
            primdata.WQ_y[0] = Wy - Qy;
#endif
#if LIBINT2_DEFINED(eri,WQ_z)
            primdata.WQ_z[0] = Wz - Qz;
#endif
#if LIBINT2_DEFINED(eri,oo2ze)
            primdata.oo2ze[0] = 0.5/(gammap+gammaq);
#endif
#if LIBINT2_DEFINED(eri,roz)
            primdata.roz[0] = gammapq*oogammap;
#endif
#if LIBINT2_DEFINED(eri,roe)
            primdata.roe[0] = gammapq*oogammaq;
#endif

            // prefactors for derivative ERI relations
            if (deriv_order_ > 0) {
#if LIBINT2_DEFINED(eri,alpha1_rho_over_zeta2)
            primdata.alpha1_rho_over_zeta2[0] = alpha0 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha2_rho_over_zeta2)
            primdata.alpha2_rho_over_zeta2[0] = alpha1 * gammapq / (gammap * gammap);
#endif
#if LIBINT2_DEFINED(eri,alpha3_rho_over_eta2)
            primdata.alpha3_rho_over_eta2[0] = alpha2 * gammapq / (gammaq * gammaq);
#endif
#if LIBINT2_DEFINED(eri,alpha4_rho_over_eta2)
            primdata.alpha4_rho_over_eta2[0] = alpha3 * gammapq / (gammaq * gammaq);
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

            const auto K1 = exp(- rhop * AB2);
            const auto K2 = exp(- rhoq * CD2);
            decltype(K1) two_times_M_PI_to_25(34.986836655249725693);
            double pfac = two_times_M_PI_to_25 * K1 * K2 / (gammap * gammaq * sqrt(gammap
                                                                                 + gammaq));
            pfac *= c0 * c1 * c2 * c3;

            const auto T = PQ2*gammapq;
            double* fm_ptr = &(primdata.LIBINT_T_SS_EREP_SS(0)[0]);
            const auto mmax = amtot + deriv_order_;
            core_eval_.eval(fm_ptr, T, mmax);

            for(auto m=0; m!=mmax+1; ++m) {
              fm_ptr[m] *= pfac;
            }

      }

    private:
      std::vector<Libint_t> primdata_;
      int lmax_;
      size_t deriv_order_;

      typename libint2::detail::twobodyengine_traits<Kernel>::core_eval_type core_eval_;

      std::vector<LIBINT2_REALTYPE> scratch_; // for transposes and/or transforming to solid harmonics

  }; // struct TwoBodyEngine
#endif // LIBINT2_SUPPORT_ERI

} // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
