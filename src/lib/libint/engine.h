
#ifndef _libint2_src_lib_libint_compute_h_
#define _libint2_src_lib_libint_compute_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>

#include <libint2.h>
//#include <libint2/boys.h>
#include <boys.h>
//#include <libint2/shell.h>
#include <shell.h>

#include <Eigen/Core>

namespace libint2 {

  class OneBodyEngine {
    public:
      enum type {
        overlap,
        kinetic,
        nuclear,
        _invalid
      };

      /// creates a default (unusable) OneBodyEngine
      OneBodyEngine() : type_(_invalid), primdata_(), lmax_(-1), fm_eval_(-1) {}

      OneBodyEngine(const OneBodyEngine& other) = default;
      OneBodyEngine(OneBodyEngine&& other) = default;

      /// Initializes a OneBodyEngine

      /// \param t integral type, see OneBodyEngine::type
      /// \param max_nprim the maximum number of primitives per contracted Gaussian shell
      /// \param max_l the maximum angular momentum of Gaussian shell
      /// \param deriv_level if not 0, will compute geometric derivatives of Gaussian integrals of order \c deriv_level
      /// \note if type == nuclear, must specify charges using set_q()
      OneBodyEngine(type t, size_t max_nprim, int max_l, int deriv_order = 0) :
        type_(t), primdata_(max_nprim * max_nprim), lmax_(max_l), deriv_order_(deriv_order),
        fm_eval_(t == nuclear ? static_cast<int>(2*lmax_ + deriv_order) : -1) {

        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;

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
        auto bra = !swap ? s1 : s2;
        auto ket = !swap ? s2 : s1;

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
                for(auto p12=0; p12 != nprimpairs; ++p12)
                  result += primdata_[p12].LIBINT_T_S_OVERLAP_S[0];
                  break;
                case kinetic:
                for(auto p12=0; p12 != nprimpairs; ++p12)
                  result += primdata_[p12].LIBINT_T_S_KINETIC_S[0];
                  break;
                case nuclear:
                for(auto p12=0; p12 != nprimpairs; ++p12)
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
                auto using_double = std::is_same<double,LIBINT2_REALTYPE>::value;
                assert(using_double);
                using namespace Eigen;
                Map<MatrixXd> braket(primdata_[0].targets[0], nbra, nket);
                Map<MatrixXd> set12(&scratch_[0], n1, n2);
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

        const auto A = s1.O;
        const auto B = s2.O;

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

  }; // struct Engine

} // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
