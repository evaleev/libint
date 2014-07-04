
#ifndef _libint2_src_lib_libint_compute_h_
#define _libint2_src_lib_libint_compute_h_

#if __cplusplus <= 199711L
# error "The simple Libint API requires C++11 support"
#endif

#include <iostream>
#include <array>
#include <vector>

#include <libint2.h>
//#include <libint2/shell.h>
#include <shell.h>

#include <Eigen/Core>

namespace libint2 {

  class OneBodyEngine {
    public:
      enum type {
        overlap,
        kinetic,
        coulomb,
        _invalid
      };

      OneBodyEngine() : type_(_invalid), primdata_(), lmax_(-1) {}
      OneBodyEngine(const OneBodyEngine& other) = default;
      OneBodyEngine(OneBodyEngine&& other) = default;
      OneBodyEngine(type t, size_t max_nprim, int max_l, int deriv_order = 0) :
        type_(t), primdata_(max_nprim * max_nprim), lmax_(max_l), deriv_order_(deriv_order) {

        const auto ncart_max = (lmax_+1)*(lmax_+2)/2;

        if (type_ == overlap) {
          assert(deriv_order_ <= LIBINT2_DERIV_ONEBODY_ORDER);
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

        assert(type_ == overlap);
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

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // make sure bra.l >= ket.l
        auto swap = (s1.contr[0].l < s2.contr[0].l);
#else // make sure bra.l <= ket.l
        auto swap = (s1.contr[0].l > s2.contr[0].l);
#endif
        auto bra = !swap ? s1 : s2;
        auto ket = !swap ? s2 : s1;

        // assert # of primitive pairs
        auto nprim_bra = bra.nprim();
        auto nprim_ket = ket.nprim();
        auto nprimpairs = nprim_bra * nprim_ket;
        assert(nprimpairs <= primdata_.size());

        // adjust max angular momentum, if needed
        auto lmax = std::max(bra.contr[0].l, ket.contr[0].l);
        assert (lmax <= lmax_);

        auto p12 = 0;
        for(auto pb=0; pb!=nprim_bra; ++pb) {
          for(auto pk=0; pk!=nprim_ket; ++pk, ++p12) {
            compute_primdata(primdata_[p12],bra,ket,pb,pk);
          }
        }
        primdata_[0].contrdepth = p12;

        if (lmax == 0) { // (s|s)
          auto& result = primdata_[0].stack[0];
          result = 0;
          for(auto p12=0; p12 != nprimpairs; ++p12)
            result += primdata_[p12]._aB_s___0___Overlap_s___0___Ab__up_[0];
          primdata_[0].targets[0] = primdata_[0].stack;
        }
        else {
          LIBINT2_PREFIXED_NAME(libint2_build_overlap)[bra.contr[0].l][ket.contr[0].l](&primdata_[0]);
        }

        if (swap) {
          const auto nbra = bra.size();
          const auto nket = ket.size();
          auto using_double = std::is_same<double,LIBINT2_REALTYPE>::value;
          assert(using_double);
          using namespace Eigen;
          Map<MatrixXd> set21(primdata_[0].targets[0], nbra, nket);
          Map<MatrixXd> set12(&scratch_[0], nket, nbra);
          set12 = set21.transpose();
          return &scratch_[0];
        }
        else
          return primdata_[0].targets[0];
      }

      void compute_primdata(Libint_t& primdata,
                            const Shell& s1, const Shell& s2,
                            size_t p1, size_t p2) {

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

        const auto PAx = Px - A[0];
        const auto PAy = Py - A[1];
        const auto PAz = Pz - A[2];
        const double PBx = Px - B[0];
        const double PBy = Py - B[1];
        const double PBz = Pz - B[2];

#if LIBINT2_SHELLQUARTET_SET == LIBINT2_SHELLQUARTET_SET_STANDARD // always VRR on bra, and HRR to bra

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

#else // always VRR on ket, HRR to ket

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

#endif

#if LIBINT2_DEFINED(eri,oo2z)
        primdata.oo2z[0] = 0.5*oogammap;
#endif

        if (deriv_order_ > 0) {
          // prefactors for derivative overlap relations
          assert(false);
        }

        const auto K1 = exp(- rhop * AB2) * oogammap;
        decltype(K1) sqrt_PI_cubed(5.56832799683170784528481798212);
        const auto ovlp_ss = sqrt_PI_cubed * sqrt(oogammap) * K1 * c1 * c2;

        primdata.LIBINT_T_S_OVERLAP_S[0] = ovlp_ss;

      }

    private:
      type type_;
      std::vector<Libint_t> primdata_;
      int lmax_;
      size_t deriv_order_;
      std::vector<LIBINT2_REALTYPE> scratch_; // for transposes and/or transforming to solid harmonics

  }; // struct Engine

} // namespace libint2

#endif /* _libint2_src_lib_libint_shell_h_ */
