/*
 *  Copyright (C) 2004-2021 Edward F. Valeev
 *
 *  This file is part of Libint.
 *
 *  Libint is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
#include "simd_wrapped_vector.hpp"
#include <iostream>
#include <cmath>
#if defined(_MSC_VER)
#  include <intrin.h>
#elif
#  include <x86intrin.h>
#endif
//#include <rr.h>
#include <iter.h>
#include <libint2/deriv_iter.h>
#include <policy_spec.h>
#include <global_macros.h>
#include <cgshell_ordering.h>

#ifndef _libint2_libint2types_h_
#define _libint2_libint2types_h_
//make sure the 4 lines below stay at line # 34-37 for timing script to work, or edit timing script
#define LIBINT2_MAX_VECLEN 1
constexpr auto am0 = 1u;
constexpr auto am1 = 1u;
constexpr auto am2 = 1u;
constexpr auto am3 = 1u;
constexpr auto am_tot = am0 + am1 + am2 + am3;
int bool_test = 0;
constexpr auto npts = am_tot / 2 + 1;

#define ALIGNED_STORAGE __attribute__((aligned(256)))

#include <libint2/util/vector.h>
#include <libint2/util/memory.h>
#include <libint2/util/intrinsic_operations.h>
#include "../../../../../libint_master/include/libint2/timer.h"

typedef struct {
    double AB_x[LIBINT2_MAX_VECLEN];
    double AB_y[LIBINT2_MAX_VECLEN];
    double AB_z[LIBINT2_MAX_VECLEN];
    double CD_x[LIBINT2_MAX_VECLEN];
    double CD_y[LIBINT2_MAX_VECLEN];
    double CD_z[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac0_0_x[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac0_0_y[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac0_0_z[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac0_1_x[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac0_1_y[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac0_1_z[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac1_0[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac1_1[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> R12kG12_pfac2[LIBINT2_MAX_VECLEN]; 
    ALIGNED_STORAGE VectorSIMD<double,npts> _00_GTG1d_00_x[LIBINT2_MAX_VECLEN]; 
    ALIGNED_STORAGE VectorSIMD<double,npts> _00_GTG1d_00_y[LIBINT2_MAX_VECLEN];
    ALIGNED_STORAGE VectorSIMD<double,npts> _00_GTG1d_00_z[LIBINT2_MAX_VECLEN];
    mutable VectorSIMD<double,npts>* stack;
    mutable VectorSIMD<double,npts>* vstack;
    mutable VectorSIMD<double,npts>* targets[LIBINT2_MAX_VECLEN];
    int veclen;
#if LIBINT2_FLOP_COUNT
    mutable LIBINT2_UINT_LEAST64* nflops;
#endif
    int contrdepth;
} Libint_t;
#endif

#include <VRR_GTG_1d_xx_xx_vec.h>
#include <libint/util.h>
#include <libint2/deriv_iter.h>
#include <test_eri/eri.h>
#ifdef LIBINT_HAVE_LIBROOTS
# include <roots/roots.hpp>
#endif

using namespace std;
using namespace libint2;

// these are used for casts
namespace libint2 {

  template <typename Output, typename Input>
  inline Output cast(Input i) {
    return Output(i);
  }

#ifdef LIBINT2_HAVE_AGNER_VECTORCLASS
  // AVX
  template <>
  inline double cast<double,Vec4d>(Vec4d i) {
    return i[0];
  }
  template <>
  inline double cast<double,Vec8f>(Vec8f i) {
    return i[0];
  }
  template <>
  inline float cast<float,Vec8f>(Vec8f i) {
    return i[0];
  }

  // SSE
  template <>
  inline double cast<double,Vec2d>(Vec2d i) {
    return i[0];
  }
  template <>
  inline double cast<double,Vec4f>(Vec4f i) {
    return i[0];
  }
  template <>
  inline float cast<float,Vec4f>(Vec4f i) {
    return i[0];
  }
#endif


  template <typename T>
  struct Tensor {
    public:
      Tensor() = default;
      Tensor(const Tensor&) = default;
 //     Tensor(Tensor&&) = default;
      ~Tensor() = default;

      template <class ... Dims> Tensor(Dims ... dims) : dims_{ std::forward<Dims>(dims)... } {
        strides_.resize(sizeof...(dims));
        // used in transform to compute strides
        struct stride {
            size_t value;
            stride() : value(1) {}
            size_t operator()(size_t dim) { auto old_value = value; value *= dim; return old_value; }
        };
        // row-major order of dimensions
        std::transform(dims_.rbegin(), dims_.rend(),
                       strides_.rbegin(), stride());
        size_t size = strides_.size() == 0 ? 0 : strides_[0] * dims_[0];
        data_.resize(size);
      }

      T* data() { return &data_[0]; }
      const T* data() const { return static_cast<const T*>(this->data()); }
      template <class ... MultiIndex> T* data(MultiIndex ... index) {
        assert(sizeof...(MultiIndex) <= dims_.size());
        auto index_list = { std::forward<MultiIndex>(index)... };
        size_t ordinal = std::inner_product(index_list.begin(), index_list.end(),
                                            strides_.begin(), 0);
        return &data_[ordinal];
      }
      template <class ... MultiIndex> const T* data(MultiIndex ... index) const {
          return static_cast<const T*>(this->data());
      }

      template <class ... MultiIndex> const T& operator()(MultiIndex ... index) {
        assert(sizeof...(MultiIndex) == dims_.size());
        auto index_list = { std::forward<MultiIndex>(index)... };
        size_t ordinal = std::inner_product(index_list.begin(), index_list.end(),
                                            strides_.begin(), 0);
        return data_[ordinal];
      }

    private:
      std::vector<size_t> dims_;
      std::vector<size_t> strides_;
      std::vector<T> data_;
  };
};

typedef unsigned int uint;

libint2::FmEval_Chebyshev7<double> fmeval_chebyshev(28);
libint2::FmEval_Taylor<double,7> fmeval_taylor(28, 1e-15);

int main(int argc, char** argv)
{
std::cout<<"number of points"<<npts<<endl;
#ifdef LIBINT_HAVE_LIBROOTS
  rysq::roots_initialize();
#endif

  libint2::Timers<3> timers; // 0 - roots/weights, 1 - 2d build, 2 - 6d build
  timers.set_now_overhead(25);

  //constexpr std::array<unsigned int, 4> am = {am0,am1,am2,am3};
	unsigned int am[] = {am0,am1,am2,am3};

  const uint veclen = 1;
  // for now hardwire contraction length to 1
  const uint contrdepth = 1;
  const uint contrdepth4 = contrdepth * contrdepth * contrdepth * contrdepth;

  const unsigned int deriv_order = 0;
  CartesianDerivIterator<4> diter(deriv_order);
  const unsigned int nderiv = diter.range_size();

  CGShell sh0(am[0]);
  CGShell sh1(am[1]);
  CGShell sh2(am[2]);
  CGShell sh3(am[3]);

  const std::array<double,3> A = {0.1, 0.4, 0.8};
  const std::array<double,3> B = {0.2, 0.5, 0.9};
  const std::array<double,3> C = { 0.3, 0.6, 1.0};
  const std::array<double,3> D = {0.4, 0.7, 1.1};
  const std::vector<double> alpha0 = {1};
  const std::vector<double> alpha1 = {2};
  const std::vector<double> alpha2 = {3};
  const std::vector<double> alpha3 = {4};
  const std::vector<double> c0{1};
  const std::vector<double> c1{1};
  const std::vector<double> c2{1};
  const std::vector<double> c3{1};
  // primitive indices are all 0 since hardwired contrdepth to 1
  const size_t p0 = 0;
  const size_t p1 = 0;
  const size_t p2 = 0;
  const size_t p3 = 0;


  LIBINT2_REF_REALTYPE Aref[4]; for(int i=0; i<4; ++i) Aref[i] = A[i];
  LIBINT2_REF_REALTYPE Bref[4]; for(int i=0; i<4; ++i) Bref[i] = B[i];
  LIBINT2_REF_REALTYPE Cref[4]; for(int i=0; i<4; ++i) Cref[i] = C[i];
  LIBINT2_REF_REALTYPE Dref[4]; for(int i=0; i<4; ++i) Dref[i] = D[i];

  typedef SubIteratorBase<CGShell> iter;
  SafePtr<iter> sh0_iter(new iter(sh0));
  SafePtr<iter> sh1_iter(new iter(sh1));
  SafePtr<iter> sh2_iter(new iter(sh2));
  SafePtr<iter> sh3_iter(new iter(sh3));

  //------------------------------------------------------
  // compute recurrence prefactors, Rys roots and weights
  //------------------------------------------------------
//  constexpr auto am_tot = am0 + am1 + am2 + am3;
//  constexpr auto npts = am_tot/2 + 1;

  std::vector<Libint_t> erieval(contrdepth4 + 1);
#if LIBINT2_FLOP_COUNT
  for(auto& v:erieval) { v.nflops = new LIBINT2_UINT_LEAST64; v.nflops[0] = 0; }
#endif

  /// prepare to compute 2-dimensional ints for quadrature point \c pt
  auto prep_data = [=,&erieval](Timers<3>& timers) {
    auto ev = erieval.begin();
    const auto a0 = alpha0[p0];
    const auto a1 = alpha1[p1];
    const auto a2 = alpha2[p2];
    const auto a3 = alpha3[p3];

    const auto gammap = a0 + a1;
    const auto oogammap = 1 / gammap;
    const auto rhop = a0 * a1 * oogammap;
    auto AB2 = 0.0; for(int xyz=0; xyz!=3; ++xyz) AB2 += (A[xyz]-B[xyz])*(A[xyz]-B[xyz]);

    const auto gammaq = a2 + a3;
    const auto oogammaq = 1 / gammaq;
    const auto rhoq = a2 * a3 * oogammaq;
    auto CD2 = 0.0; for(int xyz=0; xyz!=3; ++xyz) CD2 += (C[xyz]-D[xyz])*(C[xyz]-D[xyz]);

    LIBINT2_REALTYPE P[3]; for(auto xyz=0; xyz!=3; ++xyz) P[xyz]= (a0*A[xyz] + a1*B[xyz]) * oogammap;
    LIBINT2_REALTYPE Q[3]; for(auto xyz=0; xyz!=3; ++xyz) Q[xyz]= (a2*C[xyz] + a3*D[xyz]) * oogammaq;
    const auto PQx = P[0] - Q[0];
    const auto PQy = P[1] - Q[1];
    const auto PQz = P[2] - Q[2];
    const auto PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

    const auto gammapq = gammap + gammaq;
    const auto oogammapq = 1.0 / (gammapq);

    const auto rho = gammap * gammaq * oogammapq;
    const auto T = PQ2*rho;

    const auto k12 = exp(- rhop * AB2 - rhoq * CD2); // exp(-Gx-Gy-Gz) in RDK paper,
                                                     // DOI 10.1002/jcc.540040206
    const decltype(rho) oosqrtpi = 0.564189583547756286948079451561; // 1/sqrt(Pi)
    const auto pfac = 2 * sqrt(rho) * oosqrtpi * k12 * c0[p0] * c1[p1] * c2[p2] * c3[p3];

    // compute roots and weights for linear polynomial
    LIBINT2_REALTYPE gammas[npts]; // gamma = t^2
    LIBINT2_REALTYPE weights[npts];

    timers.start(0);
    
#ifdef LIBINT_HAVE_LIBROOTS
    {
      int32_t n = npts;
      rysq_roots(&n, const_cast<double*>(&T), &gammas[0], &weights[0]);
    }
    VectorSIMD<double,npts> gamma_vec(gammas);
    VectorSIMD<double,npts> weight_vec(weights);
#else
    assert(npts == 1);
    double Fm[2];
    fmeval_chebyshev.eval(Fm,T,1);
    gammas[0] = Fm[1]/Fm[0];
    weights[0] = Fm[0];
#endif

    timers.stop(0);

    
      std::cout  << ": " << " t=" << gamma_vec << " w=" << weight_vec << std::endl;

      VectorSIMD<double,npts> vals; 
      vals = 0.001;
      std::cout<<vals<<std::endl;
      vals = 0.0002;
      std::cout<<vals;
      std::cout<<ev->_00_GTG1d_00_x[0]<<std::endl;

      ev->_00_GTG1d_00_x[0] = M_PI * sqrt(oogammap * oogammaq);
      ev->_00_GTG1d_00_y[0] = M_PI * sqrt(oogammap * oogammaq);
      ev->_00_GTG1d_00_z[0] = M_PI * sqrt(oogammap * oogammaq);
      std::cout<<"made it even further"<<std::endl;
      // fold all extra factors into the x-axis integral
      ev->_00_GTG1d_00_x[0] = pfac * weight_vec * ev->_00_GTG1d_00_x[0];

      // compute VRR prefactors
      VectorSIMD<double,npts> gamma_o_gammapq = gamma_vec * oogammapq;
      VectorSIMD<double,npts> gammap_gamma_o_gammapq = gammap * gamma_o_gammapq;
      VectorSIMD<double,npts> gammaq_gamma_o_gammapq = gammaq * gamma_o_gammapq;

      // C_00 in RDK
      //
      VectorSIMD<double,npts> PA( P[0] - A[0]); 
      ev->R12kG12_pfac0_0_x[0] = PA - gammaq_gamma_o_gammapq * PQx;
      PA = P[1] - A[1];
      ev->R12kG12_pfac0_0_y[0] = PA - gammaq_gamma_o_gammapq * PQy;
      PA = P[2] - A[2];
      ev->R12kG12_pfac0_0_z[0] = PA - gammaq_gamma_o_gammapq * PQz;
      // C'_00 in RDK

      VectorSIMD<double,npts> QC(Q[0] - C[0]);
      ev->R12kG12_pfac0_1_x[0] = QC + gammap_gamma_o_gammapq * PQx;
      QC = Q[1] - C[1];
      ev->R12kG12_pfac0_1_y[0] = QC + gammap_gamma_o_gammapq * PQy;
      QC = Q[2] - C[2];
      ev->R12kG12_pfac0_1_z[0] = QC + gammap_gamma_o_gammapq * PQz;
      
      VectorSIMD<double,npts> ones(1);
      ev->R12kG12_pfac1_0[0] = 0.5 * oogammap * (ones - gammaq_gamma_o_gammapq); // B_10 in RDK
      ev->R12kG12_pfac1_0[0] = 0.5 * oogammap * (ones - gammaq_gamma_o_gammapq); // B_10 in RDK
      ev->R12kG12_pfac1_1[0] = 0.5 * oogammaq * (ones - gammap_gamma_o_gammapq); // B'_01 in RDK

      ev->R12kG12_pfac2[0] = 0.5 * gamma_o_gammapq; // B_00 in RDK
      std::cout <<"R12kG12: "<< ev->R12kG12_pfac2[0]<<std::endl;
      // horizontal RR prefactors
      ev->AB_x[0] = A[0]-B[0];
      ev->AB_y[0] = A[1]-B[1];
      ev->AB_z[0] = A[2]-B[2];
      ev->CD_x[0] = C[0]-D[0];
      ev->CD_y[0] = C[1]-D[1];
      ev->CD_z[0] = C[2]-D[2];
    

  };

  


  // prepare to compute 2-dimensional ints
  typedef LIBINT2_REALTYPE real_t;
  Tensor< VectorSIMD<double,npts> > gtg_x{1u,am[0]+1,am[1]+1,am[2]+1,am[3]+1};
  Tensor< VectorSIMD<double,npts> > gtg_y{1u,am[0]+1,am[1]+1,am[2]+1,am[3]+1};
  Tensor< VectorSIMD<double,npts> > gtg_z{1u,am[0]+1,am[1]+1,am[2]+1,am[3]+1};
  prep_data(timers);
  LIBINT2_UINT_LEAST64 nflops_build{0};
  timers.start(1);
  
  VRR_GTG_1d_xx_xx<CartesianAxis_X,am0,am1,am2,am3,false>::compute(&erieval[0],
                                                                     gtg_x.data(0),
                                                                     erieval[0]._00_GTG1d_00_x);
  VRR_GTG_1d_xx_xx<CartesianAxis_Y,am0,am1,am2,am3,false>::compute(&erieval[0],
                                                                     gtg_y.data(0),
                                                                     erieval[0]._00_GTG1d_00_y);
  VRR_GTG_1d_xx_xx<CartesianAxis_Z,am0,am1,am2,am3,false>::compute(&erieval[0],
                                                                     gtg_z.data(0),
                                                                     erieval[0]._00_GTG1d_00_z);

/*
  std::cout<<"hello"<<std::endl;


//#if LIBINT2_FLOP_COUNT
  //  nflops_build += erieval[pt].nflops[0];
//#endif
 // }
  timers.stop(1);

  cout << "Testing (" << sh0.label() << sh1.label()
  << "|" << sh2.label() << sh3.label() << ") ";
  if (deriv_order > 0) {
    cout << " deriv order = " << deriv_order;
  }
  cout << endl;


  bool success = true;
  int ijkl = 0;
  for(sh0_iter->init(); int(*sh0_iter); ++(*sh0_iter)) {
    for(sh1_iter->init(); int(*sh1_iter); ++(*sh1_iter)) {
      for(sh2_iter->init(); int(*sh2_iter); ++(*sh2_iter)) {
        for(sh3_iter->init(); int(*sh3_iter); ++(*sh3_iter), ijkl++) {

#if USE_BRAKET_H
          CGF bf0 = sh0_iter->elem();
          CGF bf1 = sh1_iter->elem();
          CGF bf2 = sh2_iter->elem();
          CGF bf3 = sh3_iter->elem();

          uint l0 = bf0.qn(0);
          uint m0 = bf0.qn(1);
          uint n0 = bf0.qn(2);
          uint l1 = bf1.qn(0);
          uint m1 = bf1.qn(1);
          uint n1 = bf1.qn(2);
          uint l2 = bf2.qn(0);
          uint m2 = bf2.qn(1);
          uint n2 = bf2.qn(2);
          uint l3 = bf3.qn(0);
          uint m3 = bf3.qn(1);
          uint n3 = bf3.qn(2);
#else
          SafePtr<CGF> bf0 = sh0_iter->elem();
          SafePtr<CGF> bf1 = sh1_iter->elem();
          SafePtr<CGF> bf2 = sh2_iter->elem();
          SafePtr<CGF> bf3 = sh3_iter->elem();

          uint l0 = bf0->qn(0);
          uint m0 = bf0->qn(1);
          uint n0 = bf0->qn(2);
          uint l1 = bf1->qn(0);
          uint m1 = bf1->qn(1);
          uint n1 = bf1->qn(2);
          uint l2 = bf2->qn(0);
          uint m2 = bf2->qn(1);
          uint n2 = bf2->qn(2);
          uint l3 = bf3->qn(0);
          uint m3 = bf3->qn(1);
          uint n3 = bf3->qn(2);
#endif

          for(uint v=0; v<veclen; v++) {

            std::vector<LIBINT2_REF_REALTYPE> ref_eri(nderiv, LIBINT2_REF_REALTYPE(0.0));

            uint p0123 = 0;
	    if(bool_test){

            for (uint p0 = 0; p0 < contrdepth; p0++) {
              for (uint p1 = 0; p1 < contrdepth; p1++) {
                for (uint p2 = 0; p2 < contrdepth; p2++) {
                  for (uint p3 = 0; p3 < contrdepth; p3++, p0123++) {

                    const LIBINT2_REF_REALTYPE a0 = alpha0[p0];
                    const LIBINT2_REF_REALTYPE a1 = alpha1[p1];
                    const LIBINT2_REF_REALTYPE a2 = alpha2[p2];
                    const LIBINT2_REF_REALTYPE a3 = alpha3[p3];

                    const LIBINT2_REF_REALTYPE c0123 = c0[p0] * c1[p1] * c2[p2] * c3[p3];

                    CartesianDerivIterator<4> diter(deriv_order);
                    bool last_deriv = false;
                    unsigned int di = 0;
                    do {
                      ref_eri[di++] += c0123 * eri(&(*diter)[0],
                                                   l0,m0,n0,a0,Aref,
                                                   l1,m1,n1,a1,Bref,
                                                   l2,m2,n2,a2,Cref,
                                                   l3,m3,n3,a3,Dref,
                                                   0);
                      last_deriv = diter.last();
                      if (!last_deriv) diter.next();
                    } while (!last_deriv);

                  }
                }
              }
           }

	    }
            for(unsigned int di = 0; di<nderiv; ++di) {

              // NOTE: quadrature weight was folded into gtg_x
	      LIBINT2_REALTYPE new_eri = 0;
             
              asm("#start of loop over points");
	      for(uint pt=0; pt!=npts; ++pt) { 							//I think I want to time this... ?
                new_eri += gtg_x(pt,l0,l1,l2,l3) * gtg_y(pt,m0,m1,m2,m3) * gtg_z(pt,n0,n1,n2,n3);
              }
	      asm("#end of loop over points");
              nflops_build += npts * 3;

              if ( abs(ref_eri[di] - libint2::cast<LIBINT2_REF_REALTYPE>(new_eri)) > 1.0E-10 && bool_test) {
                std::cout << "Elem " << ijkl << " di= " << di << " v=" << v
                    << " : eri.cc = " << ref_eri[di]
                    << " libint = " <<  libint2::cast<LIBINT2_REF_REALTYPE>(new_eri)
                    << " (relerr = " << abs((ref_eri[di] - libint2::cast<LIBINT2_REF_REALTYPE>(new_eri))/ref_eri[di]) << ")"
                    << endl;
                success = false;
              }
            }
	    

          } // end of vector loop
        }
      }
    }
  }


  if(bool_test)
  	cout << "test " << (success ? "ok" : "failed") << endl;
  else
	cout<< "\ntest not done. Change bool_test to 1 to complete test\n" <<endl;

  std::cout << "Rys build "<<am0<<" "<<am1<<" " <<am2<<" "<<am3<<" used " << nflops_build << " flops" << std::endl;

  std::cout << "Timers{Rys,2d-build} = {" << timers.read(0) << "," << timers.read(1) << "}\n";
#ifdef LIBINT_HAVE_LIBROOTS
  rysq::roots_finalize();
#endif
*/
  {
timers.start(2);


int int_num = 0;

VectorSIMD<double,npts> new_eri(0);
double new_eri_double = 0;
while(int_num < 1){
    int x0,y0,z0;
    FOR_CART(x0,y0,z0,am0)
      int x1,y1,z1;
      FOR_CART(x1,y1,z1,am1)
        int x2,y2,z2;
      	FOR_CART(x2,y2,z2,am2)
		int x3,y3,z3;
		FOR_CART(x3,y3,z3,am3)
		//	std::cout<<x0<<y1<<z1<<" "<<x1<<y1<<z1<<" "<<x2<<y2<<z2<<" "<<x3<<y3<<z3<<std::endl;
		    auto l0 = static_cast<unsigned int>(x0);
		    auto l1 = static_cast<unsigned int>(x1);
		    auto l2 = static_cast<unsigned int>(x2);
		    auto l3 = static_cast<unsigned int>(x3);	
				
       		    auto m0 = static_cast<unsigned int>(y0);
		    auto m1 = static_cast<unsigned int>(y1);
		    auto m2 = static_cast<unsigned int>(y2);
		    auto m3 = static_cast<unsigned int>(y3);
				
		    auto n0 = static_cast<unsigned int>(z0);
		    auto n1 = static_cast<unsigned int>(z1);
		    auto n2 = static_cast<unsigned int>(z2);
		    auto n3 = static_cast<unsigned int>(z3);
		
		    auto zero = static_cast<unsigned int>(0);
			    
		   new_eri = gtg_x(zero,l0,l1,l2,l3) * gtg_y(zero,m0,m1,m2,m3) * gtg_z(zero,n0,n1,n2,n3);
		   new_eri_double = reduce(new_eri); 
		   //Write generated function to do this...
		   //until the funciton is done for this
		   //the below has to be edited 
		   //for the number of points
		/*   if(npts <=4){
		       __m256d s = _mm256_hadd_pd(new_eri._avx0,new_eri._avx0);
		       new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];

		   }
		   
		   else if(npts >4 && npts <=8){
		        __m256d s = _mm256_hadd_pd(new_eri._avx0,new_eri._avx0);
		       new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		       s = _mm256_hadd_pd(new_eri._avx1,new_eri._avx1);
		       new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		   }*/
	/*	   
		   else if(npts > 8 && npts <= 12){		
		        __m256d s = _mm256_hadd_pd(new_eri._avx0,new_eri._avx0);
		       new_eri_double =  ((double*)&s)[0] + ((double*)&s)[2];
		       s = _mm256_hadd_pd(new_eri._avx1,new_eri._avx1);
		       new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		       s = _mm256_hadd_pd(new_eri._avx2,new_eri._avx2);
		       new_eri_double +=  ((double*)&s)[0] + ((double*)&s)[2];
		   }

	*/	
	           std::cout<<"new eri = "<<new_eri_double<<std::endl;
		END_FOR_CART
	END_FOR_CART
	END_FOR_CART
    END_FOR_CART
int_num++;	
  }
timers.stop(2);
std::cout<<new_eri<<endl;


//cout<<"6d build time per shell set: "<<timers.read(2)/max_num_ints<<std::endl;
  }

  return 0 ;
}


