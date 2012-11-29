/*
 * profile_Gm.cc
 *
 *      Author: David McDonald II && evaleev
 */

#include "boys.h"
#include <ctime>
#include <iostream>
#include <iterator>
#include <functional>
#include <boost/chrono.hpp>

using namespace std;

enum OperType {
  f12, f12_2, f12_o_r12, f12_t_f12
};
std::string to_string(OperType o) {
  std::string result;
  switch (o) {
    case f12: result = "f12"; break;
    case f12_2: result = "f12_2"; break;
    case f12_o_r12: result = "f12_o_r12"; break;
    case f12_t_f12: result = "f12_t_f12"; break;
  }
  return result;
}

void do_exp(double T, int nrepeats);
void do_chebyshev(int mmax, double T, int nrepeats);
void do_taylor(int mmax, double T, int nrepeats);
template <OperType O> void do_stg6g(int mmax, double T, double rho, int nrepeats);

int main(int argc, char* argv[]) {

  const int mmax  = atoi(argv[1]);
  const double T  = atol(argv[2]);
  const double rho = 1.0;
  const int nrepeats  = atoi(argv[3]);

  cout << "mmax = " << mmax << endl;
  cout << " T   = "<< T << endl;

  // changes precision of cout to 9 for printing doubles and such.
  cout.precision(9);

#ifndef SKIP_EXP
  do_exp(T, nrepeats);
#endif
#ifndef SKIP_CHEBYSHEV
  do_chebyshev(mmax, T, nrepeats);
#endif
#ifndef SKIP_TAYLOR
  do_taylor(mmax, T, nrepeats);
#endif
#ifndef SKIP_STG6G
  do_stg6g<f12>(mmax, T, rho, nrepeats);
  do_stg6g<f12_o_r12>(mmax, T, rho, nrepeats);
  do_stg6g<f12_2>(mmax, T, rho, nrepeats);
  do_stg6g<f12_t_f12>(mmax, T, rho, nrepeats);
#endif
}

void do_exp(double T, int nrepeats) {
  std::cout << "===================== exp(-T) ======================" << std::endl;
  double sum = 0.0;

  const auto start = boost::chrono::high_resolution_clock::now();
  for (int i = 0; i < nrepeats; ++i) {
    sum += exp(-T);
    T += 0.00001; // to ward-off unrealistic compiler optimizations
  }
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;

  std::cout << "sum of exp(-T): " << sum << std::endl;

  cout << "Time = " << fixed << elapsed.count() << endl;
  cout << "Rate = " << fixed << nrepeats / elapsed.count()  << endl;
}

void do_chebyshev(int mmax, double T, int nrepeats) {
  std::cout << "===================== Fm Cheb3 ======================" << std::endl;
  double* Fm_array = new double[mmax+1];
  double* Fm_array_sum = new double[mmax+1];
  //std::cout << "alignment of Fm = " << reinterpret_cast<unsigned long int>((char*) Fm_array) % 32ul << " bytes\n";

  libint2::FmEval_Chebyshev3 fmeval_cheb(mmax);
  std::cout << "done initialization:" << std::endl;
  std::fill(Fm_array_sum, Fm_array_sum+mmax+1, 0.0);
  const auto start = boost::chrono::high_resolution_clock::now();
  for (int i = 0; i < nrepeats; ++i) {
    // this computes all Fm for up to mmax
    fmeval_cheb.eval(Fm_array, T, mmax);
    for(int m=0; m<=mmax; ++m)
      Fm_array_sum[m] += Fm_array[m];
    T += 0.00001; // to ward-off unrealistic compiler optimizations
  }
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;

  std::cout << "sum of Fm:" << std::endl;
  std::copy(Fm_array_sum, Fm_array_sum+mmax+1, std::ostream_iterator<double>(std::cout,"\n"));

  cout << "Time = " << fixed << elapsed.count() << endl;
  cout << "Rate = " << fixed << nrepeats / elapsed.count()  << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

void do_taylor(int mmax, double T, int nrepeats) {
  std::cout << "===================== Fm Taylor6 ======================" << std::endl;
  double* Fm_array = new double[mmax+1];
  double* Fm_array_sum = new double[mmax+1];

  libint2::FmEval_Taylor<double, 6> fmeval_taylor6(mmax, 1e-15);
  std::fill(Fm_array_sum, Fm_array_sum+mmax+1, 0.0);
  const auto start = boost::chrono::high_resolution_clock::now();
  for (int i = 0; i < nrepeats; ++i)
  {
    // this computes all Fm for up to mmax
    fmeval_taylor6.eval(Fm_array, T, mmax);
    for(int m=0; m<=mmax; ++m)
      Fm_array_sum[m] += Fm_array[m];
    T += 0.00001; // to ward-off unrealistic compiler optimizations
  }
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;
  std::cout << "sum of Fm (Taylor):" << std::endl;
  std::copy(Fm_array_sum, Fm_array_sum+mmax+1, std::ostream_iterator<double>(std::cout,"\n"));

  cout << "Time = " << fixed << elapsed.count() << endl;
  cout << "Rate = " << fixed << nrepeats / elapsed.count()  << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

template<OperType O>
void do_stg6g(int mmax, double T, double rho, int nrepeats) {
  std::cout << "===================== Gm STG-6G ======================" << std::endl;
  double* Gm_array = new double[mmax+1];
  double* Gm_array_sum = new double[mmax+1];

  const size_t ng = 6;
  std::vector< std::pair<double,double> > stg_ng(ng);
#if HAVE_LAPACK
  libint2::stg_ng_fit(ng, 0.9, stg_ng);
#else
  stg_ng[0] = make_pair(0.16015391600067220691727771683865433704907890673261,
                        0.20306992259915090794062652264516576964313257462623);
  stg_ng[1] = make_pair(0.58691138376032812074703122125162923674902800850316,
                        0.29474840080158909154305767626309566520662550324323);
  stg_ng[2] = make_pair(1.9052880179050650706871766123674791603725239722860,
                        0.20652315861651088693388092033220845569017370830588);
  stg_ng[3] = make_pair(6.1508186412033182882412135545092215700186355770734,
                        0.13232619560602867340217449542493153700747744735317);
  stg_ng[4] = make_pair(22.558816746266648614747394893787336699960416307706,
                        0.084097701098685716800769730376212853566993914234229);
  stg_ng[5] = make_pair(167.12355778570626548864380047361110482628234458031,
                        0.079234606133959413896805606690618531996594605785539);
#endif

  std::vector< std::pair<double,double> > stg_ng_sq(stg_ng.size() * stg_ng.size());
  for(int b=0, bk=0; b<stg_ng.size(); ++b) {
    for(int k=0; k<stg_ng.size(); ++k, ++bk) {
      const double exp = stg_ng[b].first  + stg_ng[k].first;
      const double coef = stg_ng[b].second * stg_ng[k].second * (O == f12_t_f12 ? 4.0 * stg_ng[b].first  * stg_ng[k].first : 1.0);
      stg_ng_sq[bk] = make_pair(exp, coef);
    }
  }

  libint2::GaussianGmEval<double, (O == f12 || O == f12_2) ? 0 : (O == f12_o_r12 ? -1 : 2)>
    gtg_eval(mmax, 1e-15);
  std::fill(Gm_array_sum, Gm_array_sum+mmax+1, 0.0);
  const auto start = boost::chrono::high_resolution_clock::now();
  for (int i = 0; i < nrepeats; ++i)
  {
    // this computes all Gm for up to mmax
    gtg_eval.eval(Gm_array, rho, T, mmax,
                  (O == f12 || O == f12_o_r12) ? stg_ng : stg_ng_sq);
    for(int m=0; m<=mmax; ++m)
      Gm_array_sum[m] += Gm_array[m];
    T += 0.00001; // to ward-off unrealistic compiler optimizations
    rho += 0.000001;
  }
  const boost::chrono::duration<double> elapsed = boost::chrono::high_resolution_clock::now() - start;

  std::cout << "sum of Gm (STG-" << ng << "G O=" << to_string(O) << " ):" << std::endl;
  std::copy(Gm_array_sum, Gm_array_sum+mmax+1, std::ostream_iterator<double>(std::cout,"\n"));

  cout << "Time = " << fixed << elapsed.count() << endl;
  cout << "Rate = " << fixed << nrepeats / elapsed.count()  << endl;

  delete[] Gm_array;
  delete[] Gm_array_sum;
}

