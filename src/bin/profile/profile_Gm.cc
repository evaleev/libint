/*
 * profile_Gm.cc
 *
 *      Author: David McDonald II && evaleev
 */

#include "boys.h"
#include <ctime>
#include <iostream>
#include <iterator>

using namespace std;

void do_chebyshev(int mmax, double T, int nrepeats);
void do_taylor(int mmax, double T, int nrepeats);
void do_stg6g(int mmax, double T, int nrepeats);

int main(int argc, char* argv[]) {

  const int mmax  = atoi(argv[1]);
  const double T  = atol(argv[2]);
  const int nrepeats  = atoi(argv[3]);

  cout << "mmax = " << mmax << endl;
  cout << " T   = "<< T << endl;

  // changes precision of cout to 9 for printing doubles and such.
  cout.precision(9);

#ifndef SKIP_CHEBYSHEV
  do_chebyshev(mmax, T, nrepeats);
#endif
#ifndef SKIP_TAYLOR
  do_taylor(mmax, T, nrepeats);
#endif
#ifndef SKIP_STG6G
  do_stg6g(mmax, T, nrepeats);
#endif
}

void do_chebyshev(int mmax, double T, int nrepeats) {
  double* Fm_array = new double[mmax+1];
  double* Fm_array_sum = new double[mmax+1];
  //std::cout << "alignment of Fm = " << reinterpret_cast<unsigned long int>((char*) Fm_array) % 32ul << " bytes\n";

  libint2::FmEval_Chebyshev3 fmeval_cheb(mmax);
  std::fill(Fm_array_sum, Fm_array_sum+mmax+1, 0.0);
  const clock_t startChebyshev=clock();
  for (int i = 0; i < nrepeats; ++i) {
    // this computes all Fm for up to mmax
    fmeval_cheb.eval(Fm_array, T, mmax);
    for(int m=0; m<=mmax; ++m)
      Fm_array_sum[m] += Fm_array[m];
  }
  const clock_t stopChebyshev=clock();

  std::cout << "sum of Fm (Chebyshev):" << std::endl;
  std::copy(Fm_array_sum, Fm_array_sum+mmax+1, std::ostream_iterator<double>(std::cout,"\n"));

  const double timetotalCheb = ((double)(stopChebyshev - startChebyshev)) / (double)CLOCKS_PER_SEC;
  cout << "Time for Cheb" << endl;
  cout << fixed << timetotalCheb << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

void do_taylor(int mmax, double T, int nrepeats) {
  double* Fm_array = new double[mmax+1];
  double* Fm_array_sum = new double[mmax+1];

  libint2::FmEval_Taylor<double, 6> fmeval_taylor6(mmax, 1e-15);
  std::fill(Fm_array_sum, Fm_array_sum+mmax+1, 0.0);
  const clock_t startTaylor=clock();
  for (int i = 0; i < nrepeats; ++i)
  {
    // this computes all Fm for up to mmax
    fmeval_taylor6.eval(Fm_array, T, mmax);
    for(int m=0; m<=mmax; ++m)
      Fm_array_sum[m] += Fm_array[m];
  }
  const clock_t stopTaylor=clock();
  std::cout << "sum of Fm (Taylor):" << std::endl;
  std::copy(Fm_array_sum, Fm_array_sum+mmax+1, std::ostream_iterator<double>(std::cout,"\n"));

  const double timetotalTaylor = ((double)(stopTaylor - startTaylor)) / (double)CLOCKS_PER_SEC;
  cout << "Time for Taylor" << endl;
  cout << fixed << timetotalTaylor << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

void do_stg6g(int mmax, double T, int nrepeats) {
  double* Fm_array = new double[mmax+1];
  double* Fm_array_sum = new double[mmax+1];

  const size_t ng = 6;
  std::vector< std::pair<double,double> > stg_ng(ng);
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
  std::vector< std::pair<double,double> > stg_ng_sq(stg_ng.size() * stg_ng.size());
  std::vector< std::pair<double,double> > stg_ng_fTf(stg_ng_sq.size());
  const double rho = 10.0; // reduced exponent of primitive pairs
  libint2::GaussianGmEval<double, -1> gtg_eval_m1(mmax, 1e-15);

  const clock_t startGTGm1=clock();
  for (int i = 0; i < nrepeats; ++i)
  {
    // this computes all Fm for up to mmax
    gtg_eval_m1.eval(Fm_array, rho, T, mmax,
                     stg_ng);
  }
  const clock_t stopGTGm1=clock();

  const double timetotalGTGm1 = ((double)(stopGTGm1 - startGTGm1)) / (double)CLOCKS_PER_SEC;
  cout << "Time for GTGm1" << endl;
  cout << fixed << timetotalGTGm1 << endl;

  delete[] Fm_array;
  delete[] Fm_array_sum;
}

