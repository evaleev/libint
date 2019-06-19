#include "catch.hpp"

#include <type_traits>

#include <libint2/config.h>
#include <libint2/boys.h>

#ifdef LIBINT_HAS_MPFR
TEST_CASE("Boys reference values", "[core-ints]") {
  using scalar_type = libint2::scalar_type;

  const int mmax = 40;
  std::vector<LIBINT2_REF_REALTYPE> T_ref_values;
  std::vector<double> T_values;
  T_ref_values.push_back(LIBINT2_REF_REALTYPE(0));
  T_values.push_back(0);
  for(int i=0; i!=12; ++i) {
    using std::sqrt;
    using std::pow;
    T_ref_values.push_back(pow(sqrt(LIBINT2_REF_REALTYPE(10)), i) / 1000);
    T_values.push_back(pow(sqrt(10.),i) / 1000.);
  }

  std::cout << std::setprecision(30);
  if (!std::is_same<LIBINT2_REF_REALTYPE,scalar_type>::value) {
    auto Fm_ref_values = new LIBINT2_REF_REALTYPE[mmax+1];
    auto Fm_values = new scalar_type[mmax+1];
    for(int t=0; t!=T_ref_values.size(); ++t) {
      auto Tref = T_ref_values[t];
      auto T = T_values[t];
      libint2::FmEval_Reference<LIBINT2_REF_REALTYPE>::eval(Fm_ref_values, Tref, mmax);
      auto fm_eval = libint2::FmEval_Chebyshev7<scalar_type>::instance(mmax);
//      auto fm_eval = libint2::FmEval_Taylor<scalar_type, 3>::instance(mmax);
      fm_eval->eval(Fm_values, T, mmax);
      for(int m=0; m<=mmax; ++m) {
        const LIBINT2_REF_REALTYPE value = libint2::sstream_convert<LIBINT2_REF_REALTYPE>(Fm_values[m]);
        const LIBINT2_REF_REALTYPE abs_error = abs(Fm_ref_values[m] - value);
        const LIBINT2_REF_REALTYPE relabs_error = abs(abs_error / Fm_ref_values[m]);
//        std::cout << "m=" << m << " T=" << Tref << " ref_value=" << Fm_ref_values[m] << " value=" << Fm_values[m];
//        printf(" abs_error=%e, relabs_error=%e\n", abs_error.get_d(), relabs_error.get_d());
        // can only get about 14 digits of precision for extreme cases, but can guarantee absolute precision to epsilon
        REQUIRE(abs_error.get_d() == Approx(0.0).margin(std::numeric_limits<scalar_type>::epsilon()) );
        REQUIRE(relabs_error.get_d() == Approx(0.0).margin(125*std::numeric_limits<scalar_type>::epsilon()) );
      }
    }
  }
}

#endif  // LIBINT_HAS_MPFR

TEST_CASE("Slater/Yukawa core integral values", "[core-ints]") {
  using scalar_type = libint2::scalar_type;

  const auto zeta = 1.0;
  int mmax = 10;
  const auto T = 0.3;
  const auto rho = 0.5;
  const auto one_over_rho = 1.0/rho;

  auto yukawa = new scalar_type[mmax + 1];
  auto stg = new scalar_type[mmax + 1];
  auto tenno_eval = libint2::TennoGmEval<scalar_type>::instance(mmax);
  // eval at T=0, U
  tenno_eval->eval_yukawa(yukawa, one_over_rho, 0, mmax, zeta);
  REQUIRE(yukawa[0] == Approx(0.34432045758120183));
  REQUIRE(yukawa[1] == Approx(0.21855984747293272));
  REQUIRE(yukawa[10] == Approx(0.04525107487757222));
  // eval at T=1023, U=10^-3 => barely inside the interpolation range
  tenno_eval->eval_yukawa(yukawa, 4e-3, 1023., mmax, zeta);
  REQUIRE(yukawa[0] == Approx(0.003668769827513457));
  REQUIRE(yukawa[1] == Approx(5.420435725593596e-6));
  REQUIRE(yukawa[10] == Approx(1.2693191281731064e-26));
  // eval at T=1025, U=10^-3 => outside the interpolation range, use upward recursion
  tenno_eval->eval_yukawa(yukawa, 4e-3, 1025., mmax, zeta);
  REQUIRE(yukawa[0] == Approx(0.0036579519797749903));
  REQUIRE(yukawa[1] == Approx(5.397434252949152e-6));
  REQUIRE(yukawa[10] == Approx(1.2432947194340530e-26));
  // eval at T, U
  tenno_eval->eval_yukawa(yukawa, one_over_rho, T, mmax, zeta);
  tenno_eval->eval_slater(stg, one_over_rho, T, mmax, zeta);
  REQUIRE(yukawa[0] == Approx(0.2852744664375206));
  REQUIRE(yukawa[1] == Approx(0.17668388024501083));
  REQUIRE(yukawa[10] == Approx(0.0343734522904477255));

  // compare to ints of STG approximated by Gaussians
  auto cgtg = new scalar_type[mmax + 1];
  auto cgtg_x_coulomb = new scalar_type[mmax + 1];
  // 6-term GTG fit of exp(-r12) from DOI 10.1063/1.1999632
  std::vector<std::pair<double, double>> cgtg_params{{0.2209,0.3144},{1.004,0.3037},{3.622,0.1681},{12.16,0.09811},{45.87,0.06024},{254.4,0.03726}};
  auto cgtg_eval = libint2::GaussianGmEval<scalar_type, 0>::instance(mmax);
  auto cgtg_x_coulomb_eval = libint2::GaussianGmEval<scalar_type, -1>::instance(mmax);
  cgtg_eval->eval(cgtg, rho, T, mmax, cgtg_params);
  cgtg_x_coulomb_eval->eval(cgtg_x_coulomb, rho, T, mmax, cgtg_params);

  for(int m=0; m<=mmax; ++m) {
    REQUIRE(yukawa[m] == Approx(cgtg_x_coulomb[m]).margin(1e-2));
    REQUIRE(stg[m] == Approx(cgtg[m]).margin(1e-2));
  }
}
