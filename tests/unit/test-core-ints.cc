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
  {
    const std::vector<scalar_type> ref_values{
        0.344320457581201528456128769269,   0.2185598474729328238479570769103,
        0.1562880305054134352304085846179,  0.120530281356369509252798773626,
        0.0977188576270700545274668029304,  0.08202555839753908595204847246087,
        0.07061341858480468569599627134916, 0.06195910542767968762026691524339,
        0.05517887615131295955174900498568, 0.04972742757098352844464478921128,
        0.04525107487757221293120739098994};
    tenno_eval->eval_yukawa(yukawa, one_over_rho, 0, mmax, zeta);
    for (int m = 0; m <= mmax; ++m) {
      REQUIRE(std::abs(yukawa[m] - ref_values[m]) <= 1.4e-14);
    }
  }
  // eval at T=1023, U=10^-3 => barely inside the interpolation range
  {
    const std::vector<scalar_type> ref_values{
        0.0036687698275134500493040872,   5.420435725593595129051675e-6,
        1.153413823646514442119413472e-8, 3.34856122353435544301437e-11,
        1.25839473177094413418081e-13,    5.8627882847729072806599e-16,
        3.2750469499532682480753e-18,     2.13822913031999358167698e-20,
        1.59963080864078970551172e-22,    1.35001806319439998582769e-24,
        1.26931912817310643557324e-26};
    tenno_eval->eval_yukawa(yukawa, 4e-3, 1023., mmax, zeta);
    for (int m = 0; m <= mmax; ++m) {
      REQUIRE(std::abs(yukawa[m] - ref_values[m]) <= 1.4e-14);
    }
  }
  // eval at T=1025, U=10^-3 => outside the interpolation range, use upward recursion
  {
    const std::vector<scalar_type> ref_values{
        0.0036579519797749897405471452,   5.397434252949152394251251e-6,
        1.146741791141338373846245954e-8, 3.32351014941293773077145e-11,
        1.24673437210601174941915e-13,    5.7977128677252162394764e-16,
        3.23260050191167815283312e-18,    2.10650483406813947486468e-20,
        1.57288256640997208553838e-22,    1.32489290711137333410368e-24,
        1.2432947194340530617111e-26};
    tenno_eval->eval_yukawa(yukawa, 4e-3, 1025., mmax, zeta);
    for (int m = 0; m <= mmax; ++m) {
      REQUIRE(std::abs(yukawa[m] - ref_values[m]) <= 1.4e-14);
    }
  }
  // eval at T, U
  {
    const std::vector<scalar_type> ref_values{0.2852744664375203210290618028,
                                              0.176683880245011075178211899,
                                              0.1241798108180594674947062,
                                              0.0946078560893175776414,
                                              0.0760276379359410749,
                                              0.063397294718448975574,
                                              0.05429943192860323351,
                                              0.047452815180955241894,
                                              0.04212239826868999308,
                                              0.037858941778278764,
                                              0.03437345229044773};
    tenno_eval->eval_yukawa(yukawa, one_over_rho, T, mmax, zeta);
    for (int m = 0; m <= mmax; ++m) {
      REQUIRE(std::abs(yukawa[m] - ref_values[m]) <= 1.4e-14);
    }
  }

  // compare to ints of STG approximated by Gaussians
  {
    tenno_eval->eval_yukawa(yukawa, one_over_rho, T, mmax, zeta);
    tenno_eval->eval_slater(stg, one_over_rho, T, mmax, zeta);
    auto cgtg = new scalar_type[mmax + 1];
    auto cgtg_x_coulomb = new scalar_type[mmax + 1];
    // precise fit of exp(-r12) on [0,10)
    std::vector<std::pair<double, double>> cgtg_params = {
        {0.10535330565471572, 0.08616353459042002},
        {0.22084823136587992, 0.08653979627551414},
        {0.3543431104992702, 0.08803697599356214},
        {0.48305514665749105, 0.09192519612306953},
        {0.6550035700167584, 0.10079776426873248},
        {1.1960917050801643, 0.11666110644901695},
        {2.269278814810891, 0.14081371547404428},
        {5.953990617813977, 0.13410255216448014},
        {18.31911063199608, 0.0772095196191394},
        {66.98443868169818, 0.049343985939540556},
        {367.24137290439205, 0.03090625839896873},
        {5.655142311118115, -0.017659052507938647}};
    auto cgtg_eval = libint2::GaussianGmEval<scalar_type, 0>::instance(mmax);
    auto cgtg_x_coulomb_eval =
        libint2::GaussianGmEval<scalar_type, -1>::instance(mmax);
    cgtg_eval->eval(cgtg, rho, T, mmax, cgtg_params);
    cgtg_x_coulomb_eval->eval(cgtg_x_coulomb, rho, T, mmax, cgtg_params);

    for (int m = 0; m <= mmax; ++m) {
      REQUIRE(std::abs(yukawa[m] - cgtg_x_coulomb[m]) <= 5e-4);
      REQUIRE(std::abs(stg[m] - cgtg[m]) <= 2e-4);
    }
  }
}
