#include "../tests/unit/catch.hpp"
#include <libint2/config.h>
#include <libint2/util/generated/libint2_params.h>

#ifdef __cplusplus
extern "C" /* prevent C++ name mangling */
#endif
void compute_eri_f(int *am1, double *alpha1, double *A,
                int *am2, double *alpha2, double *B,
                int *am3, double *alpha3, double *C,
                int *am4, double *alpha4, double *D,
                double *F, int *deriv_order);

#if LIBINT2_SUPPORT_ERI

TEST_CASE("Fortran ERI", "[eri]") {

  double A[] = {0.1, 1.3, -1.5};
  double B[] = {2.0, -0.5, 1.1};
  double C[] = {-1.2, 0.6, -0.1};
  double D[] = {0.4, 1.4, 0.3};

  int n1 = std::min(1, LIBINT2_MAX_AM_eri);
  int n2 = std::min(0, LIBINT2_MAX_AM_eri);
  int n3 = std::min(2, LIBINT2_MAX_AM_eri);
  int n4 = std::min(0, LIBINT2_MAX_AM_eri);

  double alpha1 = 0.1;
  double alpha2 = 1.0;
  double alpha3 = 0.5;
  double alpha4 = 1.9;

  double F[] = {0.41608906765397796,  0.044889937015574935,
                0.013706554295511562, 0.0063780699489852013,
                0.39523364424416996,  0.038762258204098128,
                0.010936175183774838, 0.0047907366138884629};

  int deriv_order = 0;

  compute_eri_f(&n1, &alpha1, A,
  &n2, &alpha2, B,
  &n3, &alpha3, C,
  &n4, &alpha4, D, F, &deriv_order);

  // test Fortran ints here
  // REQUIRE(true);
}

#endif // LIBINT2_SUPPORT_ERI