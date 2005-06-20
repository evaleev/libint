#include <vector>

void prep_libint2(Libint_t* libint, unsigned int am1,
           const std::vector<double>& alpha1, double A[3],
	   unsigned int am2, 
           const std::vector<double>& alpha2, double B[3],
	   unsigned int am3, 
           const std::vector<double>& alpha3, double C[3],
	   unsigned int am4, 
           const std::vector<double>& alpha4, double D[3], int norm_flag,
           unsigned int veclen);

