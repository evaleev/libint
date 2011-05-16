
void calc_f(double *, int, double);       
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,       
                  double alpha1, const double* A);

/*!
  eri()

  This is a very inefficient function for computing ERIs
  over primitive Gaussian functions. The argument
  list is self-explanatory, except for norm_flag:

  \param norm_flag:  tells what kind of normalization to use,
         0 - no normalization, >0 - normalized ERI
*/

double eri(unsigned int l1, unsigned int m1, unsigned int n1,
           double alpha1, const double* A,
           unsigned int l2, unsigned int m2, unsigned int n2,
           double alpha2, const double* B,
           unsigned int l3, unsigned int m3, unsigned int n3,
           double alpha3, const double* C,
           unsigned int l4, unsigned int m4, unsigned int n4,
           double alpha4, const double* D, int norm_flag);

/// same as above, except specifies derivative order; uses the above function
/// to compute
/// \param deriv_order a set of 12 integers (3 coordinates for each of 4 basis function origins)
double eri(const unsigned int* deriv_index,
           unsigned int l1, unsigned int m1, unsigned int n1,
           double alpha1, const double* A,
           unsigned int l2, unsigned int m2, unsigned int n2,
           double alpha2, const double* B,
           unsigned int l3, unsigned int m3, unsigned int n3,
           double alpha3, const double* C,
           unsigned int l4, unsigned int m4, unsigned int n4,
           double alpha4, const double* D, int norm_flag);

double* init_array(unsigned long int size);
double** block_matrix(unsigned long int nrow, unsigned long int ncol);
void free_array(double* array);

