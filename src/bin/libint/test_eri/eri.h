
void calc_f(double *, int, double);       
double norm_const(unsigned int l1, unsigned int m1, unsigned int n1,       
                  double alpha1, double A[3]);

/*!
  eri()

  This is a very inefficient function for computing ERIs
  over primitive Gaussian functions. The argument
  list is self-explanatory, except for norm_flag:

  \param norm_flag:  tells what kind of normalization to use,
         0 - no normalization, >0 - normalized ERI
*/

double eri(unsigned int l1, unsigned int m1, unsigned int n1,
           double alpha1, double A[3],
           unsigned int l2, unsigned int m2, unsigned int n2,
           double alpha2, double B[3],
           unsigned int l3, unsigned int m3, unsigned int n3,
           double alpha3, double C[3],
           unsigned int l4, unsigned int m4, unsigned int n4,
           double alpha4, double D[3], int norm_flag);

double* init_array(unsigned long int size);
double** block_matrix(unsigned long int nrow, unsigned long int ncol);
void free_array(double* array);

