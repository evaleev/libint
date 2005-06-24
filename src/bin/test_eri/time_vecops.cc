
#include <iostream>
#include <stdexcept>
#include <cstdlib>

namespace {
  typedef
    enum {  Test_daxpy = 0,
            Test_dzexy = 1 }
    Test;
};

void usage();
Test id_to_test(const std::string& id);

int main(int argc, char* argv[])
{
  if (argc != 4)
    usage();

  const long int ncycles = atol(argv[1]);
  const long int n = atol(argv[2]);
  std::string testid(argv[3]);
  Test test = id_to_test(testid);

  std::cout << "ncycles = " << ncycles << "  n = " << n 
            << "  testid = " << testid << std::endl;
  long int pfac = 1;
  if (test == Test_daxpy)
    pfac = 2;
  const long int nflops = pfac * ncycles * n;
  std::cout << "Performing " << nflops << " FLOPs" << std::endl;

  switch (test) {
    case Test_daxpy:
    {
      const double* X = new double[n];
      double* Y = new double[n];
      const double a = 10.0;
      for(long int c=0; c<ncycles; c++)
        for(long int i=0; i<n; i++)
          Y[i] += a*X[i];
    }
    break;

    case Test_dzexy:
    {
      const double* X = new double[n];
      const double* Y = new double[n];
      double* Z = new double[n];
      for(long int c=0; c<ncycles; c++)
        for(long int i=0; i<n; i++)
          Z[i] = X[i] * Y[i];
    }
    break;
  }

}

void usage()
{
  std::cerr << "Usage: time_vecops ncycles veclen testid" << std::endl
            << "         testid can be one of the following:" << std::endl
            << "         daxpy -- y[i] += a*x[i]" << std::endl
            << "         dzexy -- z[i] = x[i]*y[i]" << std::endl;
  exit(1);
}

Test
id_to_test(const std::string& id)
{
  if (id == "daxpy") {
    return Test_daxpy;
  }
  else if (id == "dzexy") {
    return Test_dzexy;
  }
  else
    throw std::runtime_error("Invalid test requested");
}

