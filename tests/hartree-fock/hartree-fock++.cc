/*
 *  This file is a part of Libint.
 *  Copyright (C) 2004-2014 Edward F. Valeev
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 */

// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>
#include <atomic>
#include <iterator>
#include <unordered_map>
#include <mutex>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>

// have BTAS library?
#ifdef LIBINT2_HAVE_BTAS
# include <btas/btas.h>
#endif // LIBINT2_HAVE_BTAS

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/diis.h>

#if defined(_OPENMP)
# include <omp.h>
#endif

// uncomment if want to report integral timings (only useful if nthreads == 1)
// N.B. integral engine timings are controled in engine.h
//#define REPORT_INTEGRAL_TIMINGS

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic>
        DiagonalMatrix;

using libint2::Shell;
using libint2::Atom;
using libint2::BasisSet;

std::vector<Atom> read_geometry(const std::string& filename);
Matrix compute_soad(const std::vector<Atom>& atoms);
// computes norm of shell-blocks of A
Matrix compute_shellblock_norm(const BasisSet& obs,
                               const Matrix& A);

template <libint2::OneBodyEngine::operator_type obtype>
std::array<Matrix, libint2::OneBodyEngine::operator_traits<obtype>::nopers>
compute_1body_ints(const BasisSet& obs,
                   const std::vector<Atom>& atoms = std::vector<Atom>());

#if LIBINT2_DERIV_ONEBODY_ORDER
template <libint2::OneBodyEngine::operator_type obtype>
std::vector<Matrix>
compute_1body_deriv_ints(unsigned deriv_order,
                         const BasisSet& obs,
                         const std::vector<Atom>& atoms);
#endif

Matrix compute_schwartz_ints(const BasisSet& bs1,
                             const BasisSet& bs2 = BasisSet(),
                             bool use_2norm = false // use infty norm by default
                            );
Matrix compute_do_ints(const BasisSet& bs1,
                       const BasisSet& bs2 = BasisSet(),
                       bool use_2norm = false // use infty norm by default
                      );

using shellpair_list_t = std::unordered_map<size_t,std::vector<size_t>>;
shellpair_list_t obs_shellpair_list; // shellpair list for OBS

/// computes non-negligible shell pair list; shells \c i and \c j form a non-negligible
/// pair if they share a center or the Frobenius norm of their overlap is greater than threshold
shellpair_list_t
compute_shellpair_list(const BasisSet& bs1,
                       const BasisSet& bs2 = BasisSet(),
                       double threshold = 1e-12
                      );

Matrix compute_2body_fock(const BasisSet& obs,
                          const Matrix& D,
                          double precision = std::numeric_limits<double>::epsilon(), // discard contributions smaller than this
                          const Matrix& Schwartz = Matrix() // K_ij = sqrt(||(ij|ij)||_\infty); if empty, do not Schwartz screen
                         );
// an Fock builder that can accept densities expressed a separate basis
Matrix compute_2body_fock_general(const BasisSet& obs,
                                  const Matrix& D,
                                  const BasisSet& D_bs,
                                  bool D_is_sheldiagonal = false, // set D_is_shelldiagonal if doing SOAD
                                  double precision = std::numeric_limits<double>::epsilon() // discard contributions smaller than this
                                 );

// returns {X,X^{-1},S_condition_number_after_conditioning}, where
// X is the generalized square-root-inverse such that X.transpose() * S * X = I
// columns of Xinv is the basis conditioned such that
// the condition number of its metric (Xinv.transpose . Xinv) < S_condition_number_threshold
std::tuple<Matrix,Matrix,double>
conditioning_orthogonalizer(const Matrix& S, double S_condition_number_threshold);

#ifdef LIBINT2_HAVE_BTAS
# define HAVE_DENSITY_FITTING 1
  struct DFFockEngine {
    const BasisSet& obs;
    const BasisSet& dfbs;
    DFFockEngine(const BasisSet& _obs, const BasisSet& _dfbs) :
      obs(_obs), dfbs(_dfbs)
    {
    }

    typedef btas::RangeNd<CblasRowMajor, std::array<long, 3> > Range3d;
    typedef btas::Tensor<double, Range3d> Tensor3d;
    Tensor3d xyK;

    // a DF-based builder, using coefficients of occupied MOs
    Matrix compute_2body_fock_dfC(const Matrix& Cocc);
  };
#endif // HAVE_DENSITY_FITTING

namespace libint2 {
  int nthreads;

  /// fires off \c nthreads instances of lambda in parallel
  template <typename Lambda>
  void parallel_do(Lambda& lambda) {
#ifdef _OPENMP
  #pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else // use C++11 threads
  std::vector<std::thread> threads;
  for(int thread_id=0; thread_id != libint2::nthreads; ++thread_id) {
    if(thread_id != nthreads-1)
      threads.push_back(std::thread(lambda,
                                    thread_id));
    else
      lambda(thread_id);
  } // threads_id
  for(int thread_id=0; thread_id<nthreads-1; ++thread_id)
    threads[thread_id].join();
#endif
  }

}

int main(int argc, char *argv[]) {

  using std::cout;
  using std::cerr;
  using std::endl;

  try {

    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    const auto basisname = (argc > 2) ? argv[2] : "aug-cc-pVDZ";
    bool do_density_fitting = false;
#ifdef HAVE_DENSITY_FITTING
    do_density_fitting = (argc > 3);
    const auto dfbasisname = do_density_fitting ? argv[3] : "";
#endif
    std::vector<Atom> atoms = read_geometry(filename);

    // set up thread pool
    {
      using libint2::nthreads;
      auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
      nthreads = 1;
      if (nthreads_cstr && strcmp(nthreads_cstr,"")) {
        std::istringstream iss(nthreads_cstr);
        iss >> nthreads;
        if (nthreads > 1<<16 || nthreads <= 0)
          nthreads = 1;
      }
#if defined(_OPENMP)
      omp_set_num_threads(nthreads);
#endif
      std::cout << "Will scale over " << nthreads
#if defined(_OPENMP)
                << " OpenMP"
#else
                << " C++11"
#endif
                << " threads" << std::endl;
    }

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;
    cout << "# of electrons = " << nelectron << endl;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij*xij + yij*yij + zij*zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << std::setprecision(15) << enuc << endl;

    libint2::Shell::do_enforce_unit_normalization(false);

    cout << "Atomic Cartesian coordinates (a.u.):" << endl;
    for(const auto& a: atoms)
      std::cout << a.atomic_number << " " << a.x << " " << a.y << " " << a.z << std::endl;

    BasisSet obs(basisname, atoms);
    cout << "orbital basis set rank = " << obs.nbf() << endl;

#ifdef HAVE_DENSITY_FITTING
    BasisSet dfbs;
    if (do_density_fitting) {
      dfbs = BasisSet(dfbasisname, atoms);
      cout << "density-fitting basis set rank = " << dfbs.nbf() << endl;
    }
    DFFockEngine dffockengine(obs,dfbs);
#endif // HAVE_DENSITY_FITTING


    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::init();

    // compute OBS non-negligible shell-pair list
    {
      obs_shellpair_list = compute_shellpair_list(obs);
      size_t nsp = 0;
      for(auto& sp: obs_shellpair_list) {
        nsp += sp.second.size();
      }
      std::cout << "# of {all,non-negligible} shell-pairs = {"
                << obs.size()*(obs.size()+1)/2 << ","
                << nsp << "}" << std::endl;
    }

    // compute one-body integrals
    auto S = compute_1body_ints<libint2::OneBodyEngine::overlap>(obs)[0];
    auto T = compute_1body_ints<libint2::OneBodyEngine::kinetic>(obs)[0];
    auto V = compute_1body_ints<libint2::OneBodyEngine::nuclear>(obs, atoms)[0];
    Matrix H = T + V;
    T.resize(0,0);
    V.resize(0,0);

    // compute orthogonalizer X such that X.transpose() . S . X = I
    Matrix X, Xinv;
    double XtX_condition_number; // condition number of "re-conditioned" overlap obtained as Xinv.transpose() . Xinv
                                 // one should think of columns of Xinv as the conditioned basis
                                 // Re: name ... cond # (Xinv.transpose() . Xinv) = cond # (X.transpose() . X)
    // by default assume can manage to compute with condition number of S <= 1/eps
    // this is probably too optimistic, but in well-behaved cases even 10^11 is OK
    double S_condition_number_threshold = 1.0/std::numeric_limits<double>::epsilon();
    std::tie(X, Xinv, XtX_condition_number) = conditioning_orthogonalizer(S,
        S_condition_number_threshold);

    Matrix D;
    Matrix C_occ;
    Matrix evals;
    {  // use SOAD as the guess density
      const auto tstart = std::chrono::high_resolution_clock::now();

      auto D_minbs = compute_soad(atoms); // compute guess in minimal basis
      BasisSet minbs("STO-3G", atoms);
      if (minbs == obs)
        D = D_minbs;
      else { // if basis != minimal basis, map non-representable SOAD guess into the AO basis
             // by diagonalizing a Fock matrix
        std::cout << "projecting SOAD into AO basis ... ";
        auto F = H;
        F += compute_2body_fock_general(obs, D_minbs, minbs,
                                        true /* SOAD_D_is_shelldiagonal */,
                                        std::numeric_limits<double>::epsilon() // this is cheap, no reason to be cheaper
                                       );

        // solve F C = e S C by (conditioned) transformation to F' C' = e C', where
        // F' = X.transpose() . F . X; the original C is obtained as C = X . C'
        Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F * X);
        auto C = X * eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        C_occ = C.leftCols(ndocc);
        D = C_occ * C_occ.transpose();

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;
        std::cout << "done (" << time_elapsed.count() << " s)" << std::endl;

      }
    }

    // pre-compute data for Schwartz bounds
    auto K = compute_schwartz_ints(obs);

    /*** =========================== ***/
    /***          SCF loop           ***/
    /*** =========================== ***/

    const auto maxiter = 100;
    const auto conv = 1e-12;
    auto iter = 0;
    auto rms_error = 1.0;
    auto ediff_rel = 0.0;
    auto ehf = 0.0;
    auto n2 = D.cols() * D.rows();
    libint2::DIIS<Matrix> diis(2); // start DIIS on second iteration

    // prepare for incremental Fock build ...
    Matrix D_diff = D;
    Matrix F = H;
    bool reset_incremental_fock_formation = false;
    bool incremental_Fbuild_started = false;
    double start_incremental_F_threshold = 1e-5;
    double next_reset_threshold = 0.0;
    size_t last_reset_iteration = 0;
    // ... unless doing DF, then use MO coefficients, hence not "incremental"
    if (do_density_fitting) start_incremental_F_threshold = 0.0;

    do {
      const auto tstart = std::chrono::high_resolution_clock::now();
      ++iter;

      // Last iteration's energy and density
      auto ehf_last = ehf;
      Matrix D_last = D;

      if (not incremental_Fbuild_started && rms_error < start_incremental_F_threshold) {
        incremental_Fbuild_started = true;
        reset_incremental_fock_formation = false;
        last_reset_iteration = iter - 1;
        next_reset_threshold = rms_error / 1e1;
        std::cout << "== started incremental fock build" << std::endl;
      }
      if (reset_incremental_fock_formation || not incremental_Fbuild_started) {
        F = H;
        D_diff = D;
      }
      if (reset_incremental_fock_formation && incremental_Fbuild_started) {
          reset_incremental_fock_formation = false;
          last_reset_iteration = iter;
          next_reset_threshold = rms_error / 1e1;
          std::cout << "== reset incremental fock build" << std::endl;
      }

      // build a new Fock matrix
      if (not do_density_fitting) {
        // totally empirical precision variation, involves the condition number
        const auto precision_F = std::min(std::min(1e-3 / XtX_condition_number, 1e-7),
            std::max(rms_error/1e4,std::numeric_limits<double>::epsilon()));
        F += compute_2body_fock(obs, D_diff, precision_F, K);
      }
#if HAVE_DENSITY_FITTING
      else { // do DF
        F = H + dffockengine.compute_2body_fock_dfC(C_occ);
      }
#else
      else { assert(false); } // do_density_fitting is true but HAVE_DENSITY_FITTING is not defined! should not happen
#endif // HAVE_DENSITY_FITTING

      // compute HF energy with the non-extrapolated Fock matrix
      ehf = D.cwiseProduct(H+F).sum();
      ediff_rel = std::abs((ehf - ehf_last)/ehf);

      // compute SCF error
      Matrix FD_comm = F*D*S - S*D*F;
      rms_error = FD_comm.norm()/n2;
      if (rms_error < next_reset_threshold || iter - last_reset_iteration >= 8)
        reset_incremental_fock_formation = true;

      // DIIS extrapolate F
      Matrix F_diis = F; // extrapolated F cannot be used in incremental Fock build; only used to produce the density
                         // make a copy of the unextrapolated matrix
      diis.extrapolate(F_diis,FD_comm);

      // solve F C = e S C by (conditioned) transformation to F' C' = e C', where
      // F' = X.transpose() . F . X; the original C is obtained as C = X . C'
      Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F_diis * X);
      evals = eig_solver.eigenvalues();
      auto C = X * eig_solver.eigenvectors();


      // compute density, D = C(occ) . C(occ)T
      C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
      D_diff = D - D_last;

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        std::cout <<
        "\n\nIter         E(HF)                 D(E)/E         RMS([F,D])/nn       Time(s)\n";
      printf(" %02d %20.12f %20.12e %20.12e %10.5lf\n", iter, ehf + enuc,
             ediff_rel, rms_error, time_elapsed.count());

    } while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));

    auto Mu = compute_1body_ints<libint2::OneBodyEngine::emultipole2>(obs);
    std::array<double,3> mu;
    std::array<double,6> qu;
    for(int xyz=0; xyz!=3; ++xyz)
      mu[xyz] = -2 * D.cwiseProduct(Mu[xyz+1]).sum(); // 2 = alpha + beta, -1 = electron charge
    for(int k=0; k!=6; ++k)
      qu[k] = -2 * D.cwiseProduct(Mu[k+4]).sum(); // 2 = alpha + beta, -1 = electron charge
    std::cout << "** edipole = "; std::copy(mu.begin(), mu.end(), std::ostream_iterator<double>(std::cout, " ")); std::cout << std::endl;
    std::cout << "** equadrupole = ";std::copy(qu.begin(), qu.end(), std::ostream_iterator<double>(std::cout, " ")); std::cout << std::endl;

#if LIBINT2_DERIV_ONEBODY_ORDER
    // compute forces
    {
      Matrix F1 = Matrix::Zero(atoms.size(), 3);
      Matrix F_Pulay = Matrix::Zero(atoms.size(), 3);
      //////////
      // one-body contributions to the forces
      //////////
      auto T1 = compute_1body_deriv_ints<libint2::OneBodyEngine::kinetic>(1, obs, atoms);
      auto V1 = compute_1body_deriv_ints<libint2::OneBodyEngine::nuclear>(1, obs, atoms);
      for(auto atom=0, i=0; atom!=atoms.size(); ++atom) {
        for(auto xyz=0; xyz!=3; ++xyz, ++i) {
          auto force = 2 * (T1[i]+V1[i]).cwiseProduct(D).sum();
          F1(atom, xyz) += force;
//        std::cout << "one-body force=" << force << std::endl;
//        std::cout << "derivative nuclear ints:\n" << V1[i] << std::endl;
        }
      }

      //////////
      // Pulay force
      //////////
      // orbital energy density
      DiagonalMatrix evals_occ(evals.topRows(ndocc));
      Matrix W = C_occ * evals_occ * C_occ.transpose();
      auto S1 = compute_1body_deriv_ints<libint2::OneBodyEngine::overlap>(1, obs, atoms);
      for(auto atom=0, i=0; atom!=atoms.size(); ++atom) {
        for(auto xyz=0; xyz!=3; ++xyz, ++i) {
          auto force = 2 * S1[i].cwiseProduct(W).sum();
          F_Pulay(atom, xyz) -= force;
//        std::cout << "Pulay force=" << force << std::endl;
//        std::cout << "derivative overlap ints:\n" << S1[i] << std::endl;
        }
      }

      std::cout << "** 1-body forces = ";
      for(int atom=0; atom!=atoms.size(); ++atom)
        for(int xyz=0; xyz!=3; ++xyz)
          std::cout << F1(atom,xyz) << " ";
      std::cout << std::endl;
      std::cout << "** Pulay forces = ";
      for(int atom=0; atom!=atoms.size(); ++atom)
        for(int xyz=0; xyz!=3; ++xyz)
          std::cout << F_Pulay(atom,xyz) << " ";
      std::cout << std::endl;
    }
#endif

    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);

    libint2::cleanup(); // done with libint

  } // end of try block; if any exceptions occurred, report them and exit cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
  catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }

  return 0;
}

std::vector<Atom> read_geometry(const std::string& filename) {

  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good());

  // to prepare for MPI parallelization, we will read the entire file into a string that can be
  // broadcast to everyone, then converted to an std::istringstream object that can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
  if ( filename.rfind(".xyz") != std::string::npos)
    return libint2::read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

// computes Superposition-Of-Atomic-Densities guess for the molecular density matrix
// in minimal basis; occupies subshells by smearing electrons evenly over the orbitals
Matrix compute_soad(const std::vector<Atom>& atoms) {

  // compute number of atomic orbitals
  size_t nao = 0;
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) // H, He
      nao += 1;
    else if (Z <= 10) // Li - Ne
      nao += 5;
    else
      throw "SOAD with Z > 10 is not yet supported";
  }

  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0; // first AO of this atom
  for(const auto& atom: atoms) {
    const auto Z = atom.atomic_number;
    if (Z == 1 || Z == 2) { // H, He
      D(ao_offset, ao_offset) = Z; // all electrons go to the 1s
      ao_offset += 1;
    }
    else if (Z <= 10) {
      D(ao_offset, ao_offset) = 2; // 2 electrons go to the 1s
      D(ao_offset+1, ao_offset+1) = (Z == 3) ? 1 : 2; // Li? only 1 electron in 2s, else 2 electrons
      // smear the remaining electrons in 2p orbitals
      const double num_electrons_per_2p = (Z > 4) ? (double)(Z - 4)/3 : 0;
      for(auto xyz=0; xyz!=3; ++xyz)
        D(ao_offset+2+xyz, ao_offset+2+xyz) = num_electrons_per_2p;
      ao_offset += 5;
    }
  }

  return D * 0.5; // we use densities normalized to # of electrons/2
}

Matrix compute_shellblock_norm(const BasisSet& obs,
                               const Matrix& A) {
  const auto nsh = obs.size();
  Matrix Ash(nsh, nsh);

  auto shell2bf = obs.shell2bf();
  for(size_t s1=0; s1!=nsh; ++s1) {
    const auto& s1_first = shell2bf[s1];
    const auto& s1_size = obs[s1].size();
    for(size_t s2=0; s2!=nsh; ++s2) {
      const auto& s2_first = shell2bf[s2];
      const auto& s2_size = obs[s2].size();

      Ash(s1, s2) = A.block(s1_first, s2_first, s1_size, s2_size).lpNorm<Eigen::Infinity>();
    }
  }

  return Ash;
}

template <libint2::OneBodyEngine::operator_type obtype>
std::array<Matrix, libint2::OneBodyEngine::operator_traits<obtype>::nopers>
compute_1body_ints(const BasisSet& obs,
                   const std::vector<Atom>& atoms)
{
  const auto n = obs.nbf();
  const auto nshells = obs.size();
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif
  typedef std::array<Matrix, libint2::OneBodyEngine::operator_traits<obtype>::nopers> result_type;
  const unsigned int nopers = libint2::OneBodyEngine::operator_traits<obtype>::nopers;
  result_type result; for(auto& r: result) r = Matrix::Zero(n,n);

  // construct the 1-body integrals engine
  std::vector<libint2::OneBodyEngine> engines(nthreads);
  engines[0] = libint2::OneBodyEngine(obtype, obs.max_nprim(), obs.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == libint2::OneBodyEngine::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engines[0].set_params(q);
  }
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for(auto s1=0l, s12=0l; s1!=nshells; ++s1) {

      auto bf1 = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();

      for(auto s2=0; s2<=s1; ++s2) {

        if (s12 % nthreads != thread_id)
          continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();

        auto n12 = n1 * n2;

        // compute shell pair; return is the pointer to the buffer
        const auto* buf = engines[thread_id].compute(obs[s1], obs[s2]);

        for(unsigned int op=0; op!=nopers; ++op, buf+=n12) {
          // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
          Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
          result[op].block(bf1, bf2, n1, n2) = buf_mat;
          if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
            result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }

      }
    }
  } // omp parallel

  return result;
}

#if LIBINT2_DERIV_ONEBODY_ORDER
template <libint2::OneBodyEngine::operator_type obtype>
std::vector<Matrix>
compute_1body_deriv_ints(unsigned deriv_order,
                         const BasisSet& obs,
                         const std::vector<Atom>& atoms)
{
  const auto n = obs.nbf();
  const auto nshells = obs.size();
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif
  constexpr auto nopers = libint2::OneBodyEngine::operator_traits<obtype>::nopers;
  const auto nresults = nopers * libint2::num_geometrical_derivatives(atoms.size(),deriv_order);
  typedef std::vector<Matrix> result_type;
  result_type result(nresults); for(auto& r: result) r = Matrix::Zero(n,n);

  // construct the 1-body integrals engine
  std::vector<libint2::OneBodyEngine> engines(nthreads);
  engines[0] = libint2::OneBodyEngine(obtype, obs.max_nprim(), obs.max_l(), deriv_order);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == libint2::OneBodyEngine::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engines[0].set_params(q);
  }
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for(auto s1=0l, s12=0l; s1!=nshells; ++s1) {

      auto bf1 = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();
      auto atom1 = shell2atom[s1];
      assert(atom1 != -1);

      for(auto s2=0; s2<=s1; ++s2) {

        if (s12 % nthreads != thread_id)
          continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();
        auto atom2 = shell2atom[s2];

        auto n12 = n1 * n2;

        // compute shell pair; return is the pointer to the buffer
        const auto* buf = engines[thread_id].compute(obs[s1], obs[s2]);

        assert(deriv_order == 1); // the loop structure below needs to be generalized for higher-order derivatives
        // 1. process derivatives with respect to the Gaussian origins first ...
        for(unsigned int d=0; d!=6; ++d) { // 2 centers x 3 axes = 6 cartesian geometric derivatives
          auto atom = d < 3 ? atom1 : atom2;
          auto op_start = (3*atom+d%3) * nopers;
          auto op_fence = op_start + nopers;
          for(unsigned int op=op_start; op!=op_fence; ++op, buf+=n12) {
            // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
            Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
//            std::cout << "s1=" << s1 << " s2=" << s2 << " x=" << (3*atom+d%3) << " op=" << op-op_start << ":\n";
//            std::cout << buf_mat << std::endl;
            result[op].block(bf1, bf2, n1, n2) += buf_mat;
            if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
              result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
          }
        }
        // 2. process derivatives of nuclear Coulomb operator, if needed
        if (obtype == libint2::OneBodyEngine::nuclear) {
          for(unsigned int atom=0; atom!=atoms.size(); ++atom) {
            for(unsigned int xyz=0; xyz!=3; ++xyz) {
              auto op_start = (3*atom+xyz) * nopers;
              auto op_fence = op_start + nopers;
              for(unsigned int op=op_start; op!=op_fence; ++op, buf+=n12) {
                Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
                result[op].block(bf1, bf2, n1, n2) += buf_mat;
                if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
                  result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
              }
            }
          }
        }

      }
    }
  } // omp parallel

  return result;
}
#endif

Matrix compute_schwartz_ints(const BasisSet& bs1,
                             const BasisSet& _bs2,
                             bool use_2norm) {
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  Matrix K = Matrix::Zero(nsh1,nsh2);
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(bs1.max_nprim(), bs2.max_l(), 0);
  engines[0].set_precision(0.); // !!! very important: cannot screen primitives in Schwartz computation !!!
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::cout << "computing Schwartz bound prerequisites ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s12=0l; s1!=nsh1; ++s1) {

      auto n1 = bs1[s1].size();// number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2-1;
      for(auto s2=0; s2<=s2_max; ++s2, ++s12) {

        if (s12 % nthreads != thread_id)
          continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1*n2;

        const auto* buf = engines[thread_id].compute(bs1[s1], bs2[s2], bs1[s1], bs2[s2]);

        // the diagonal elements are the Schwartz ints ... use Map.diagonal()
        Eigen::Map<const Matrix> buf_mat(buf, n12, n12);
        auto norm2 = use_2norm ? buf_mat.diagonal().norm() : buf_mat.diagonal().lpNorm<Eigen::Infinity>();
        K(s1,s2) = std::sqrt(norm2);
        if (bs1_equiv_bs2) K(s2,s1) = K(s1,s2);

      }
    }
  }

  timer.stop(0);
  std::cout << "done (" << timer.read(0) << " s)"<< std::endl;

  return K;
}

Matrix compute_do_ints(const BasisSet& bs1,
                       const BasisSet& _bs2,
                       bool use_2norm) {
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  Matrix K = Matrix::Zero(nsh1,nsh2);
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Delta> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(bs1.max_nprim(), bs2.max_l(), 0);
  engines[0].set_precision(0.); // !!! very important: cannot screen primitives in Schwartz computation !!!
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::cout << "computing DOIs ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s12=0l; s1!=nsh1; ++s1) {

      auto n1 = bs1[s1].size();// number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2-1;
      for(auto s2=0; s2<=s2_max; ++s2, ++s12) {

        if (s12 % nthreads != thread_id)
          continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1*n2;

        const auto* buf = engines[thread_id].compute(bs1[s1], bs2[s2], bs1[s1], bs2[s2]);

        // the diagonal elements are the Schwartz ints ... use Map.diagonal()
        Eigen::Map<const Matrix> buf_mat(buf, n12, n12);
        auto norm2 = use_2norm ? buf_mat.diagonal().norm() : buf_mat.diagonal().lpNorm<Eigen::Infinity>();
        K(s1,s2) = std::sqrt(norm2);
        if (bs1_equiv_bs2) K(s2,s1) = K(s1,s2);

      }
    }
  }

  timer.stop(0);
  std::cout << "done (" << timer.read(0) << " s)"<< std::endl;

  return K;
}

shellpair_list_t
compute_shellpair_list(const BasisSet& bs1,
                       const BasisSet& _bs2,
                       const double threshold) {
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif

  // construct the 2-electron repulsion integrals engine
  using libint2::OneBodyEngine;
  std::vector<OneBodyEngine> engines; engines.reserve(nthreads);
  engines.emplace_back(OneBodyEngine::overlap,
                       std::max(bs1.max_nprim(),bs2.max_nprim()),
                       std::max(bs1.max_l(),bs2.max_l()),
                       0);
  for(size_t i=1; i!=nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  std::cout << "computing non-negligible shell-pair list ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  shellpair_list_t result;

  std::mutex mx;

  auto compute = [&] (int thread_id) {

    auto& engine = engines[thread_id];

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s12=0l; s1!=nsh1; ++s1) {

      mx.lock();
      if (result.find(s1) == result.end())
        result.insert(std::make_pair(s1,std::vector<size_t>()));
      mx.unlock();

      auto n1 = bs1[s1].size();// number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2-1;
      for(auto s2=0; s2<=s2_max; ++s2, ++s12) {

        if (s12 % nthreads != thread_id)
          continue;

        auto on_same_center = (bs1[s1].O == bs2[s2].O);
        bool significant = on_same_center;
        if (not on_same_center) {
          auto n2 = bs2[s2].size();
          const auto* buf = engines[thread_id].compute(bs1[s1], bs2[s2]);
          Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
          auto norm = buf_mat.norm();
          significant = (norm >= threshold);
        }

        if (significant) {
          mx.lock();
          result[s1].emplace_back(s2);
          mx.unlock();
        }
      }
    }
  }; // end of compute

  libint2::parallel_do(compute);

  // resort shell list in increasing order, i.e. result[s][s1] < result[s][s2] if s1 < s2
  auto sort = [&] (int thread_id) {
    for(auto s1=0l; s1!=nsh1; ++s1) {
      if (s1%nthreads == thread_id) {
        auto& list = result[s1];
        std::sort(list.begin(), list.end());
      }
    }
  }; // end of sort

  libint2::parallel_do(sort);

  timer.stop(0);
  std::cout << "done (" << timer.read(0) << " s)"<< std::endl;

  return result;
}

// returns {X,X^{-1},rank,A_condition_number,result_A_condition_number}, where
// X is the generalized square-root-inverse such that X.transpose() * A * X = I
//
// if symmetric is true, produce "symmetric" sqrtinv: X = U . A_evals_sqrtinv . U.transpose()),
// else produce "canonical" sqrtinv: X = U . A_evals_sqrtinv
// where U are eigenvectors of A
// rows and cols of symmetric X are equivalent; for canonical X the rows are original basis (AO),
// cols are transformed basis ("orthogonal" AO)
//
// A is conditioned to max_condition_number
std::tuple<Matrix,Matrix,size_t,double,double>
gensqrtinv(const Matrix& S, bool symmetric=false, double max_condition_number=1e8) {
  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(S);
  auto U = eig_solver.eigenvectors();
  auto s = eig_solver.eigenvalues();
  auto s_max = s.maxCoeff();
  auto condition_number = std::min(s_max / std::max(s.minCoeff(),std::numeric_limits<double>::min()),
                                   1.0 / std::numeric_limits<double>::epsilon());
  auto threshold = s_max / max_condition_number;
  long n = s.rows();
  long n_cond=0;
  for(long i=n-1; i>=0; --i) {
    if (s(i) >= threshold) {
      ++n_cond;
    }
    else
      i = 0; // skip rest since eigenvalues are in ascending order
  }

  auto sigma = s.bottomRows(n_cond);
  auto result_condition_number = sigma.maxCoeff() / sigma.minCoeff();
  auto sigma_sqrt = sigma.array().sqrt().matrix().asDiagonal();
  auto sigma_invsqrt = sigma.array().sqrt().inverse().matrix().asDiagonal();

  // make canonical X/Xinv
  auto U_cond = U.block(0,n-n_cond,n,n_cond);
  Matrix X =  U_cond * sigma_invsqrt;
  Matrix Xinv = U_cond * sigma_sqrt;
  // convert to symmetric, if needed
  if (symmetric) {
    X = X * U_cond.transpose();
    Xinv = Xinv * U_cond.transpose();
  }
  return std::make_tuple(X,Xinv,size_t(n_cond),condition_number,result_condition_number);
}

std::tuple<Matrix,Matrix,double>
conditioning_orthogonalizer(const Matrix& S, double S_condition_number_threshold) {
  size_t obs_rank;
  double S_condition_number;
  double XtX_condition_number;
  Matrix X, Xinv;

  assert(S.rows() == S.cols());

  std::tie(X, Xinv, obs_rank, S_condition_number, XtX_condition_number) =
      gensqrtinv(S, false, S_condition_number_threshold);
  auto obs_nbf_omitted = (long)S.rows() - (long)obs_rank;
  std::cout << "overlap condition number = " << S_condition_number;
  if (obs_nbf_omitted > 0)
    std::cout << " (dropped " << obs_nbf_omitted
              << " " << (obs_nbf_omitted > 1 ? "fns" : "fn")
              << " to reduce to " << XtX_condition_number << ")";
  std::cout << std::endl;

  if (obs_nbf_omitted > 0) {
    Matrix should_be_I = X.transpose() * S * X;
    Matrix I = Matrix::Identity(should_be_I.rows(),should_be_I.cols());
    std::cout << "||X^t * S * X - I||_2 = " << (should_be_I - I).norm()
              << " (should be 0)" << std::endl;
  }

  return std::make_tuple(X,Xinv,XtX_condition_number);
}

Matrix compute_2body_2index_ints(const BasisSet& bs)
{
  const auto n = bs.nbf();
  const auto nshells = bs.size();
#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif
  Matrix result = Matrix::Zero(n,n);

  // build engines for each thread
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(bs.max_nprim(), bs.max_l(), 0);
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = bs.shell2bf();
  auto unitshell = Shell::unit();

#ifdef _OPENMP
  #pragma omp parallel
#endif
  {
#ifdef _OPENMP
    auto thread_id = omp_get_thread_num();
#else
    auto thread_id = 0;
#endif

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for(auto s1=0l, s12=0l; s1!=nshells; ++s1) {

      auto bf1 = shell2bf[s1]; // first basis function in this shell
      auto n1 = bs[s1].size();

      for(auto s2=0; s2<=s1; ++s2) {

        if (s12 % nthreads != thread_id)
          continue;

        auto bf2 = shell2bf[s2];
        auto n2 = bs[s2].size();

        // compute shell pair; return is the pointer to the buffer
        const auto* buf = engines[thread_id].compute(bs[s1], unitshell, bs[s2], unitshell);

        // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
        Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
        result.block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

      }
    }
  } // omp parallel

  return result;
}

Matrix compute_2body_fock(const BasisSet& obs,
                          const Matrix& D,
                          double precision,
                          const Matrix& Schwartz) {

  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  std::vector<Matrix> G(nthreads, Matrix::Zero(n,n));

  const auto do_schwartz_screen = Schwartz.cols() != 0 && Schwartz.rows() != 0;
  Matrix D_shblk_norm; // matrix of infty-norms of shell blocks
  if (do_schwartz_screen) {
    D_shblk_norm = compute_shellblock_norm(obs, D);
  }

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim*max_nprim*max_nprim*max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) / max_nprim4;

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(obs.max_nprim(), obs.max_l(), 0);
  engines[0].set_precision(engine_precision); // shellset-dependent precision control will likely break positive definiteness
                                       // stick with this simple recipe
  std::cout << "compute_2body_fock:precision = " << precision << std::endl;
  std::cout << "TwoBodyEngine::precision = " << engines[0].precision() << std::endl;
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();

  auto lambda = [&] (int thread_id) {

    auto& engine = engines[thread_id];
    auto& g = G[thread_id];

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s1234=0l; s1!=nshells; ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();// number of basis functions in this shell

      for(const auto& s2: obs_shellpair_list[s1]) {

        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        const auto Dnorm12 = do_schwartz_screen ? D_shblk_norm(s1,s2) : 0.;

        for(auto s3=0; s3<=s1; ++s3) {

          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();

          const auto Dnorm123 = do_schwartz_screen ? std::max(D_shblk_norm(s1,s3),
                                                              std::max(D_shblk_norm(s2,s3),Dnorm12)
                                                             )
                                                   : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for(const auto& s4: obs_shellpair_list[s3]) {
            if (s4 > s4_max) break; // for each s3, s4 are stored in monotonically increasing order

            if ((s1234++) % nthreads != thread_id)
              continue;

            const auto Dnorm1234 = do_schwartz_screen ? std::max(D_shblk_norm(s1,s4),
                                                                 std::max(D_shblk_norm(s2,s4),
                                                                          std::max(D_shblk_norm(s3,s4),Dnorm123)
                                                                         )
                                                                )
                                                      : 0.;

            if (do_schwartz_screen && Dnorm1234 * Schwartz(s1,s2) * Schwartz(s3,s4) < fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();

            num_ints_computed += n1*n2*n3*n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            const auto* buf = engine.compute(obs[s1], obs[s2], obs[s3], obs[s4]);

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            // 1) each shell set of integrals contributes up to 6 shell sets of the Fock matrix:
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    g(bf1,bf2) += D(bf3,bf4) * value_scal_by_deg;
                    g(bf3,bf4) += D(bf1,bf2) * value_scal_by_deg;
                    g(bf1,bf3) -= 0.25 * D(bf2,bf4) * value_scal_by_deg;
                    g(bf2,bf4) -= 0.25 * D(bf1,bf3) * value_scal_by_deg;
                    g(bf1,bf4) -= 0.25 * D(bf2,bf3) * value_scal_by_deg;
                    g(bf2,bf3) -= 0.25 * D(bf1,bf4) * value_scal_by_deg;
                  }
                }
              }
            }

          }
        }
      }
    }

  }; // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for(size_t i=1; i!=nthreads; ++i) {
    G[0] += G[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for(auto& t: timers) {
    time_for_ints += t.read(0);
  }
  std::cout << "time for integrals = " << time_for_ints << std::endl;
  for(int t=0; t!=nthreads; ++t)
    engines[t].print_timers();
#endif

  Matrix GG = 0.5 * (G[0] + G[0].transpose());

  std::cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  return GG;
}

Matrix compute_2body_fock_general(const BasisSet& obs,
                                  const Matrix& D,
                                  const BasisSet& D_bs,
                                  bool D_is_shelldiagonal,
                                  double precision) {

  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto n_D = D_bs.nbf();
  assert(D.cols() == D.rows() && D.cols() == n_D);

  using libint2::nthreads;
  std::vector<Matrix> G(nthreads, Matrix::Zero(n,n));

  // construct the 2-electron repulsion integrals engine
  typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
  std::vector<coulomb_engine_type> engines(nthreads);
  engines[0] = coulomb_engine_type(std::max(obs.max_nprim(),D_bs.max_nprim()),
                                   std::max(obs.max_l(), D_bs.max_l()), 0);
  engines[0].set_precision(precision); // shellset-dependent precision control will likely break positive definiteness
                                       // stick with this simple recipe
  for(size_t i=1; i!=nthreads; ++i) {
    engines[i] = engines[0];
  }
  auto shell2bf = obs.shell2bf();
  auto shell2bf_D = D_bs.shell2bf();

  auto lambda = [&] (int thread_id) {

    auto& engine = engines[thread_id];
    auto& g = G[thread_id];

    // loop over permutationally-unique set of shells
    for(auto s1=0l, s1234=0l; s1!=nshells; ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = obs[s1].size();   // number of basis functions in this shell

      for(auto s2=0; s2<=s1; ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        for(auto s3=0; s3<D_bs.size(); ++s3) {

          auto bf3_first = shell2bf_D[s3];
          auto n3 = D_bs[s3].size();

          auto s4_begin = D_is_shelldiagonal ? s3 : 0;
          auto s4_fence = D_is_shelldiagonal ? s3+1 : D_bs.size();

          for(auto s4=s4_begin; s4!=s4_fence; ++s4, ++s1234) {

            if (s1234 % nthreads != thread_id)
              continue;

            auto bf4_first = shell2bf_D[s4];
            auto n4 = D_bs[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

            if (s3 >= s4) {
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s1234_deg = s12_deg * s34_deg;
              //auto s1234_deg = s12_deg;
              const auto* buf_J = engine.compute(obs[s1], obs[s2], D_bs[s3], D_bs[s4]);

              for(auto f1=0, f1234=0; f1!=n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for(auto f2=0; f2!=n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for(auto f3=0; f3!=n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      const auto value = buf_J[f1234];
                      const auto value_scal_by_deg = value * s1234_deg;
                      g(bf1,bf2) += 2.0 * D(bf3,bf4) * value_scal_by_deg;
                    }
                  }
                }
              }
            }

            const auto* buf_K = engine.compute(obs[s1], D_bs[s3], obs[s2], D_bs[s4]);

            for(auto f1=0, f1324=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f3=0; f3!=n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for(auto f2=0; f2!=n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1324) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_K[f1324];
                    const auto value_scal_by_deg = value * s12_deg;
                    g(bf1,bf2) -= D(bf3,bf4) * value_scal_by_deg;
                  }
                }
              }
            }

          }
        }
      }
    }

  }; // thread lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for(size_t i=1; i!=nthreads; ++i) {
    G[0] += G[i];
  }

  // symmetrize the result and return
  return 0.5 * (G[0] + G[0].transpose());
}

#ifdef HAVE_DENSITY_FITTING

// uncomment if want to compute statistics of shell blocks
//#define COMPUTE_DF_INTS_STATS 1

Matrix
DFFockEngine::compute_2body_fock_dfC(const Matrix& Cocc) {

#ifdef _OPENMP
  const auto nthreads = omp_get_max_threads();
#else
  const auto nthreads = 1;
#endif

  const auto n          =  obs.nbf();
  const auto ndf        = dfbs.nbf();
#ifndef _OPENMP
  libint2::Timers<5> timer;
  timer.set_now_overhead(25);
#endif // not defined _OPENMP

  typedef btas::RangeNd<CblasRowMajor, std::array<long, 1> > Range1d;
  typedef btas::RangeNd<CblasRowMajor, std::array<long, 2> > Range2d;
  typedef btas::Tensor<double, Range1d> Tensor1d;
  typedef btas::Tensor<double, Range2d> Tensor2d;

  // using first time? compute 3-center ints and transform to inv sqrt repreentation
  if (xyK.size() == 0) {

    const auto nshells    =  obs.size();
    const auto nshells_df = dfbs.size();
    const auto unitshell = libint2::Shell::unit();

    // construct the 2-electron 3-center repulsion integrals engine
    typedef libint2::TwoBodyEngine<libint2::Coulomb> coulomb_engine_type;
    std::vector<coulomb_engine_type> engines(nthreads);
    engines[0] = coulomb_engine_type(std::max(obs.max_nprim(), dfbs.max_nprim()),
                                     std::max(obs.max_l(), dfbs.max_l()), 0);
    for(size_t i=1; i!=nthreads; ++i) {
      engines[i] = engines[0];
    }

    auto shell2bf    =  obs.shell2bf();
    auto shell2bf_df = dfbs.shell2bf();

    Tensor3d xyZ{n, n, ndf};

#if COMPUTE_DF_INTS_STATS
    const double count_threshold = 1e-7;
    std::atomic<size_t> nints{0};  // number of ints in shell blocks with Frobenius norm > count_threshold
    std::atomic<size_t> nints2{0}; // number of ints in shell blocks with Frobenius norm/elem > count_threshold
    typedef btas::Tensor<char, Range3d> Tensor3dBool;
    const size_t nobs_per_mol = 24; // cc-pVDZ
    const size_t ndfbs_per_mol = 84; // cc-pVDZ-RI
    const size_t nmols = n/nobs_per_mol;
    Tensor3dBool nonzero_mol_blocks{nmols, nmols, nmols}; // number of molecule blocks that contain shell blocks with Frobenius norm > count_threshold
    Tensor3d mol_blocks_frob2{nmols, nmols, nmols}; // sum of Frobenius norms squared in each molecule block
    std::fill(nonzero_mol_blocks.begin(), nonzero_mol_blocks.end(), 0);
    std::fill(mol_blocks_frob2.begin(), mol_blocks_frob2.end(), 0.);
#endif // COMPUTE_DF_INTS_STATS

  #ifdef _OPENMP
    #pragma omp parallel
  #endif
    {
  #ifdef _OPENMP
      auto thread_id = omp_get_thread_num();
  #else
      auto thread_id = 0;
  #endif

      // loop over permutationally-unique set of shells
      for(auto s1=0l, s123=0l; s1!=nshells; ++s1) {

        auto bf1_first = shell2bf[s1]; // first basis function in this shell
        auto n1 = obs[s1].size();// number of basis functions in this shell

        for(auto s2=0; s2!=nshells; ++s2) {

          auto bf2_first = shell2bf[s2];
          auto n2 = obs[s2].size();
          const auto n12 = n1*n2;

          for(auto s3=0; s3!=nshells_df; ++s3, ++s123) {

            if (s123 % nthreads != thread_id)
              continue;

            auto bf3_first = shell2bf_df[s3];
            auto n3 = dfbs[s3].size();
            const auto n123 = n12*n3;

  #ifndef _OPENMP
            timer.start(0);
  #endif

            const auto* buf = engines[thread_id].compute(obs[s1], obs[s2], dfbs[s3], unitshell);

  #ifndef _OPENMP
            timer.stop(0);
  #endif

#if COMPUTE_DF_INTS_STATS
            const auto buf_2norm = sqrt(std::inner_product(buf, buf+n123, buf, 0.));
            const auto buf_2norm_scaled = sqrt(buf_2norm * buf_2norm/ n123);
            const auto buf_infnorm = std::abs(
                  *std::max_element(buf, buf+n123, [](double a, double b){return std::abs(a) < std::abs(b);})
            );
            if (buf_2norm > count_threshold) {
              nints += n123;

              nonzero_mol_blocks(bf1_first/nobs_per_mol,
                                 bf2_first/nobs_per_mol,
                                 bf3_first/ndfbs_per_mol) = 1;
            }
            if (buf_2norm_scaled > count_threshold) {
              nints2 += n123;
            }
            mol_blocks_frob2(bf1_first/nobs_per_mol,
                             bf2_first/nobs_per_mol,
                             bf3_first/ndfbs_per_mol) += buf_2norm*buf_2norm;
#endif // COMPUTE_DF_INTS_STATS

  #ifndef _OPENMP
            timer.start(1);
  #endif

            auto lower_bound = {bf1_first, bf2_first, bf3_first};
            auto upper_bound = {bf1_first+n1, bf2_first+n2, bf3_first+n3};
            auto view = btas::make_view( xyZ.range().slice(lower_bound, upper_bound),
                                         xyZ.storage());
            std::copy(buf, buf+n123, view.begin());

  #ifndef _OPENMP
            timer.stop(1);
  #endif

          } // s3
        } // s2
      } // s1

    } // omp parallel

  #ifndef _OPENMP
    std::cout << "time for integrals = " << timer.read(0) << std::endl;
    std::cout << "time for copying into BTAS = " << timer.read(1) << std::endl;
    engines[0].print_timers();
  #endif // not defined _OPENMP

#if COMPUTE_DF_INTS_STATS
    {
      const auto nints_per_molblock = nobs_per_mol * nobs_per_mol * ndfbs_per_mol;
      const size_t nints_mols = std::count(nonzero_mol_blocks.begin(), nonzero_mol_blocks.end(), 1) *
                                nints_per_molblock;
      const size_t nints_mols2= std::count_if(mol_blocks_frob2.begin(), mol_blocks_frob2.end(),
                                              [&](double a) { return sqrt(a) > count_threshold; }
                                             ) *
                                                 nints_per_molblock;
      const size_t nints_mols3= std::count_if(mol_blocks_frob2.begin(), mol_blocks_frob2.end(),
                                              [&](double a) { return sqrt(a)/nints_per_molblock > count_threshold; }
                                             ) *
                                                 nints_per_molblock;
      std::cout << "# of ints in shell blocks with norm greater than " << count_threshold << " = " << nints << std::endl;
      std::cout << "# of ints in shell blocks with scaled norm greater than " << count_threshold << " = " << nints2 << std::endl;
      std::cout << "# of ints in molecule blocks whose any shell-block norm was greater than " << count_threshold << " = " << nints_mols << std::endl;
      std::cout << "# of ints in molecule blocks with norm greater than " << count_threshold << " = " << nints_mols2 << std::endl;
      std::cout << "# of ints in molecule blocks with scaled norm greater than " << count_threshold << " = " << nints_mols3 << std::endl;
      std::cout << "# of total ints = " << n*n*ndf << std::endl;
    }
#endif // COMPUTE_DF_INTS_STATS

    timer.start(2);

    Matrix V = compute_2body_2index_ints(dfbs);
    Eigen::LLT<Matrix> V_LLt(V);
    Matrix I = Matrix::Identity(ndf, ndf);
    auto L = V_LLt.matrixL();
    Matrix V_L = L;
    Matrix Linv = L.solve(I).transpose();
    // check
  //  std::cout << "||V - L L^t|| = " << (V - V_L * V_L.transpose()).norm() << std::endl;
  //  std::cout << "||I - L L^-1^t|| = " << (I - V_L * Linv.transpose()).norm() << std::endl;
  //  std::cout << "||V^-1 - L^-1 L^-1^t|| = " << (V.inverse() - Linv * Linv.transpose()).norm() << std::endl;

    Tensor2d K{ndf, ndf};
    std::copy(Linv.data(), Linv.data()+ndf*ndf, K.begin());

    xyK = Tensor3d{n, n, ndf};
    btas::contract(1.0, xyZ, {1,2,3}, K, {3,4}, 0.0, xyK, {1,2,4});
    xyZ = Tensor3d{0,0,0}; // release memory

    timer.stop(2);
    std::cout << "time for integrals tform = " << timer.read(2) << std::endl;
  } // if (xyK.size() == 0)

  // compute exchange
  timer.start(3);

  const auto nocc = Cocc.cols();
  Tensor2d Co{n, nocc};
  std::copy(Cocc.data(), Cocc.data()+n*nocc, Co.begin());
  Tensor3d xiK{n, nocc, ndf};
  btas::contract(1.0, xyK, {1,2,3}, Co, {2,4}, 0.0, xiK, {1,4,3});

  Tensor2d G{n, n};
  btas::contract(1.0, xiK, {1,2,3}, xiK, {4,2,3}, 0.0, G, {1,4});

  timer.stop(3);
  std::cout << "time for exchange = " << timer.read(3) << std::endl;

  // compute Coulomb
  timer.start(4);

  Tensor1d Jtmp{ndf};
  btas::contract(1.0, xiK, {1,2,3}, Co, {1,2}, 0.0, Jtmp, {3});
  xiK = Tensor3d{0,0,0};
  btas::contract(2.0, xyK, {1,2,3}, Jtmp, {3}, -1.0, G, {1,2});

  timer.stop(4);
  std::cout << "time for coulomb = " << timer.read(4) << std::endl;

  // copy result to an Eigen::Matrix
  Matrix result(n, n);
  std::copy(G.cbegin(), G.cend(), result.data());
  return result;
}
#endif // HAVE_DENSITY_FITTING

// should be a unit test somewhere
void api_basic_compile_test(const BasisSet& obs) {

  libint2::OneBodyEngine onebody_engine(libint2::OneBodyEngine::overlap, // will compute overlap ints
                                        obs.max_nprim(), // max # of primitives in shells this engine will accept
                                        obs.max_l()      // max angular momentum of shells this engine will accept
                                       );
  auto shell2bf = obs.shell2bf();
  for(auto s1=0; s1!=obs.size(); ++s1) {
    for(auto s2=0; s2!=obs.size(); ++s2) {

      std::cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
      const auto* ints_shellset = onebody_engine.compute(obs[s1], obs[s2]);
      std::cout << "done" << std::endl;

      auto bf1 = shell2bf[s1];  // first basis function in first shell
      auto n1 = obs[s1].size(); // number of basis functions in first shell
      auto bf2 = shell2bf[s2];  // first basis function in second shell
      auto n2 = obs[s2].size(); // number of basis functions in second shell

      // this iterates over integrals in the order they are packed in array ints_shellset
      for(auto f1=0; f1!=n1; ++f1)
        for(auto f2=0; f2!=n2; ++f2)
          std::cout << "  " << bf1+f1 << " " << bf2+f2 << " " << ints_shellset[f1*n2+f2] << std::endl;
    }
  }

  typedef libint2::TwoBodyEngine<libint2::cGTG> contracted_gtg_engine_t;
  std::vector<std::pair<double,double>> cgtg_params{{0.1,0.2}, {-0.3,0.4}, {0.5, 0.6}};
  contracted_gtg_engine_t twobody_engine(obs.max_nprim(), obs.max_l(),
                                         0, // level of geometrical derivatives; 0 = regular integrals
                                         std::numeric_limits<double>::epsilon(), // target precision
                                         cgtg_params // geminal parameters
                                        );
  auto buf = twobody_engine.compute(obs[0], obs[1], obs[2], obs[3]);
  typedef Eigen::Map<const Matrix> ConstMap;
  ConstMap buf_mat(buf, obs[0].size() * obs[1].size(), obs[2].size() * obs[3].size());
  std::cout << buf_mat << std::endl;

  typedef libint2::TwoBodyEngine<libint2::cGTG_times_Coulomb> cgtg_coulomb_engine_t;
  cgtg_coulomb_engine_t twobody_engine3(obs.max_nprim(), obs.max_l(),
                                        0, std::numeric_limits<double>::epsilon(),
                                        cgtg_params);
  typedef libint2::TwoBodyEngine<libint2::DelcGTG_square> delcgtg2_engine_t;
  delcgtg2_engine_t twobody_engine4(obs.max_nprim(), obs.max_l(),
                                    0, std::numeric_limits<double>::epsilon(),
                                    cgtg_params);
}
