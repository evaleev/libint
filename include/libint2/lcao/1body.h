/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint library.
 *
 *  Libint library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libint library.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// standard C++ headers
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

// have BTAS library?
#ifdef LIBINT2_HAVE_BTAS
#include <btas/btas.h>
#else  // LIBINT2_HAVE_BTAS
#error "libint2::lcao requires BTAS"
#endif

// Libint Gaussian integrals library
#include <libint2/chemistry/sto3g_atomic_density.h>
#include <libint2/diis.h>
#include <libint2/util/intpart_iter.h>

#include <libint2.hpp>

#if defined(_OPENMP)
#include <omp.h>
#endif

typedef btas::RangeNd<CblasRowMajor, std::array<long, 2>> Range2;
typedef btas::RangeNd<CblasRowMajor, std::array<long, 3>> Range3;
typedef btas::RangeNd<CblasRowMajor, std::array<long, 4>> Range4;
typedef btas::Tensor<double, Range2> Tensor2d;
typedef btas::Tensor<double, Range3> Tensor3d;
typedef btas::Tensor<double, Range3> Tensor4d;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;  // import dense, dynamically sized Matrix type from Eigen;
             // this is a matrix with row-major storage
             // (http://en.wikipedia.org/wiki/Row-major_order)
// to meet the layout of the integrals returned by the Libint integral library
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic>
    DiagonalMatrix;

using libint2::BasisSet;
using libint2::BraKet;
using libint2::Operator;
using libint2::Shell;
using libint2::libint2::Atom;

std::vector<libint2::Atom> read_geometry(const std::string& filename);
Matrix compute_soad(const std::vector<libint2::Atom>& atoms);
// computes norm of shell-blocks of A
Matrix compute_shellblock_norm(const BasisSet& obs, const Matrix& A);

template <Operator obtype>
std::array<Matrix, libint2::operator_traits<obtype>::nopers> compute_1body_ints(
    const BasisSet& obs,
    const std::vector<libint2::Atom>& atoms = std::vector<libint2::Atom>());

#if LIBINT2_DERIV_ONEBODY_ORDER
template <Operator obtype>
std::vector<Matrix> compute_1body_ints_deriv(
    unsigned deriv_order, const BasisSet& obs,
    const std::vector<libint2::Atom>& atoms);
#endif  // LIBINT2_DERIV_ONEBODY_ORDER

template <libint2::Operator Kernel = libint2::Operator::coulomb>
Matrix compute_schwarz_ints(
    const BasisSet& bs1, const BasisSet& bs2 = BasisSet(),
    bool use_2norm = false,  // use infty norm by default
    typename libint2::operator_traits<Kernel>::oper_params_type params =
        libint2::operator_traits<Kernel>::default_params());
Matrix compute_do_ints(const BasisSet& bs1, const BasisSet& bs2 = BasisSet(),
                       bool use_2norm = false  // use infty norm by default
);

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;
shellpair_list_t obs_shellpair_list;  // shellpair list for OBS

/// computes non-negligible shell pair list; shells \c i and \c j form a
/// non-negligible
/// pair if they share a center or the Frobenius norm of their overlap is
/// greater than threshold
shellpair_list_t compute_shellpair_list(const BasisSet& bs1,
                                        const BasisSet& bs2 = BasisSet(),
                                        double threshold = 1e-12);

Matrix compute_2body_fock(
    const BasisSet& obs, const Matrix& D,
    double precision = std::numeric_limits<
        double>::epsilon(),           // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                      // empty, do not Schwarz screen
);
// an Fock builder that can accept densities expressed a separate basis
Matrix compute_2body_fock_general(
    const BasisSet& obs, const Matrix& D, const BasisSet& D_bs,
    bool D_is_sheldiagonal = false,  // set D_is_shelldiagonal if doing SOAD
    double precision = std::numeric_limits<
        double>::epsilon()  // discard contributions smaller than this
);

#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<Matrix> compute_2body_fock_deriv(
    const BasisSet& obs, const std::vector<libint2::Atom>& atoms,
    const Matrix& D,
    double precision = std::numeric_limits<
        double>::epsilon(),           // discard contributions smaller than this
    const Matrix& Schwarz = Matrix()  // K_ij = sqrt(||(ij|ij)||_\infty); if
                                      // empty, do not Schwarz screen
);
#endif  // LIBINT2_DERIV_ERI_ORDER

// returns {X,X^{-1},S_condition_number_after_conditioning}, where
// X is the generalized square-root-inverse such that X.transpose() * S * X = I
// columns of Xinv is the basis conditioned such that
// the condition number of its metric (Xinv.transpose . Xinv) <
// S_condition_number_threshold
std::tuple<Matrix, Matrix, double> conditioning_orthogonalizer(
    const Matrix& S, double S_condition_number_threshold);

#ifdef LIBINT2_HAVE_BTAS
#define HAVE_DENSITY_FITTING 1
struct DFFockEngine {
  const BasisSet& obs;
  const BasisSet& dfbs;
  DFFockEngine(const BasisSet& _obs, const BasisSet& _dfbs)
      : obs(_obs), dfbs(_dfbs) {}

  typedef btas::RangeNd<CblasRowMajor, std::array<long, 3>> Range3d;
  typedef btas::Tensor<double, Range3d> Tensor3d;
  Tensor3d xyK;

  // a DF-based builder, using coefficients of occupied MOs
  Matrix compute_2body_fock_dfC(const Matrix& Cocc);
};
#endif  // HAVE_DENSITY_FITTING

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
#else  // use C++11 threads
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id != libint2::nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  }  // threads_id
  for (int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}
}  // namespace libint2

int main(int argc, char* argv[]) {
  using std::cerr;
  using std::cout;
  using std::endl;

  try {
    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take
    // filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";
    const auto basisname = (argc > 2) ? argv[2] : "aug-cc-pVDZ";
    bool do_density_fitting = false;
#ifdef HAVE_DENSITY_FITTING
    do_density_fitting = (argc > 3);
    const auto dfbasisname = do_density_fitting ? argv[3] : "";
#endif
    std::vector<libint2::Atom> atoms = read_geometry(filename);

    // set up thread pool
    {
      using libint2::nthreads;
      auto nthreads_cstr = getenv("LIBINT_NUM_THREADS");
      nthreads = 1;
      if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
        std::istringstream iss(nthreads_cstr);
        iss >> nthreads;
        if (nthreads > 1 << 16 || nthreads <= 0) nthreads = 1;
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
    for (auto i = 0; i < atoms.size(); ++i) nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;
    cout << "# of electrons = " << nelectron << endl;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij * xij + yij * yij + zij * zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << std::setprecision(15) << enuc
         << endl;

    libint2::Shell::do_enforce_unit_normalization(false);

    cout << "libint2::Atomic Cartesian coordinates (a.u.):" << endl;
    for (const auto& a : atoms)
      std::cout << a.atomic_number << " " << a.x << " " << a.y << " " << a.z
                << std::endl;

    BasisSet obs(basisname, atoms);
    cout << "orbital basis set rank = " << obs.nbf() << endl;

#ifdef HAVE_DENSITY_FITTING
    BasisSet dfbs;
    if (do_density_fitting) {
      dfbs = BasisSet(dfbasisname, atoms);
      cout << "density-fitting basis set rank = " << dfbs.nbf() << endl;
    }
#endif  // HAVE_DENSITY_FITTING

    /*** =========================== ***/
    /*** compute 1-e integrals       ***/
    /*** =========================== ***/

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute OBS non-negligible shell-pair list
    {
      obs_shellpair_list = compute_shellpair_list(obs);
      size_t nsp = 0;
      for (auto& sp : obs_shellpair_list) {
        nsp += sp.second.size();
      }
      std::cout << "# of {all,non-negligible} shell-pairs = {"
                << obs.size() * (obs.size() + 1) / 2 << "," << nsp << "}"
                << std::endl;
    }

    // compute one-body integrals
    auto S = compute_1body_ints<Operator::overlap>(obs)[0];
    auto T = compute_1body_ints<Operator::kinetic>(obs)[0];
    auto V = compute_1body_ints<Operator::nuclear>(
        obs, libint2::make_point_charges(atoms))[0];
    Matrix H = T + V;
    T.resize(0, 0);
    V.resize(0, 0);

    // compute orthogonalizer X such that X.transpose() . S . X = I
    Matrix X, Xinv;
    double XtX_condition_number;  // condition number of "re-conditioned"
                                  // overlap obtained as Xinv.transpose() . Xinv
    // one should think of columns of Xinv as the conditioned basis
    // Re: name ... cond # (Xinv.transpose() . Xinv) = cond # (X.transpose() .
    // X)
    // by default assume can manage to compute with condition number of S <=
    // 1/eps
    // this is probably too optimistic, but in well-behaved cases even 10^11 is
    // OK
    double S_condition_number_threshold =
        1.0 / std::numeric_limits<double>::epsilon();
    std::tie(X, Xinv, XtX_condition_number) =
        conditioning_orthogonalizer(S, S_condition_number_threshold);

    Matrix D;
    Matrix C_occ;
    Matrix evals;
    {  // use SOAD as the guess density
      const auto tstart = std::chrono::high_resolution_clock::now();

      auto D_minbs = compute_soad(atoms);  // compute guess in minimal basis
      BasisSet minbs("STO-3G", atoms);
      if (minbs == obs)
        D = D_minbs;
      else {  // if basis != minimal basis, map non-representable SOAD guess
              // into the AO basis
              // by diagonalizing a Fock matrix
        std::cout << "projecting SOAD into AO basis ... ";
        auto F = H;
        F += compute_2body_fock_general(
            obs, D_minbs, minbs, true /* SOAD_D_is_shelldiagonal */,
            std::numeric_limits<double>::epsilon()  // this is cheap, no reason
                                                    // to be cheaper
        );

        // solve F C = e S C by (conditioned) transformation to F' C' = e C',
        // where
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

    // pre-compute data for Schwarz bounds
    auto K = compute_schwarz_ints<>(obs, BasisSet(), true);

// prepare for density fitting
#ifdef HAVE_DENSITY_FITTING
    std::unique_ptr<DFFockEngine> dffockengine(
        do_density_fitting ? new DFFockEngine(obs, dfbs) : nullptr);
#endif  // HAVE_DENSITY_FITTING

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
    libint2::DIIS<Matrix> diis(2);  // start DIIS on second iteration

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

      if (not incremental_Fbuild_started &&
          rms_error < start_incremental_F_threshold) {
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
        const auto precision_F = std::min(
            std::min(1e-3 / XtX_condition_number, 1e-7),
            std::max(rms_error / 1e4, std::numeric_limits<double>::epsilon()));
        F += compute_2body_fock(obs, D_diff, precision_F, K);
      }
#if HAVE_DENSITY_FITTING
      else {  // do DF
        F = H + dffockengine->compute_2body_fock_dfC(C_occ);
      }
#else
      else {
        assert(false);
      }  // do_density_fitting is true but HAVE_DENSITY_FITTING is not defined!
         // should not happen
#endif  // HAVE_DENSITY_FITTING

      // compute HF energy with the non-extrapolated Fock matrix
      ehf = D.cwiseProduct(H + F).sum();
      ediff_rel = std::abs((ehf - ehf_last) / ehf);

      // compute SCF error
      Matrix FD_comm = F * D * S - S * D * F;
      rms_error = FD_comm.norm() / n2;
      if (rms_error < next_reset_threshold || iter - last_reset_iteration >= 8)
        reset_incremental_fock_formation = true;

      // DIIS extrapolate F
      Matrix F_diis = F;  // extrapolated F cannot be used in incremental Fock
                          // build; only used to produce the density
                          // make a copy of the unextrapolated matrix
      diis.extrapolate(F_diis, FD_comm);

      // solve F C = e S C by (conditioned) transformation to F' C' = e C',
      // where
      // F' = X.transpose() . F . X; the original C is obtained as C = X . C'
      Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(X.transpose() * F_diis *
                                                       X);
      evals = eig_solver.eigenvalues();
      auto C = X * eig_solver.eigenvectors();

      // compute density, D = C(occ) . C(occ)T
      C_occ = C.leftCols(ndocc);
      D = C_occ * C_occ.transpose();
      D_diff = D - D_last;

      const auto tstop = std::chrono::high_resolution_clock::now();
      const std::chrono::duration<double> time_elapsed = tstop - tstart;

      if (iter == 1)
        std::cout << "\n\nIter         E(HF)                 D(E)/E         "
                     "RMS([F,D])/nn       Time(s)\n";
      printf(" %02d %20.12f %20.12e %20.12e %10.5lf\n", iter, ehf + enuc,
             ediff_rel, rms_error, time_elapsed.count());

    } while (((ediff_rel > conv) || (rms_error > conv)) && (iter < maxiter));

    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);

    auto Mu = compute_1body_ints<Operator::emultipole2>(obs);
    std::array<double, 3> mu;
    std::array<double, 6> qu;
    for (int xyz = 0; xyz != 3; ++xyz)
      mu[xyz] = -2 * D.cwiseProduct(Mu[xyz + 1])
                         .sum();  // 2 = alpha + beta, -1 = electron charge
    for (int k = 0; k != 6; ++k)
      qu[k] = -2 * D.cwiseProduct(Mu[k + 4])
                       .sum();  // 2 = alpha + beta, -1 = electron charge
    std::cout << "** edipole = ";
    std::copy(mu.begin(), mu.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "** equadrupole = ";
    std::copy(qu.begin(), qu.end(),
              std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    {  // compute force
#if LIBINT2_DERIV_ONEBODY_ORDER
       // compute 1-e forces
      Matrix F1 = Matrix::Zero(atoms.size(), 3);
      Matrix F_Pulay = Matrix::Zero(atoms.size(), 3);
      //////////
      // one-body contributions to the forces
      //////////
      auto T1 = compute_1body_ints_deriv<Operator::kinetic>(1, obs, atoms);
      auto V1 = compute_1body_ints_deriv<Operator::nuclear>(1, obs, atoms);
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          auto force = 2 * (T1[i] + V1[i]).cwiseProduct(D).sum();
          F1(atom, xyz) += force;
        }
      }

      //////////
      // Pulay force
      //////////
      // orbital energy density
      DiagonalMatrix evals_occ(evals.topRows(ndocc));
      Matrix W = C_occ * evals_occ * C_occ.transpose();
      auto S1 = compute_1body_ints_deriv<Operator::overlap>(1, obs, atoms);
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          auto force = 2 * S1[i].cwiseProduct(W).sum();
          F_Pulay(atom, xyz) -= force;
        }
      }

      std::cout << "** 1-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) std::cout << F1(atom, xyz) << " ";
      std::cout << std::endl;
      std::cout << "** Pulay forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz)
          std::cout << F_Pulay(atom, xyz) << " ";
      std::cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER

#if LIBINT2_DERIV_ERI_ORDER
      // compute 2-e forces
      Matrix F2 = Matrix::Zero(atoms.size(), 3);

      //////////
      // two-body contributions to the forces
      //////////
      auto G1 = compute_2body_fock_deriv<1>(obs, atoms, D);
      for (auto atom = 0, i = 0; atom != atoms.size(); ++atom) {
        for (auto xyz = 0; xyz != 3; ++xyz, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          auto force = G1[i].cwiseProduct(D).sum();
          F2(atom, xyz) += force;
        }
      }

      std::cout << "** 2-body forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) std::cout << F2(atom, xyz) << " ";
      std::cout << std::endl;
#endif

// if support 1-e and 2-e derivatives compute nuclear repulsion force and the
// total force
#if LIBINT2_DERIV_ONEBODY_ORDER && LIBINT2_DERIV_ERI_ORDER
      // compute nuclear repulsion forces
      Matrix FN = Matrix::Zero(atoms.size(), 3);
      //////////
      // nuclear repulsion contribution to the forces
      //////////
      for (auto a1 = 1; a1 != atoms.size(); ++a1) {
        const auto& atom1 = atoms[a1];
        for (auto a2 = 0; a2 < a1; ++a2) {
          const auto& atom2 = atoms[a2];

          auto x12 = atom1.x - atom2.x;
          auto y12 = atom1.y - atom2.y;
          auto z12 = atom1.z - atom2.z;
          auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
          auto r12 = sqrt(r12_2);
          auto r12_3 = r12 * r12_2;

          auto z1z2_over_r12_3 =
              atom1.atomic_number * atom2.atomic_number / r12_3;

          auto fx = -x12 * z1z2_over_r12_3;
          auto fy = -y12 * z1z2_over_r12_3;
          auto fz = -z12 * z1z2_over_r12_3;
          FN(a1, 0) += fx;
          FN(a1, 1) += fy;
          FN(a1, 2) += fz;
          FN(a2, 0) -= fx;
          FN(a2, 1) -= fy;
          FN(a2, 2) -= fz;
        }
      }

      std::cout << "** nuclear repulsion forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) std::cout << FN(atom, xyz) << " ";
      std::cout << std::endl;

      auto F = F1 + F_Pulay + F2 + FN;
      std::cout << "** Hartree-Fock forces = ";
      for (int atom = 0; atom != atoms.size(); ++atom)
        for (int xyz = 0; xyz != 3; ++xyz) std::cout << F(atom, xyz) << " ";
      std::cout << std::endl;
#endif
    }

    {  // compute hessian
      const auto ncoords = 3 * atoms.size();
      // # of elems in upper triangle
      const auto nelem = ncoords * (ncoords + 1) / 2;
#if LIBINT2_DERIV_ONEBODY_ORDER > 1
      // compute 1-e hessian
      Matrix H1 = Matrix::Zero(ncoords, ncoords);
      Matrix H_Pulay = Matrix::Zero(ncoords, ncoords);
      //////////
      // one-body contributions to the hessian
      //////////
      auto T2 = compute_1body_ints_deriv<Operator::kinetic>(2, obs, atoms);
      auto V2 = compute_1body_ints_deriv<Operator::nuclear>(2, obs, atoms);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          auto hess = 2 * (T2[i] + V2[i]).cwiseProduct(D).sum();
          H1(row, col) += hess;
        }
      }

      //////////
      // Pulay hessian
      //////////
      // orbital energy density
      DiagonalMatrix evals_occ(evals.topRows(ndocc));
      Matrix W = C_occ * evals_occ * C_occ.transpose();
      auto S2 = compute_1body_ints_deriv<Operator::overlap>(2, obs, atoms);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          auto hess = 2 * S2[i].cwiseProduct(W).sum();
          H_Pulay(row, col) -= hess;
        }
      }

      std::cout << "** 1-body hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          std::cout << H1(row, col) << " ";
        }
      }
      std::cout << std::endl;

      std::cout << "** Pulay hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          std::cout << H_Pulay(row, col) << " ";
        }
      }
      std::cout << std::endl;
#endif  // LIBINT2_DERIV_ONEBODY_ORDER > 1

#if LIBINT2_DERIV_ERI_ORDER > 1
      // compute 2-e forces
      Matrix H2 = Matrix::Zero(ncoords, ncoords);

      //////////
      // two-body contributions to the forces
      //////////
      auto G2 = compute_2body_fock_deriv<2>(obs, atoms, D);
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col, ++i) {
          // identity prefactor since E(HF) = trace(H + F, D) = trace(2H + G, D)
          auto hess = G2[i].cwiseProduct(D).sum();
          H2(row, col) += hess;
        }
      }

      std::cout << "** 2-body hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          std::cout << H2(row, col) << " ";
        }
      }
      std::cout << std::endl;
#endif

// if support 1-e and 2-e 2nd derivatives compute nuclear repulsion hessian and
// the total hessian
#if LIBINT2_DERIV_ONEBODY_ORDER > 1 && LIBINT2_DERIV_ERI_ORDER > 1
      // compute nuclear repulsion hessian
      // NB only the upper triangle is computed!!!
      Matrix HN = Matrix::Zero(ncoords, ncoords);

      //////////
      // nuclear repulsion contribution to the hessian
      //////////
      for (auto a1 = 1; a1 != atoms.size(); ++a1) {
        const auto& atom1 = atoms[a1];
        for (auto a2 = 0; a2 < a1; ++a2) {
          const auto& atom2 = atoms[a2];

          auto x12 = atom1.x - atom2.x;
          auto y12 = atom1.y - atom2.y;
          auto z12 = atom1.z - atom2.z;
          auto x12_2 = x12 * x12;
          auto y12_2 = y12 * y12;
          auto z12_2 = z12 * z12;
          auto r12_2 = x12 * x12 + y12 * y12 + z12 * z12;
          auto r12 = sqrt(r12_2);
          auto r12_5 = r12 * r12_2 * r12_2;

          auto z1z2_over_r12_5 =
              atom1.atomic_number * atom2.atomic_number / r12_5;

          HN(3 * a1 + 0, 3 * a1 + 0) += z1z2_over_r12_5 * (3 * x12_2 - r12_2);
          HN(3 * a1 + 1, 3 * a1 + 1) += z1z2_over_r12_5 * (3 * y12_2 - r12_2);
          HN(3 * a1 + 2, 3 * a1 + 2) += z1z2_over_r12_5 * (3 * z12_2 - r12_2);
          HN(3 * a1 + 0, 3 * a1 + 1) += z1z2_over_r12_5 * (3 * x12 * y12);
          HN(3 * a1 + 0, 3 * a1 + 2) += z1z2_over_r12_5 * (3 * x12 * z12);
          HN(3 * a1 + 1, 3 * a1 + 2) += z1z2_over_r12_5 * (3 * y12 * z12);

          HN(3 * a2 + 0, 3 * a2 + 0) += z1z2_over_r12_5 * (3 * x12_2 - r12_2);
          HN(3 * a2 + 1, 3 * a2 + 1) += z1z2_over_r12_5 * (3 * y12_2 - r12_2);
          HN(3 * a2 + 2, 3 * a2 + 2) += z1z2_over_r12_5 * (3 * z12_2 - r12_2);
          HN(3 * a2 + 0, 3 * a2 + 1) += z1z2_over_r12_5 * (3 * x12 * y12);
          HN(3 * a2 + 0, 3 * a2 + 2) += z1z2_over_r12_5 * (3 * x12 * z12);
          HN(3 * a2 + 1, 3 * a2 + 2) += z1z2_over_r12_5 * (3 * y12 * z12);

          HN(3 * a2 + 0, 3 * a1 + 0) -= z1z2_over_r12_5 * (3 * x12_2 - r12_2);
          HN(3 * a2 + 1, 3 * a1 + 1) -= z1z2_over_r12_5 * (3 * y12_2 - r12_2);
          HN(3 * a2 + 2, 3 * a1 + 2) -= z1z2_over_r12_5 * (3 * z12_2 - r12_2);
          HN(3 * a2 + 1, 3 * a1 + 0) -= z1z2_over_r12_5 * (3 * y12 * x12);
          HN(3 * a2 + 2, 3 * a1 + 0) -= z1z2_over_r12_5 * (3 * z12 * x12);
          HN(3 * a2 + 2, 3 * a1 + 1) -= z1z2_over_r12_5 * (3 * z12 * y12);
          HN(3 * a2 + 0, 3 * a1 + 1) -= z1z2_over_r12_5 * (3 * x12 * y12);
          HN(3 * a2 + 0, 3 * a1 + 2) -= z1z2_over_r12_5 * (3 * x12 * z12);
          HN(3 * a2 + 1, 3 * a1 + 2) -= z1z2_over_r12_5 * (3 * y12 * z12);
        }
      }

      std::cout << "** nuclear repulsion hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          std::cout << HN(row, col) << " ";
        }
      }
      std::cout << std::endl;

      auto H = H1 + H_Pulay + H2 + HN;
      std::cout << "** Hartree-Fock hessian = ";
      for (auto row = 0, i = 0; row != ncoords; ++row) {
        for (auto col = row; col != ncoords; ++col) {
          std::cout << H(row, col) << " ";
        }
      }
      std::cout << std::endl;
#endif
    }

    libint2::finalize();  // done with libint

  }  // end of try block; if any exceptions occurred, report them and exit
     // cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  } catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  } catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  } catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }

  return 0;
}

std::vector<libint2::Atom> read_geometry(const std::string& filename) {
  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  if (not is.good()) {
    char errmsg[256] = "Could not open file ";
    strncpy(errmsg + 20, filename.c_str(), 235);
    errmsg[255] = '\0';
    throw std::ios_base::failure(errmsg);
  }

  // to prepare for MPI parallelization, we will read the entire file into a
  // string that can be
  // broadcast to everyone, then converted to an std::istringstream object that
  // can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise
  // throw an exception
  if (filename.rfind(".xyz") != std::string::npos)
    return libint2::read_dotxyz(iss);
  else
    throw std::invalid_argument{"only .xyz files are accepted"};
}

// computes Superposition-Of-libint2::Atomic-Densities guess for the molecular
// density matrix in minimal basis; occupies subshells by smearing electrons
// evenly over the orbitals
Matrix compute_soad(const std::vector<libint2::Atom>& atoms) {
  // compute number of atomic orbitals
  size_t nao = 0;
  for (const auto& atom : atoms) {
    const auto Z = atom.atomic_number;
    nao += libint2::sto3g_num_ao(Z);
  }

  // compute the minimal basis density
  Matrix D = Matrix::Zero(nao, nao);
  size_t ao_offset = 0;  // first AO of this atom
  for (const auto& atom : atoms) {
    const auto Z = atom.atomic_number;
    const auto& occvec = libint2::sto3g_ao_occupation_vector(Z);
    for (const auto& occ : occvec) {
      D(ao_offset, ao_offset) = occ;
      ++ao_offset;
    }
  }

  return D * 0.5;  // we use densities normalized to # of electrons/2
}

Matrix compute_shellblock_norm(const BasisSet& obs, const Matrix& A) {
  const auto nsh = obs.size();
  Matrix Ash(nsh, nsh);

  auto shell2bf = obs.shell2bf();
  for (size_t s1 = 0; s1 != nsh; ++s1) {
    const auto& s1_first = shell2bf[s1];
    const auto& s1_size = obs[s1].size();
    for (size_t s2 = 0; s2 != nsh; ++s2) {
      const auto& s2_first = shell2bf[s2];
      const auto& s2_size = obs[s2].size();

      Ash(s1, s2) = A.block(s1_first, s2_first, s1_size, s2_size)
                        .lpNorm<Eigen::Infinity>();
    }
  }

  return Ash;
}

template <Operator obtype, typename OperatorParams>
std::array<Matrix, libint2::operator_traits<obtype>::nopers> compute_1body_ints(
    const BasisSet& obs, OperatorParams params) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  typedef std::array<Matrix, libint2::operator_traits<obtype>::nopers>
      result_type;
  const unsigned int nopers = libint2::operator_traits<obtype>::nopers;
  result_type result;
  for (auto& r : result) r = Matrix::Zero(n, n);

  // construct the 1-body integrals engine
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(obtype, obs.max_nprim(), obs.max_l(), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical
  // charges
  if (obtype == Operator::nuclear || obtype == Operator::erf_nuclear ||
      obtype == Operator::erfc_nuclear || obtype == Operator::erfx_nuclear) {
    engines[0].set_params(params);
  }
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();

  auto compute = [&](int thread_id) {
    const auto& buf = engines[thread_id].results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();
      auto s1_offset = s1 * (s1 + 1) / 2;

      for (auto s2 : obs_shellpair_list[s1]) {
        auto s12 = s1_offset + s2;
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();

        auto n12 = n1 * n2;

        // compute shell pair; return is the pointer to the buffer
        engines[thread_id].compute(obs[s1], obs[s2]);

        for (unsigned int op = 0; op != nopers; ++op) {
          // "map" buffer to a const Eigen Matrix, and copy it to the
          // corresponding blocks of the result
          Eigen::Map<const Matrix> buf_mat(buf[op], n1, n2);
          result[op].block(bf1, bf2, n1, n2) = buf_mat;
          if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding
                         // {s2,s1} block, note the transpose!
            result[op].block(bf2, bf1, n2, n1) = buf_mat.transpose();
        }
      }
    }
  };  // compute lambda

  libint2::parallel_do(compute);

  return result;
}

#if LIBINT2_DERIV_ONEBODY_ORDER
template <Operator obtype>
std::vector<Matrix> compute_1body_ints_deriv(
    unsigned deriv_order, const BasisSet& obs,
    const std::vector<libint2::Atom>& atoms) {
  using libint2::nthreads;
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  constexpr auto nopers = libint2::operator_traits<obtype>::nopers;
  const auto nresults =
      nopers * libint2::num_geometrical_derivatives(atoms.size(), deriv_order);
  typedef std::vector<Matrix> result_type;
  result_type result(nresults);
  for (auto& r : result) r = Matrix::Zero(n, n);

  // construct the 1-body integrals engine
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] =
      libint2::Engine(obtype, obs.max_nprim(), obs.max_l(), deriv_order);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical
  // charges
  if (obtype == Operator::nuclear) {
    std::vector<std::pair<double, std::array<double, 3>>> q;
    for (const auto& atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }
    engines[0].set_params(q);
  }
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

  const auto natoms = atoms.size();
  const auto two_times_ncoords = 6 * natoms;
  const auto nderivcenters_shset =
      2 + ((obtype == Operator::nuclear) ? natoms : 0);

  auto compute = [&](int thread_id) {
    const auto& buf = engines[thread_id].results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();
      auto atom1 = shell2atom[s1];
      assert(atom1 != -1);

      auto s1_offset = s1 * (s1 + 1) / 2;
      for (auto s2 : obs_shellpair_list[s1]) {
        auto s12 = s1_offset + s2;
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = obs[s2].size();
        auto atom2 = shell2atom[s2];

        auto n12 = n1 * n2;

        // compute shell pair; return is the pointer to the buffer
        engines[thread_id].compute(obs[s1], obs[s2]);

        // "copy" lambda copies shell set \c idx to the operator matrix with
        // index \c op
        auto add_shellset_to_dest = [&](std::size_t op, std::size_t idx,
                                        double scale = 1.0) {
          // "map" buffer to a const Eigen Matrix, and copy it to the
          // corresponding blocks of the result
          Eigen::Map<const Matrix> buf_mat(buf[idx], n1, n2);
          if (scale == 1.0)
            result[op].block(bf1, bf2, n1, n2) += buf_mat;
          else
            result[op].block(bf1, bf2, n1, n2) += scale * buf_mat;
          if (s1 != s2) {  // if s1 >= s2, copy {s1,s2} to the corresponding
                           // {s2,s1} block, note the transpose!
            if (scale == 1.0)
              result[op].block(bf2, bf1, n2, n1) += buf_mat.transpose();
            else
              result[op].block(bf2, bf1, n2, n1) += scale * buf_mat.transpose();
          }
        };

        switch (deriv_order) {
          case 0:
            for (std::size_t op = 0; op != nopers; ++op) {
              add_shellset_to_dest(op, op);
            }
            break;

            // map deriv quanta for this shell pair to the overall deriv quanta
            //
            // easiest to explain with example:
            // in sto-3g water shells 0 1 2 sit on atom 0, shells 3 and 4 on
            // atoms 1 and 2 respectively each call to engine::compute for
            // nuclear ints will return derivatives with respect to 15
            // coordinates, obtained as 3 (x,y,z) times 2 + 3 = 5 centers (2
            // centers on which shells sit + 3 nuclear charges) (for overlap,
            // kinetic, and emultipole ints we there are only 6 coordinates
            //  since the operator is coordinate-independent, or derivatives
            //  with respect to the operator coordinates are not computed)
            //

          case 1: {
            std::size_t shellset_idx = 0;
            for (auto c = 0; c != nderivcenters_shset; ++c) {
              auto atom = (c == 0) ? atom1 : ((c == 1) ? atom2 : c - 2);
              auto op_start = 3 * atom * nopers;
              auto op_fence = op_start + nopers;
              for (auto xyz = 0; xyz != 3;
                   ++xyz, op_start += nopers, op_fence += nopers) {
                for (unsigned int op = op_start; op != op_fence;
                     ++op, ++shellset_idx) {
                  add_shellset_to_dest(op, shellset_idx);
                }
              }
            }
          } break;

          case 2: {
            //
            // must pay attention to symmetry when computing 2nd and
            // higher-order derivs e.g. d2 (s1|s2) / dX dY involves several
            // cases:
            // 1. only s1 (or only s2) depends on X AND Y (i.e. X and Y refer to
            // same atom) =>
            //    d2 (s1|s2) / dX dY = (d2 s1 / dX dY | s2)
            // 2. s1 depends on X only, s2 depends on Y only (or vice versa) =>
            //    d2 (s1|s2) / dX dY = (d s1 / dX | d s2 / dY)
            // 3. s1 AND s2 depend on X AND Y (i.e. X and Y refer to same atom)
            // =>
            //    case A: X != Y
            //    d2 (s1|s2) / dX dY = (d2 s1 / dX dY | s2) + (d s1 / dX | d s2
            //    / dY)
            //      + (d s1 / dY | d s2 / dX) + (s1| d2 s2 / dX dY )
            //    case B: X == Y
            //    d2 (s1|s2) / dX2 = (d2 s1 / dX2 | s2) + 2 (d s1 / dX | d s2 /
            //    dX)
            //      + (s1| d2 s2 / dX2 )

            // computes upper triangle index
            // n2 = matrix size times 2
            // i,j = (unordered) indices
            auto upper_triangle_index = [](size_t n2, long i, long j) {
              return std::min(i, j) * (n2 - std::min(i, j) - 1) / 2 +
                     std::max(i, j);
            };

            // look over shellsets in the order in which they appear
            std::size_t shellset_idx = 0;
            for (auto c1 = 0; c1 != nderivcenters_shset; ++c1) {
              auto a1 = (c1 == 0) ? atom1 : ((c1 == 1) ? atom2 : c1 - 2);
              auto coord1 = 3 * a1;
              for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
                for (auto c2 = c1; c2 != nderivcenters_shset; ++c2) {
                  auto a2 = (c2 == 0) ? atom1 : ((c2 == 1) ? atom2 : c2 - 2);
                  auto xyz2_start = (c1 == c2) ? xyz1 : 0;
                  auto coord2 = 3 * a2 + xyz2_start;
                  for (auto xyz2 = xyz2_start; xyz2 != 3; ++xyz2, ++coord2) {
                    double scale = (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

                    const auto coord12 =
                        upper_triangle_index(two_times_ncoords, coord1, coord2);
                    auto op_start = coord12 * nopers;
                    auto op_fence = op_start + nopers;
                    for (auto op = op_start; op != op_fence;
                         ++op, ++shellset_idx) {
                      add_shellset_to_dest(op, shellset_idx, scale);
                    }
                  }
                }
              }
            }
          } break;

          default: {
            assert(false && "not yet implemented");

            using ShellSetDerivIterator =
                libint2::FixedOrderedIntegerPartitionIterator<
                    std::vector<unsigned int>>;
            ShellSetDerivIterator shellset_diter(deriv_order,
                                                 nderivcenters_shset);
            while (shellset_diter) {
              const auto& deriv = *shellset_diter;
            }
          }
        }  // copy shell block switch

      }  // s2 <= s1
    }    // s1
  };     // compute lambda

  libint2::parallel_do(compute);

  return result;
}
#endif

template <libint2::Operator Kernel>
Matrix compute_schwarz_ints(
    const BasisSet& bs1, const BasisSet& _bs2, bool use_2norm,
    typename libint2::operator_traits<Kernel>::oper_params_type params) {
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  Matrix K = Matrix::Zero(nsh1, nsh2);

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  using libint2::nthreads;
  std::vector<Engine> engines(nthreads);

  // !!! very important: cannot screen primitives in Schwarz computation !!!
  auto epsilon = 0.;
  engines[0] = Engine(Kernel, std::max(bs1.max_nprim(), bs2.max_nprim()),
                      std::max(bs1.max_l(), bs2.max_l()), 0, epsilon, params);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::cout << "computing Schwarz bound prerequisites (kernel=" << (int)Kernel
            << ") ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  auto compute = [&](int thread_id) {
    const auto& buf = engines[thread_id].results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      auto n1 = bs1[s1].size();  // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto n2 = bs2[s2].size();
        auto n12 = n1 * n2;

        engines[thread_id].compute2<Kernel, BraKet::xx_xx, 0>(bs1[s1], bs2[s2],
                                                              bs1[s1], bs2[s2]);
        assert(buf[0] != nullptr &&
               "to compute Schwarz ints turn off primitive screening");

        // the diagonal elements are the Schwarz ints ... use Map.diagonal()
        Eigen::Map<const Matrix> buf_mat(buf[0], n12, n12);
        auto norm2 =
            use_2norm ? buf_mat.norm() : buf_mat.lpNorm<Eigen::Infinity>();
        K(s1, s2) = std::sqrt(norm2);
        if (bs1_equiv_bs2) K(s2, s1) = K(s1, s2);
      }
    }
  };  // thread lambda

  libint2::parallel_do(compute);

  timer.stop(0);
  std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return K;
}

Matrix compute_do_ints(const BasisSet& bs1, const BasisSet& bs2,
                       bool use_2norm) {
  return compute_schwarz_ints<libint2::Operator::delta>(bs1, bs2, use_2norm);
}

shellpair_list_t compute_shellpair_list(const BasisSet& bs1,
                                        const BasisSet& _bs2,
                                        const double threshold) {
  const BasisSet& bs2 = (_bs2.empty() ? bs1 : _bs2);
  const auto nsh1 = bs1.size();
  const auto nsh2 = bs2.size();
  const auto bs1_equiv_bs2 = (&bs1 == &bs2);

  using libint2::nthreads;

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  std::vector<Engine> engines;
  engines.reserve(nthreads);
  engines.emplace_back(Operator::overlap,
                       std::max(bs1.max_nprim(), bs2.max_nprim()),
                       std::max(bs1.max_l(), bs2.max_l()), 0);
  for (size_t i = 1; i != nthreads; ++i) {
    engines.push_back(engines[0]);
  }

  std::cout << "computing non-negligible shell-pair list ... ";

  libint2::Timers<1> timer;
  timer.set_now_overhead(25);
  timer.start(0);

  shellpair_list_t result;

  std::mutex mx;

  auto compute = [&](int thread_id) {
    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
      mx.lock();
      if (result.find(s1) == result.end())
        result.insert(std::make_pair(s1, std::vector<size_t>()));
      mx.unlock();

      auto n1 = bs1[s1].size();  // number of basis functions in this shell

      auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
      for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto on_same_center = (bs1[s1].O == bs2[s2].O);
        bool significant = on_same_center;
        if (not on_same_center) {
          auto n2 = bs2[s2].size();
          engines[thread_id].compute(bs1[s1], bs2[s2]);
          Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
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
  };  // end of compute

  libint2::parallel_do(compute);

  // resort shell list in increasing order, i.e. result[s][s1] < result[s][s2]
  // if s1 < s2
  auto sort = [&](int thread_id) {
    for (auto s1 = 0l; s1 != nsh1; ++s1) {
      if (s1 % nthreads == thread_id) {
        auto& list = result[s1];
        std::sort(list.begin(), list.end());
      }
    }
  };  // end of sort

  libint2::parallel_do(sort);

  timer.stop(0);
  std::cout << "done (" << timer.read(0) << " s)" << std::endl;

  return result;
}

// returns {X,X^{-1},rank,A_condition_number,result_A_condition_number}, where
// X is the generalized square-root-inverse such that X.transpose() * A * X = I
//
// if symmetric is true, produce "symmetric" sqrtinv: X = U . A_evals_sqrtinv .
// U.transpose()),
// else produce "canonical" sqrtinv: X = U . A_evals_sqrtinv
// where U are eigenvectors of A
// rows and cols of symmetric X are equivalent; for canonical X the rows are
// original basis (AO),
// cols are transformed basis ("orthogonal" AO)
//
// A is conditioned to max_condition_number
std::tuple<Matrix, Matrix, size_t, double, double> gensqrtinv(
    const Matrix& S, bool symmetric = false,
    double max_condition_number = 1e8) {
  Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(S);
  auto U = eig_solver.eigenvectors();
  auto s = eig_solver.eigenvalues();
  auto s_max = s.maxCoeff();
  auto condition_number = std::min(
      s_max / std::max(s.minCoeff(), std::numeric_limits<double>::min()),
      1.0 / std::numeric_limits<double>::epsilon());
  auto threshold = s_max / max_condition_number;
  long n = s.rows();
  long n_cond = 0;
  for (long i = n - 1; i >= 0; --i) {
    if (s(i) >= threshold) {
      ++n_cond;
    } else
      i = 0;  // skip rest since eigenvalues are in ascending order
  }

  auto sigma = s.bottomRows(n_cond);
  auto result_condition_number = sigma.maxCoeff() / sigma.minCoeff();
  auto sigma_sqrt = sigma.array().sqrt().matrix().asDiagonal();
  auto sigma_invsqrt = sigma.array().sqrt().inverse().matrix().asDiagonal();

  // make canonical X/Xinv
  auto U_cond = U.block(0, n - n_cond, n, n_cond);
  Matrix X = U_cond * sigma_invsqrt;
  Matrix Xinv = U_cond * sigma_sqrt;
  // convert to symmetric, if needed
  if (symmetric) {
    X = X * U_cond.transpose();
    Xinv = Xinv * U_cond.transpose();
  }
  return std::make_tuple(X, Xinv, size_t(n_cond), condition_number,
                         result_condition_number);
}

std::tuple<Matrix, Matrix, double> conditioning_orthogonalizer(
    const Matrix& S, double S_condition_number_threshold) {
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
    std::cout << " (dropped " << obs_nbf_omitted << " "
              << (obs_nbf_omitted > 1 ? "fns" : "fn") << " to reduce to "
              << XtX_condition_number << ")";
  std::cout << std::endl;

  if (obs_nbf_omitted > 0) {
    Matrix should_be_I = X.transpose() * S * X;
    Matrix I = Matrix::Identity(should_be_I.rows(), should_be_I.cols());
    std::cout << "||X^t * S * X - I||_2 = " << (should_be_I - I).norm()
              << " (should be 0)" << std::endl;
  }

  return std::make_tuple(X, Xinv, XtX_condition_number);
}

Matrix compute_2body_2index_ints(const BasisSet& bs) {
  using libint2::nthreads;
  const auto n = bs.nbf();
  const auto nshells = bs.size();
  Matrix result = Matrix::Zero(n, n);

  // build engines for each thread
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] =
      Engine(libint2::Operator::coulomb, bs.max_nprim(), bs.max_l(), 0);
  engines[0].set(BraKet::xs_xs);
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }

  auto shell2bf = bs.shell2bf();
  auto unitshell = Shell::unit();

  auto compute = [&](int thread_id) {
    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0l, s12 = 0l; s1 != nshells; ++s1) {
      auto bf1 = shell2bf[s1];  // first basis function in this shell
      auto n1 = bs[s1].size();

      for (auto s2 = 0; s2 <= s1; ++s2, ++s12) {
        if (s12 % nthreads != thread_id) continue;

        auto bf2 = shell2bf[s2];
        auto n2 = bs[s2].size();

        // compute shell pair; return is the pointer to the buffer
        engine.compute(bs[s1], bs[s2]);
        if (buf[0] == nullptr)
          continue;  // if all integrals screened out, skip to next shell set

        // "map" buffer to a const Eigen Matrix, and copy it to the
        // corresponding blocks of the result
        Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
        result.block(bf1, bf2, n1, n2) = buf_mat;
        if (s1 != s2)  // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                       // block, note the transpose!
          result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
      }
    }
  };  // compute lambda

  libint2::parallel_do(compute);

  return result;
}

Matrix compute_2body_fock(const BasisSet& obs, const Matrix& D,
                          double precision, const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  using libint2::nthreads;
  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
  Matrix D_shblk_norm =
      compute_shellblock_norm(obs, D);  // matrix of infty-norms of shell blocks

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 0);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  std::cout << "compute_2body_fock:precision = " << precision << std::endl;
  std::cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();

  auto lambda = [&](int thread_id) {
    auto& engine = engines[thread_id];
    auto& g = G[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();

          const auto Dnorm123 =
              do_schwarz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dnorm1234 =
                do_schwarz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();

            num_ints_computed += n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            const auto* buf_1234 = buf[0];
            if (buf_1234 == nullptr)
              continue;  // if all integrals screened out, skip to next quartet

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            Eigen::Map<MatrixXd> buf_1234_map(buf_1234, n12, n34);
            assert(buf_1234_map.norm() < Schwarz(s1, s2) * Schwarz(s3, s4));

            // 1) each shell set of integrals contributes up to 6 shell sets of
            // the Fock matrix:
            //    F(a,b) += (ab|cd) * D(c,d)
            //    F(c,d) += (ab|cd) * D(a,b)
            //    F(b,d) -= 1/4 * (ab|cd) * D(a,c)
            //    F(b,c) -= 1/4 * (ab|cd) * D(a,d)
            //    F(a,c) -= 1/4 * (ab|cd) * D(b,d)
            //    F(a,d) -= 1/4 * (ab|cd) * D(b,c)
            // 2) each permutationally-unique integral (shell set) must be
            // scaled by its degeneracy,
            //    i.e. the number of the integrals/sets equivalent to it
            // 3) the end result must be symmetrized
            for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f2 = 0; f2 != n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for (auto f3 = 0; f3 != n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1234[f1234];

                    const auto value_scal_by_deg = value * s1234_deg;

                    g(bf1, bf2) += D(bf3, bf4) * value_scal_by_deg;
                    g(bf3, bf4) += D(bf1, bf2) * value_scal_by_deg;
                    g(bf1, bf3) -= 0.25 * D(bf2, bf4) * value_scal_by_deg;
                    g(bf2, bf4) -= 0.25 * D(bf1, bf3) * value_scal_by_deg;
                    g(bf1, bf4) -= 0.25 * D(bf2, bf3) * value_scal_by_deg;
                    g(bf2, bf3) -= 0.25 * D(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }
  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  Matrix GG = 0.5 * (G[0] + G[0].transpose());

  std::cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  return GG;
}

#if LIBINT2_DERIV_ERI_ORDER
template <unsigned deriv_order>
std::vector<Matrix> compute_2body_fock_deriv(
    const BasisSet& obs, const std::vector<libint2::Atom>& atoms,
    const Matrix& D, double precision, const Matrix& Schwarz) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto nderiv_shellset = libint2::num_geometrical_derivatives(
      4, deriv_order);  // # of derivs for each shell quartet
  const auto nderiv = libint2::num_geometrical_derivatives(
      atoms.size(), deriv_order);  // total # of derivs
  const auto ncoords_times_two = (atoms.size() * 3) * 2;
  using libint2::nthreads;
  std::vector<Matrix> G(nthreads * nderiv, Matrix::Zero(n, n));

  const auto do_schwarz_screen = Schwarz.cols() != 0 && Schwarz.rows() != 0;
  Matrix D_shblk_norm =
      compute_shellblock_norm(obs, D);  // matrix of infty-norms of shell blocks

  auto fock_precision = precision;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  auto max_nprim = obs.max_nprim();
  auto max_nprim4 = max_nprim * max_nprim * max_nprim * max_nprim;
  auto engine_precision = std::min(fock_precision / D_shblk_norm.maxCoeff(),
                                   std::numeric_limits<double>::epsilon()) /
                          max_nprim4;

  // construct the 2-electron repulsion integrals engine pool
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] =
      Engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), deriv_order);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  std::cout << "compute_2body_fock:precision = " << precision << std::endl;
  std::cout << "Engine::precision = " << engines[0].precision() << std::endl;
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  std::atomic<size_t> num_ints_computed{0};

#if defined(REPORT_INTEGRAL_TIMINGS)
  std::vector<libint2::Timers<1>> timers(nthreads);
#endif

  auto shell2bf = obs.shell2bf();
  auto shell2atom = obs.shell2atom(atoms);

  auto lambda = [&](int thread_id) {
    auto& engine = engines[thread_id];
    const auto& buf = engine.results();

#if defined(REPORT_INTEGRAL_TIMINGS)
    auto& timer = timers[thread_id];
    timer.clear();
    timer.set_now_overhead(25);
#endif

    size_t shell_atoms[4];

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell
      shell_atoms[0] = shell2atom[s1];

      for (const auto& s2 : obs_shellpair_list[s1]) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();
        shell_atoms[1] = shell2atom[s2];

        const auto Dnorm12 = do_schwarz_screen ? D_shblk_norm(s1, s2) : 0.;

        for (auto s3 = 0; s3 <= s1; ++s3) {
          auto bf3_first = shell2bf[s3];
          auto n3 = obs[s3].size();
          shell_atoms[2] = shell2atom[s3];

          const auto Dnorm123 =
              do_schwarz_screen
                  ? std::max(D_shblk_norm(s1, s3),
                             std::max(D_shblk_norm(s2, s3), Dnorm12))
                  : 0.;

          const auto s4_max = (s1 == s3) ? s2 : s3;
          for (const auto& s4 : obs_shellpair_list[s3]) {
            if (s4 > s4_max)
              break;  // for each s3, s4 are stored in monotonically increasing
                      // order

            if ((s1234++) % nthreads != thread_id) continue;

            const auto Dnorm1234 =
                do_schwarz_screen
                    ? std::max(
                          D_shblk_norm(s1, s4),
                          std::max(D_shblk_norm(s2, s4),
                                   std::max(D_shblk_norm(s3, s4), Dnorm123)))
                    : 0.;

            if (do_schwarz_screen &&
                Dnorm1234 * Schwarz(s1, s2) * Schwarz(s3, s4) < fock_precision)
              continue;

            auto bf4_first = shell2bf[s4];
            auto n4 = obs[s4].size();
            shell_atoms[3] = shell2atom[s4];

            const auto n1234 = n1 * n2 * n3 * n4;

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;
            auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
            auto s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            auto s1234_deg = s12_deg * s34_deg * s12_34_deg;

            // computes contribution from shell set \c idx to the operator
            // matrix with index \c op
            auto add_shellset_to_dest = [&](std::size_t op, std::size_t idx,
                                            int coord1, int coord2,
                                            double scale = 1.0) {
              auto& g = G[op];
              auto shset = buf[idx];
              const auto weight = scale * s1234_deg;

              for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                const auto bf1 = f1 + bf1_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f3 = 0; f3 != n3; ++f3) {
                    const auto bf3 = f3 + bf3_first;
                    for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                      const auto bf4 = f4 + bf4_first;

                      const auto value = shset[f1234];
                      const auto wvalue = value * weight;

                      g(bf1, bf2) += D(bf3, bf4) * wvalue;
                      g(bf3, bf4) += D(bf1, bf2) * wvalue;
                      g(bf1, bf3) -= 0.25 * D(bf2, bf4) * wvalue;
                      g(bf2, bf4) -= 0.25 * D(bf1, bf3) * wvalue;
                      g(bf1, bf4) -= 0.25 * D(bf2, bf3) * wvalue;
                      g(bf2, bf3) -= 0.25 * D(bf1, bf4) * wvalue;
                    }
                  }
                }
              }
            };

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.start(0);
#endif

            engine.compute2<Operator::coulomb, BraKet::xx_xx, deriv_order>(
                obs[s1], obs[s2], obs[s3], obs[s4]);
            if (buf[0] == nullptr)
              continue;  // if all integrals screened out, skip to next quartet
            num_ints_computed += nderiv_shellset * n1234;

#if defined(REPORT_INTEGRAL_TIMINGS)
            timer.stop(0);
#endif

            switch (deriv_order) {
              case 0: {
                int coord1 = 0, coord2 = 0;
                add_shellset_to_dest(thread_id, 0, coord1, coord2);
              } break;

              case 1: {
                for (auto d = 0; d != 12; ++d) {
                  const int a = d / 3;
                  const int xyz = d % 3;

                  auto coord = shell_atoms[a] * 3 + xyz;
                  auto& g = G[thread_id * nderiv + coord];

                  int coord1 = 0, coord2 = 0;

                  add_shellset_to_dest(thread_id * nderiv + coord, d, coord1,
                                       coord2);

                }  // d \in [0,12)
              } break;

              case 2: {
                // computes upper triangle index
                // n2 = matrix size times 2
                // i,j = (unordered) indices
                auto upper_triangle_index = [](size_t n2, size_t i, size_t j) {
                  return std::min(i, j) * (n2 - std::min(i, j) - 1) / 2 +
                         std::max(i, j);
                };
                // look over shellsets in the order in which they appear
                std::size_t shellset_idx = 0;
                for (auto c1 = 0; c1 != 4; ++c1) {
                  auto a1 = shell_atoms[c1];
                  auto coord1 = 3 * a1;
                  for (auto xyz1 = 0; xyz1 != 3; ++xyz1, ++coord1) {
                    for (auto c2 = c1; c2 != 4; ++c2) {
                      auto a2 = shell_atoms[c2];
                      auto xyz2_start = (c1 == c2) ? xyz1 : 0;
                      auto coord2 = 3 * a2 + xyz2_start;
                      for (auto xyz2 = xyz2_start; xyz2 != 3;
                           ++xyz2, ++coord2) {
                        double scale =
                            (coord1 == coord2 && c1 != c2) ? 2.0 : 1.0;

                        const auto coord12 = upper_triangle_index(
                            ncoords_times_two, coord1, coord2);
                        auto op = thread_id * nderiv + coord12;
                        add_shellset_to_dest(op, shellset_idx, coord1, coord2,
                                             scale);
                        ++shellset_idx;
                      }
                    }
                  }
                }
              } break;

              default:
                assert(deriv_order <= 2 &&
                       "support for 3rd and higher derivatives of the Fock "
                       "matrix not yet implemented");
            }
          }
        }
      }
    }
  };  // end of lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t t = 1; t != nthreads; ++t) {
    for (auto d = 0; d != nderiv; ++d) {
      G[d] += G[t * nderiv + d];
    }
  }

#if defined(REPORT_INTEGRAL_TIMINGS)
  double time_for_ints = 0.0;
  for (auto& t : timers) {
    time_for_ints += t.read(0);
  }
  std::cout << "time for integrals = " << time_for_ints << std::endl;
  for (int t = 0; t != nthreads; ++t) engines[t].print_timers();
#endif

  std::vector<Matrix> GG(nderiv);
  for (auto d = 0; d != nderiv; ++d) {
    GG[d] = 0.5 * (G[d] + G[d].transpose());
  }

  std::cout << "# of integrals = " << num_ints_computed << std::endl;

  // symmetrize the result and return
  return GG;
}

#endif

Matrix compute_2body_fock_general(const BasisSet& obs, const Matrix& D,
                                  const BasisSet& D_bs, bool D_is_shelldiagonal,
                                  double precision) {
  const auto n = obs.nbf();
  const auto nshells = obs.size();
  const auto n_D = D_bs.nbf();
  assert(D.cols() == D.rows() && D.cols() == n_D);

  using libint2::nthreads;
  std::vector<Matrix> G(nthreads, Matrix::Zero(n, n));

  // construct the 2-electron repulsion integrals engine
  using libint2::Engine;
  std::vector<Engine> engines(nthreads);
  engines[0] = Engine(libint2::Operator::coulomb,
                      std::max(obs.max_nprim(), D_bs.max_nprim()),
                      std::max(obs.max_l(), D_bs.max_l()), 0);
  engines[0].set_precision(precision);  // shellset-dependent precision control
                                        // will likely break positive
                                        // definiteness
                                        // stick with this simple recipe
  for (size_t i = 1; i != nthreads; ++i) {
    engines[i] = engines[0];
  }
  auto shell2bf = obs.shell2bf();
  auto shell2bf_D = D_bs.shell2bf();

  auto lambda = [&](int thread_id) {
    auto& engine = engines[thread_id];
    auto& g = G[thread_id];
    const auto& buf = engine.results();

    // loop over permutationally-unique set of shells
    for (auto s1 = 0l, s1234 = 0l; s1 != nshells; ++s1) {
      auto bf1_first = shell2bf[s1];  // first basis function in this shell
      auto n1 = obs[s1].size();       // number of basis functions in this shell

      for (auto s2 = 0; s2 <= s1; ++s2) {
        auto bf2_first = shell2bf[s2];
        auto n2 = obs[s2].size();

        for (auto s3 = 0; s3 < D_bs.size(); ++s3) {
          auto bf3_first = shell2bf_D[s3];
          auto n3 = D_bs[s3].size();

          auto s4_begin = D_is_shelldiagonal ? s3 : 0;
          auto s4_fence = D_is_shelldiagonal ? s3 + 1 : D_bs.size();

          for (auto s4 = s4_begin; s4 != s4_fence; ++s4, ++s1234) {
            if (s1234 % nthreads != thread_id) continue;

            auto bf4_first = shell2bf_D[s4];
            auto n4 = D_bs[s4].size();

            // compute the permutational degeneracy (i.e. # of equivalents) of
            // the given shell set
            auto s12_deg = (s1 == s2) ? 1.0 : 2.0;

            if (s3 >= s4) {
              auto s34_deg = (s3 == s4) ? 1.0 : 2.0;
              auto s1234_deg = s12_deg * s34_deg;
              // auto s1234_deg = s12_deg;
              engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                  obs[s1], obs[s2], D_bs[s3], D_bs[s4]);
              const auto* buf_1234 = buf[0];
              if (buf_1234 != nullptr) {
                for (auto f1 = 0, f1234 = 0; f1 != n1; ++f1) {
                  const auto bf1 = f1 + bf1_first;
                  for (auto f2 = 0; f2 != n2; ++f2) {
                    const auto bf2 = f2 + bf2_first;
                    for (auto f3 = 0; f3 != n3; ++f3) {
                      const auto bf3 = f3 + bf3_first;
                      for (auto f4 = 0; f4 != n4; ++f4, ++f1234) {
                        const auto bf4 = f4 + bf4_first;

                        const auto value = buf_1234[f1234];
                        const auto value_scal_by_deg = value * s1234_deg;
                        g(bf1, bf2) += 2.0 * D(bf3, bf4) * value_scal_by_deg;
                      }
                    }
                  }
                }
              }
            }

            engine.compute2<Operator::coulomb, BraKet::xx_xx, 0>(
                obs[s1], D_bs[s3], obs[s2], D_bs[s4]);
            const auto* buf_1324 = buf[0];
            if (buf_1324 == nullptr)
              continue;  // if all integrals screened out, skip to next quartet

            for (auto f1 = 0, f1324 = 0; f1 != n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for (auto f3 = 0; f3 != n3; ++f3) {
                const auto bf3 = f3 + bf3_first;
                for (auto f2 = 0; f2 != n2; ++f2) {
                  const auto bf2 = f2 + bf2_first;
                  for (auto f4 = 0; f4 != n4; ++f4, ++f1324) {
                    const auto bf4 = f4 + bf4_first;

                    const auto value = buf_1324[f1324];
                    const auto value_scal_by_deg = value * s12_deg;
                    g(bf1, bf2) -= D(bf3, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }
  };  // thread lambda

  libint2::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i != nthreads; ++i) {
    G[0] += G[i];
  }

  // symmetrize the result and return
  return 0.5 * (G[0] + G[0].transpose());
}

#ifdef HAVE_DENSITY_FITTING

Matrix DFFockEngine::compute_2body_fock_dfC(const Matrix& Cocc) {
  using libint2::nthreads;

  const auto n = obs.nbf();
  const auto ndf = dfbs.nbf();

  libint2::Timers<1> wall_timer;
  wall_timer.set_now_overhead(25);
  std::vector<libint2::Timers<5>> timers(nthreads);
  for (auto& timer : timers) timer.set_now_overhead(25);

  typedef btas::RangeNd<CblasRowMajor, std::array<long, 1>> Range1d;
  typedef btas::RangeNd<CblasRowMajor, std::array<long, 2>> Range2d;
  typedef btas::Tensor<double, Range1d> Tensor1d;
  typedef btas::Tensor<double, Range2d> Tensor2d;

  // using first time? compute 3-center ints and transform to inv sqrt
  // representation
  if (xyK.size() == 0) {
    wall_timer.start(0);

    const auto nshells = obs.size();
    const auto nshells_df = dfbs.size();
    const auto& unitshell = libint2::Shell::unit();

    // construct the 2-electron 3-center repulsion integrals engine
    // since the code assumes (xx|xs) braket, and Engine/libint only produces
    // (xs|xx), use 4-center engine
    std::vector<libint2::Engine> engines(nthreads);
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
                                 std::max(obs.max_nprim(), dfbs.max_nprim()),
                                 std::max(obs.max_l(), dfbs.max_l()), 0);
    engines[0].set(BraKet::xs_xx);
    for (size_t i = 1; i != nthreads; ++i) {
      engines[i] = engines[0];
    }

    auto shell2bf = obs.shell2bf();
    auto shell2bf_df = dfbs.shell2bf();

    Tensor3d Zxy{ndf, n, n};

    auto lambda = [&](int thread_id) {
      auto& engine = engines[thread_id];
      auto& timer = timers[thread_id];
      const auto& results = engine.results();

      // loop over permutationally-unique set of shells
      for (auto s1 = 0l, s123 = 0l; s1 != nshells_df; ++s1) {
        auto bf1_first = shell2bf_df[s1];  // first basis function in this shell
        auto n1 = dfbs[s1].size();  // number of basis functions in this shell

        for (auto s2 = 0; s2 != nshells; ++s2) {
          auto bf2_first = shell2bf[s2];
          auto n2 = obs[s2].size();
          const auto n12 = n1 * n2;

          for (auto s3 = 0; s3 != nshells; ++s3, ++s123) {
            if (s123 % nthreads != thread_id) continue;

            auto bf3_first = shell2bf[s3];
            auto n3 = obs[s3].size();
            const auto n123 = n12 * n3;

            timer.start(0);

            engine.compute2<Operator::coulomb, BraKet::xs_xx, 0>(
                dfbs[s1], unitshell, obs[s2], obs[s3]);
            const auto* buf = results[0];
            if (buf == nullptr) continue;

            timer.stop(0);
            timer.start(1);

            auto lower_bound = {bf1_first, bf2_first, bf3_first};
            auto upper_bound = {bf1_first + n1, bf2_first + n2, bf3_first + n3};
            auto view = btas::make_view(
                Zxy.range().slice(lower_bound, upper_bound), Zxy.storage());
            std::copy(buf, buf + n123, view.begin());

            timer.stop(1);
          }  // s3
        }    // s2
      }      // s1
    };       // lambda

    libint2::parallel_do(lambda);

    wall_timer.stop(0);

    double ints_time = 0;
    for (const auto& timer : timers) ints_time += timer.read(0);
    std::cout << "time for Zxy integrals = " << ints_time
              << " (total from all threads)" << std::endl;
    double copy_time = 0;
    for (const auto& timer : timers) copy_time += timer.read(1);
    std::cout << "time for copying into BTAS = " << copy_time
              << " (total from all threads)" << std::endl;
    std::cout << "wall time for Zxy integrals + copy = " << wall_timer.read(0)
              << std::endl;

    timers[0].start(2);

    Matrix V = compute_2body_2index_ints(dfbs);
    Eigen::LLT<Matrix> V_LLt(V);
    Matrix I = Matrix::Identity(ndf, ndf);
    auto L = V_LLt.matrixL();
    Matrix V_L = L;
    Matrix Linv = L.solve(I).transpose();
    // check
    //  std::cout << "||V - L L^t|| = " << (V - V_L * V_L.transpose()).norm() <<
    //  std::endl;
    //  std::cout << "||I - L L^-1^t|| = " << (I - V_L *
    //  Linv.transpose()).norm() << std::endl;
    //  std::cout << "||V^-1 - L^-1 L^-1^t|| = " << (V.inverse() - Linv *
    //  Linv.transpose()).norm() << std::endl;

    Tensor2d K{ndf, ndf};
    std::copy(Linv.data(), Linv.data() + ndf * ndf, K.begin());

    xyK = Tensor3d{n, n, ndf};
    btas::contract(1.0, Zxy, {1, 2, 3}, K, {1, 4}, 0.0, xyK, {2, 3, 4});
    Zxy = Tensor3d{0, 0, 0};  // release memory

    timers[0].stop(2);
    std::cout << "time for integrals metric tform = " << timers[0].read(2)
              << std::endl;
  }  // if (xyK.size() == 0)

  // compute exchange
  timers[0].start(3);

  const auto nocc = Cocc.cols();
  Tensor2d Co{n, nocc};
  std::copy(Cocc.data(), Cocc.data() + n * nocc, Co.begin());
  Tensor3d xiK{n, nocc, ndf};
  btas::contract(1.0, xyK, {1, 2, 3}, Co, {2, 4}, 0.0, xiK, {1, 4, 3});

  Tensor2d G{n, n};
  btas::contract(1.0, xiK, {1, 2, 3}, xiK, {4, 2, 3}, 0.0, G, {1, 4});

  timers[0].stop(3);
  std::cout << "time for exchange = " << timers[0].read(3) << std::endl;

  // compute Coulomb
  timers[0].start(4);

  Tensor1d Jtmp{ndf};
  btas::contract(1.0, xiK, {1, 2, 3}, Co, {1, 2}, 0.0, Jtmp, {3});
  xiK = Tensor3d{0, 0, 0};
  btas::contract(2.0, xyK, {1, 2, 3}, Jtmp, {3}, -1.0, G, {1, 2});

  timers[0].stop(4);
  std::cout << "time for coulomb = " << timers[0].read(4) << std::endl;

  // copy result to an Eigen::Matrix
  Matrix result(n, n);
  std::copy(G.cbegin(), G.cend(), result.data());
  return result;
}
#endif  // HAVE_DENSITY_FITTING

// should be a unit test somewhere
void api_basic_compile_test(const BasisSet& obs) {
  using namespace libint2;
  Engine onebody_engine(
      Operator::overlap,  // will compute overlap ints
      obs.max_nprim(),    // max # of primitives in shells this engine will
                          // accept
      obs.max_l()  // max angular momentum of shells this engine will accept
  );
  auto shell2bf = obs.shell2bf();
  const auto& results = onebody_engine.results();
  for (auto s1 = 0; s1 != obs.size(); ++s1) {
    for (auto s2 = 0; s2 != obs.size(); ++s2) {
      std::cout << "compute shell set {" << s1 << "," << s2 << "} ... ";
      onebody_engine.compute(obs[s1], obs[s2]);
      const auto* ints_shellset = results[0];
      std::cout << "done" << std::endl;

      auto bf1 = shell2bf[s1];   // first basis function in first shell
      auto n1 = obs[s1].size();  // number of basis functions in first shell
      auto bf2 = shell2bf[s2];   // first basis function in second shell
      auto n2 = obs[s2].size();  // number of basis functions in second shell

      // this iterates over integrals in the order they are packed in array
      // ints_shellset
      for (auto f1 = 0; f1 != n1; ++f1)
        for (auto f2 = 0; f2 != n2; ++f2)
          std::cout << "  " << bf1 + f1 << " " << bf2 + f2 << " "
                    << ints_shellset[f1 * n2 + f2] << std::endl;
    }
  }

  using libint2::Operator;

  std::vector<std::pair<double, double>> cgtg_params{
      {0.1, 0.2}, {0.3, 0.4}, {0.5, 0.6}};
  {
    auto K = compute_schwarz_ints<Operator::cgtg>(obs, obs, false, cgtg_params);
    std::cout << "cGTG Schwarz ints\n" << K << std::endl;
  }
  {
    auto K = compute_schwarz_ints<Operator::cgtg_x_coulomb>(obs, obs, false,
                                                            cgtg_params);
    std::cout << "cGTG/r12 Schwarz ints\n" << K << std::endl;
  }
  {
    auto K =
        compute_schwarz_ints<Operator::delcgtg2>(obs, obs, false, cgtg_params);
    std::cout << "||Del.cGTG||^2 Schwarz ints\n" << K << std::endl;
  }

  {  // test 2-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());
    Engine eri2_engine = eri4_engine;
    eri2_engine.set(BraKet::xs_xs);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results2 = eri2_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());
        eri2_engine.compute(obs[s1], obs[s2]);

        auto bf1 = shell2bf[s1];   // first basis function in first shell
        auto n1 = obs[s1].size();  // number of basis functions in first shell
        auto bf2 = shell2bf[s2];   // first basis function in second shell
        auto n2 = obs[s2].size();  // number of basis functions in second shell

        const auto* buf4 = results4[0];
        const auto* buf2 = results2[0];

        // this iterates over integrals in the order they are packed in array
        // ints_shellset
        for (auto f1 = 0, f12 = 0; f1 != n1; ++f1)
          for (auto f2 = 0; f2 != n2; ++f2, ++f12)
            assert(std::abs(buf4[f12] - buf2[f12]) < 1e-12 &&
                   "2-center ints test failed");
      }
    }
  }
  {  // test 3-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l());
    Engine eri3_engine = eri4_engine;
    eri3_engine.set(BraKet::xs_xx);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results3 = eri3_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        for (auto s3 = 0; s3 != obs.size(); ++s3) {
          eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], obs[s3]);
          eri3_engine.compute(obs[s1], obs[s2], obs[s3]);

          auto bf1 = shell2bf[s1];   // first basis function in first shell
          auto n1 = obs[s1].size();  // number of basis functions in first shell
          auto bf2 = shell2bf[s2];   // first basis function in second shell
          auto n2 =
              obs[s2].size();       // number of basis functions in second shell
          auto bf3 = shell2bf[s3];  // first basis function in third shell
          auto n3 = obs[s3].size();  // number of basis functions in third shell

          const auto* buf4 = results4[0];
          const auto* buf3 = results3[0];

          // this iterates over integrals in the order they are packed in array
          // ints_shellset
          for (auto f1 = 0, f123 = 0; f1 != n1; ++f1)
            for (auto f2 = 0; f2 != n2; ++f2)
              for (auto f3 = 0; f3 != n3; ++f3, ++f123)
                assert(std::abs(buf4[f123] - buf3[f123]) < 1e-12 &&
                       "3-center ints test failed");
        }
      }
    }
  }

#if LIBINT2_DERIV_ERI_ORDER
  {  // test deriv 2-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 1);
    Engine eri2_engine = eri4_engine;
    eri2_engine.set(BraKet::xs_xs);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results2 = eri2_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());
        eri2_engine.compute(obs[s1], obs[s2]);

        auto bf1 = shell2bf[s1];   // first basis function in first shell
        auto n1 = obs[s1].size();  // number of basis functions in first shell
        auto bf2 = shell2bf[s2];   // first basis function in second shell
        auto n2 = obs[s2].size();  // number of basis functions in second shell

        // loop over derivative shell sets
        for (auto d = 0; d != 6; ++d) {
          const auto* buf4 = results4[d < 3 ? d : d + 3];
          const auto* buf2 = results2[d];

          // this iterates over integrals in the order they are packed in array
          // ints_shellset
          for (auto f1 = 0, f12 = 0; f1 != n1; ++f1)
            for (auto f2 = 0; f2 != n2; ++f2, ++f12)
              assert(std::abs(buf4[f12] - buf2[f12]) < 1e-12 &&
                     "deriv 2-center ints test failed");
        }
      }
    }
  }
#endif

#if LIBINT2_DERIV_ERI_ORDER > 1
  {  // test 2nd deriv 2-index ints
    Engine eri4_engine(Operator::coulomb, obs.max_nprim(), obs.max_l(), 2);
    Engine eri2_engine = eri4_engine;
    eri2_engine.set(BraKet::xs_xs);
    auto shell2bf = obs.shell2bf();
    const auto& results4 = eri4_engine.results();
    const auto& results2 = eri2_engine.results();
    for (auto s1 = 0; s1 != obs.size(); ++s1) {
      for (auto s2 = 0; s2 != obs.size(); ++s2) {
        eri4_engine.compute(obs[s1], Shell::unit(), obs[s2], Shell::unit());
        eri2_engine.compute(obs[s1], obs[s2]);

        auto bf1 = shell2bf[s1];   // first basis function in first shell
        auto n1 = obs[s1].size();  // number of basis functions in first shell
        auto bf2 = shell2bf[s2];   // first basis function in second shell
        auto n2 = obs[s2].size();  // number of basis functions in second shell

        // loop over derivative shell sets
        for (auto d1 = 0, d12 = 0; d1 != 6; ++d1) {
          const auto dd1 = d1 < 3 ? d1 : d1 + 3;
          for (auto d2 = d1; d2 != 6; ++d2, ++d12) {
            const auto dd2 = d2 < 3 ? d2 : d2 + 3;
            const auto dd12 = dd1 * (24 - dd1 - 1) / 2 + dd2;
            const auto* buf4 = results4[dd12];
            const auto* buf2 = results2[d12];

            // this iterates over integrals in the order they are packed in
            // array
            // ints_shellset
            for (auto f1 = 0, f12 = 0; f1 != n1; ++f1)
              for (auto f2 = 0; f2 != n2; ++f2, ++f12)
                assert(std::abs(buf4[f12] - buf2[f12]) < 1e-12 &&
                       "2nd deriv 2-center ints test failed");
          }
        }
      }
    }
  }
#endif
}
