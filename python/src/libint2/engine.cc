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

#include <libint2.hpp>

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#include <atomic>
#include <new>
#include <thread>
#include <tuple>
#include <vector>

#if defined(_MSC_VER)
#include <BaseTsd.h>
// handles ssize_t in pybind11/numpy.h
typedef SSIZE_T ssize_t;
#endif

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "basis.h"

namespace py = pybind11;

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

namespace libint2::python::engine {

static std::optional<size_t> num_threads_;

size_t num_threads() {
  if (!num_threads_) {
    return std::thread::hardware_concurrency();
  }
  return *num_threads_;
}

void set_num_threads(size_t num_threads) { num_threads_ = num_threads; }

std::vector<size_t> leading_dimensions(size_t bytes,
                                       std::vector<size_t> sizes) {
  std::vector<size_t> dims;
  size_t n = 1;
  for (auto it = sizes.rbegin(); it != sizes.rend(); ++it) {
    dims.push_back(bytes * n);
    n *= *it;
  }
  std::reverse(dims.begin(), dims.end());
  return dims;
}

template <typename... Shells>
py::object compute(libint2::Engine &engine, const Shells &...shells) {
  using real_type = libint2::value_type;
  engine.compute((const libint2::Shell &)shells...);
  const auto &buf = engine.results();
  // py::print("engine::compute:", buf[0], shells...);
  if (!buf[0]) return py::none();
  // return buf[0][0];
  py::buffer_info buffer(
      const_cast<real_type *>(buf[0]), sizeof(real_type),
      py::format_descriptor<real_type>::format(), sizeof...(Shells),
      {shells.size()...},
      leading_dimensions(sizeof(real_type), {shells.size()...}));
  return py::array_t<real_type>(buffer);
}

void parallel_for(std::function<void(std::atomic<int> &)> task,
                  size_t num_threads) {
  assert(num_threads);
  std::atomic<int> counter(0);
  std::vector<std::thread> threads;
  for (size_t k = 0; k < num_threads; ++k) {
    threads.emplace_back([&, task]() { task(counter); });
  }
  for (auto &thread : threads) {
    thread.join();
  }
}

template <class Task>
void parallel_for_symmetric(Task &&task, size_t num_threads, size_t n) {
  auto thread_task = [&, task](std::atomic<int> &counter) mutable {
    size_t next = counter++;
    for (size_t i = 0, ij = 0; i < n; ++i) {
      for (size_t j = 0; j <= i; ++j, ++ij) {
        if (ij != next) continue;
        next = counter++;
        task(i, j);
      }
    }
  };
  parallel_for(thread_task, num_threads);
}

template <class Task, class... Args>
void parallel_for_eval(Task &&task, std::tuple<Args...> args) {
  std::apply(task, args);
}

template <class Task, class... Args, class Range, class... Ranges>
void parallel_for_eval(Task &&task, std::tuple<Args...> args, Range &&R,
                       Ranges &&...Rs) {
  for (auto r : R) {
    parallel_for_eval(task, std::tuple_cat(args, std::forward_as_tuple(r)),
                      Rs...);
  }
}

template <class Task, class Range, class... CD>
void parallel_for(Task &&task, int num_threads, const Range &A, const Range &B,
                  CD &&...cd) {
  auto thread_task = [&, task](std::atomic<int> &counter) mutable {
    size_t next = counter++;
    size_t idx = 0;
    for (auto a : A) {
      for (auto b : B) {
        if (idx++ != next) continue;
        next = counter++;
        parallel_for_eval(task, std::forward_as_tuple(a, b), cd...);
      }
    }
  };
  parallel_for(thread_task, num_threads);
}

template <class... Args>
py::object parallel_compute(libint2::Engine &engine, const Args &...args) {
  constexpr auto N = sizeof...(Args);

  // Create Python array and grab a Eigen::Tensor view into its data
  const std::array<size_t, N> dims = {nbf(args)...};
  py::array_t<double> result(dims);
  auto result_ptr = result.mutable_data();
  Eigen::TensorMap<Eigen::Tensor<double, N, Eigen::RowMajor>> V(result_ptr,
                                                                dims);

  auto task = [&, engine](auto... braket) mutable {
    // Compute the integral - this does not use Python API
    engine.compute(std::get<0>(braket)...);
    const auto &buf = engine.results();
    if (!buf[0]) return;

    // Get shell sizes and basis function offsets
    std::array<size_t, N> shell_sizes = {std::get<0>(braket).size()...};
    std::array<size_t, N> bf_offsets = {
        static_cast<size_t>(std::get<1>(braket))...};

    auto V_sh = V.slice(bf_offsets, shell_sizes);
    // N.B. Libint returns shells sets in row-major layout
    V_sh = Eigen::TensorMap<Eigen::Tensor<const double, N, Eigen::RowMajor>>(
        buf[0], shell_sizes);
  };

  parallel_for(task, num_threads(), basis::enumerate(args)...);

  return result;
}

template <class Arg>
py::object compute2(libint2::Engine &engine, const Arg &, const Arg &);

template <class Arg>
py::object compute3(libint2::Engine &engine, const Arg &, const Arg &,
                    const Arg &);

template <class Arg>
py::object compute4(libint2::Engine &engine, const Arg &, const Arg &,
                    const Arg &, const Arg &);

template <>
py::object compute2(libint2::Engine &engine, const Shell &A, const Shell &B) {
  return compute(engine, A, B);
}

template <>
py::object compute3(libint2::Engine &engine, const Shell &A, const Shell &B,
                    const Shell &C) {
  return compute(engine, A, B, C);
}

template <>
py::object compute4(libint2::Engine &engine, const Shell &A, const Shell &B,
                    const Shell &C, const Shell &D) {
  return compute(engine, A, B, C, D);
}

template <>
py::object compute2(libint2::Engine &engine, const BasisSet &A,
                    const BasisSet &B) {
  return parallel_compute(engine, A, B);
}

template <>
py::object compute3(libint2::Engine &engine, const BasisSet &A,
                    const BasisSet &B, const BasisSet &C) {
  return parallel_compute(engine, A, B, C);
}

template <>
py::object compute4(libint2::Engine &engine, const BasisSet &A,
                    const BasisSet &B, const BasisSet &C, const BasisSet &D) {
  return parallel_compute(engine, A, B, C, D);
}

Matrix compute_1body_ints(libint2::Engine &engine,
                          const libint2::BasisSet &basis) {
  Matrix h = Matrix::Zero(basis.nbf(), basis.nbf());
  auto task = [&, engine](size_t i, size_t j) mutable {
    const auto &P = basis[i];
    const auto &Q = basis[j];
    Matrix h_ij = Matrix::Zero(P.size(), Q.size());
    engine.compute(P, Q);
    const auto *buf = engine.results()[0];
    if (!buf) return;
    double sym = 1.0;
    if (i == j) sym /= 2.0;
    for (size_t p = 0; p < P.size(); ++p) {
      for (size_t q = 0; q < Q.size(); ++q) {
        auto buf_val = sym * (*buf++);
        h_ij(p, q) = buf_val;
      }
    }
    const auto &bf = basis.shell2bf();
    h.block(bf[i], bf[j], h_ij.rows(), h_ij.cols()) = h_ij;
  };
  parallel_for_symmetric(task, num_threads(), basis.size());
  return h + h.transpose();
}

Matrix compute_2body_fock(libint2::Engine &engine,
                          const Eigen::Ref<const Matrix> &D,
                          const libint2::BasisSet &basis) {
  size_t nbf = basis.nbf();
  Matrix f = Matrix::Zero(nbf, nbf);

  const auto &shell = basis;
  const auto &bf = basis.shell2bf();

  std::mutex mutex;

  auto task = [&, engine](size_t P, size_t Q) mutable {
    Matrix f_i = Matrix::Zero(shell[P].size(), nbf);
    Matrix f_j = Matrix::Zero(shell[Q].size(), nbf);
    // notice r/s symmetry not implemented
    for (size_t R = 0; R < basis.size(); ++R) {
      for (size_t S = 0; S <= R; ++S) {
        engine.compute(shell[P], shell[Q], shell[R], shell[S]);
        const auto *buf = engine.results()[0];
        if (!buf) continue;
        double sym = 1.0;
        if (P == Q) sym /= 2.0;
        if (R == S) sym /= 2.0;
        for (size_t p = 0; p < shell[P].size(); ++p) {
          for (size_t q = 0; q < shell[Q].size(); ++q) {
            for (size_t r = 0; r < shell[R].size(); ++r) {
              for (size_t s = 0; s < shell[S].size(); ++s) {
                double v = sym * (*buf++);
                f_i(p, bf[Q] + q) += 4 * v * D(bf[R] + r, bf[S] + s);
                f_i(p, bf[S] + s) -= v / 2 * D(bf[Q] + q, bf[R] + r);
                f_i(p, bf[R] + r) -= v / 2 * D(bf[Q] + q, bf[S] + s);
                f_j(q, bf[R] + r) -= v / 2 * D(bf[P] + p, bf[S] + s);
                f_j(q, bf[S] + s) -= v / 2 * D(bf[P] + p, bf[R] + r);
              }
            }
          }
        }
      }
    }

    {
      std::unique_lock<std::mutex> lock(mutex);
      f.block(bf[P], 0, f_i.rows(), f_i.cols()) += f_i;
      f.block(bf[Q], 0, f_j.rows(), f_j.cols()) += f_j;
    }
  };

  parallel_for_symmetric(task, num_threads(), basis.size());

  return f + f.transpose();
}

}  // namespace libint2::python::engine
