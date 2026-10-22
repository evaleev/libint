// Independent offline Chebyshev-table generator. LGPL-3.0-or-later.
// Build this optional tool with a compiler supported by the installed Boost.
#include <algorithm>
#include <array>
#include <atomic>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <queue>
#include <string>
#include <thread>
#include <utility>
#include <vector>
using Big =
    boost::multiprecision::number<boost::multiprecision::cpp_bin_float<100>>;
int order = 64, degree = 16, n = degree + 1, jobs = 8;
const long double tol = 2.220446049250313080847263336181640625e-16L / 32;
const long double rscale = 8;
long double width;

// Reference arithmetic is confined to this offline generator. The independent
// mpmath fixtures additionally validate the generated runtime interpolator.
std::vector<long double> reference(int region, long double u, long double v) {
  std::vector<long double> result(order + 1);
  Big R, T;
  if (region == 0) {
    R = Big(rscale) * Big(u);
    Big s = (Big(rscale) + Big(width)) * Big(v);
    T = s * s;
  } else {
    const Big d = Big(width) * (2 * Big(v) - 1);
    if (u == 0) {
      result[0] = static_cast<long double>(sqrt(acos(Big(-1))) / 4 * erfc(d));
      return result;
    }
    R = Big(rscale) / Big(u);
    Big s = R + d;
    T = s * s;
  }
  if (R == 0) return result;
  const Big factor = region ? R : Big(1);
  if (T <= order / 2 + 4) {
    Big w = 1;
    int k = 0;
    do {
      ++k;
      w *= T / k;
    } while (k < 2 * T + 2 || abs(w) > Big(tol) / 1000);
    const int last = order + k;
    std::vector<Big> g(last + 1);
    const Big x = R * R, e = exp(-x);
    Big prev = 0, cur = x < Big(.5) ? Big(-expm1(-x)) : e;
    const bool small = x < Big(.5);
    for (int m = 0; m <= last; ++m) {
      g[m] = (small ? cur : Big(1) - cur) / (2 * m + 1);
      Big next =
          ((2 * m + Big(.5) - x) * cur - m * prev + (small ? x : Big(0))) /
          (m + Big(.5));
      prev = cur;
      cur = next;
    }
    for (int m = 0; m <= order; ++m) {
      Big sum = g[m], w = 1;
      for (int j = 1; j <= k; ++j) {
        w *= -T / j;
        sum += w * g[m + j];
      }
      result[m] = static_cast<long double>(sum * factor);
    }
  } else {
    const Big s = sqrt(T), x = R * R, et = exp(-T);
    const Big a = exp(-(R - s) * (R - s)), b = exp(-(R + s) * (R + s));
    Big h;
    if (R <= s)
      h = erfc(s - R) - 2 * erfc(s) + erfc(s + R);
    else
      h = 2 * erf(s) - erfc(R - s) + erfc(R + s);
    Big g = sqrt(acos(Big(-1))) * h / (4 * s), prev = 0, cur = (a + b) / 2,
        next = cur - R / s * (a - b) / 2;
    for (int m = 0; m <= order; ++m) {
      result[m] = static_cast<long double>(g * factor);
      g = ((2 * m + 1) * g - et + cur) / (2 * T);
      const Big after = ((2 * T + m + Big(.5)) * next -
                         (T + 2 * m + Big(.5) - x) * cur + m * prev) /
                        T;
      prev = cur;
      cur = next;
      next = after;
    }
  }
  return result;
}
struct Cell {
  int region;
  long double a, b, c, d;
  int left = -1, right = -1, axis = -1;
  long double split = 0;
};
struct Fit {
  bool accept = false;
  int axis = 0;
  long double error = 0;
  std::vector<long double> coeff;
};
std::vector<long double> cosine;
long double evaluate(const long double* p, long double x, long double y) {
  long double rows[65];
  for (int i = 0; i < n; ++i) {
    long double b1 = 0, b2 = 0;
    for (int j = n - 1; j > 0; --j) {
      auto b = 2 * y * b1 - b2 + p[i * n + j];
      b2 = b1;
      b1 = b;
    }
    rows[i] = y * b1 - b2 + p[i * n];
  }
  long double b1 = 0, b2 = 0;
  for (int i = n - 1; i > 0; --i) {
    auto b = 2 * x * b1 - b2 + rows[i];
    b2 = b1;
    b1 = b;
  }
  return x * b1 - b2 + rows[0];
}
Fit fit(const Cell& cell) {
  Fit out;
  out.coeff.resize((order + 1) * n * n);
  std::vector<long double> samples((order + 1) * n * n), temp(n * n);
  const auto pi = acosl(-1.L);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      const auto x = cosl(pi * (i + .5L) / n), y = cosl(pi * (j + .5L) / n);
      const auto u = (cell.a + cell.b + (cell.b - cell.a) * x) / 2;
      const auto v = (cell.c + cell.d + (cell.d - cell.c) * y) / 2;
      const auto values = reference(cell.region, u, v);
      for (int m = 0; m <= order; ++m) samples[(m * n + i) * n + j] = values[m];
    }
  long double tx = 0, ty = 0;
  for (int m = 0; m <= order; ++m) {
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j) {
        long double sum = 0;
        for (int i = 0; i < n; ++i)
          sum += samples[(m * n + i) * n + j] * cosine[k * n + i];
        temp[k * n + j] = sum * (k == 0 ? 1.L / n : 2.L / n);
      }
    long double ex = 0, ey = 0;
    for (int k = 0; k < n; ++k)
      for (int l = 0; l < n; ++l) {
        long double sum = 0;
        for (int j = 0; j < n; ++j) sum += temp[k * n + j] * cosine[l * n + j];
        sum *= l == 0 ? 1.L / n : 2.L / n;
        out.coeff[(m * n + k) * n + l] = sum;
        if (k >= n - 3) ex += std::abs(sum);
        if (l >= n - 3) ey += std::abs(sum);
      }
    tx = std::max(tx, ex);
    ty = std::max(ty, ey);
  }
  out.axis = tx >= ty ? 0 : 1;
  // Staggered verification grid, including all edges and corners.
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 9; ++j) {
      const auto x = cosl(pi * i / 8), y = cosl(pi * j / 8);
      const auto u = (cell.a + cell.b + (cell.b - cell.a) * x) / 2;
      const auto v = (cell.c + cell.d + (cell.d - cell.c) * y) / 2;
      const auto values = reference(cell.region, u, v);
      const auto scale = cell.region ? u / rscale : 1;
      for (int m = 0; m <= order; ++m)
        out.error = std::max(
            out.error,
            std::abs(evaluate(&out.coeff[m * n * n], x, y) - values[m]) *
                scale);
    }
  out.accept = out.error <= tol;
  return out;
}
int main(int argc, char** argv) {
  if (argc < 2 || argc > 5) {
    std::cerr
        << "usage: generator output-prefix [mmax=64] [degree=16] [jobs=8]\n";
    return 1;
  }
  if (argc > 2) order = std::stoi(argv[2]);
  if (argc > 3) degree = std::stoi(argv[3]);
  if (argc > 4) jobs = std::stoi(argv[4]);
  if (order < 0 || degree < 4 || degree > 64 || jobs < 1) return 1;
  n = degree + 1;
  cosine.resize(n * n);
  width = std::sqrt(4 * (order + std::log(1.L / tol)));
  const auto pi = acosl(-1.L);
  for (int k = 0; k < n; ++k)
    for (int i = 0; i < n; ++i)
      cosine[k * n + i] = cosl(pi * k * (i + .5L) / n);
  std::vector<Cell> cells{{0, 0, 1, 0, 1}, {1, 0, 1, 0, 1}};
  std::vector<int> pending{0, 1};
  std::vector<std::vector<long double>> coefficients(2);
  int level = 0, accepted = 0;
  while (!pending.empty()) {
    if (level++ > 24) {
      std::cerr << "too many refinement levels\n";
      return 2;
    }
    std::vector<Fit> fits(pending.size());
    std::atomic<std::size_t> next{0};
    std::vector<std::thread> workers;
    for (int k = 0; k < jobs; ++k)
      workers.emplace_back([&]() {
        for (;;) {
          auto i = next.fetch_add(1);
          if (i >= pending.size()) break;
          fits[i] = fit(cells[pending[i]]);
        }
      });
    for (auto& t : workers) t.join();
    std::vector<int> following;
    long double worst = 0;
    for (std::size_t i = 0; i < pending.size(); ++i) {
      const int index = pending[i];
      auto& f = fits[i];
      worst = std::max(worst, f.error);
      if (f.accept) {
        coefficients[index] = std::move(f.coeff);
        ++accepted;
        continue;
      }
      Cell a = cells[index], b = a;
      const auto midpoint = f.axis == 0 ? (a.a + a.b) / 2 : (a.c + a.d) / 2;
      if (f.axis == 0) {
        a.b = midpoint;
        b.a = midpoint;
      } else {
        a.d = midpoint;
        b.c = midpoint;
      }
      cells[index].axis = f.axis;
      cells[index].split = midpoint;
      cells[index].left = cells.size();
      cells[index].right = cells.size() + 1;
      following.push_back(cells.size());
      cells.push_back(a);
      coefficients.emplace_back();
      following.push_back(cells.size());
      cells.push_back(b);
      coefficients.emplace_back();
    }
    std::cerr << "level " << level << " tested " << pending.size()
              << " accepted-total " << accepted << " worst " << double(worst)
              << '\n';
    pending.swap(following);
  }
  const unsigned block_size = 1u << 19;
  struct Polynomial {
    unsigned offset;
    std::vector<unsigned> lengths;
  };
  std::vector<Polynomial> descriptors;
  std::vector<double> data;
  std::vector<int> polynomial_index(cells.size(), -1);
  for (std::size_t k = 0; k < cells.size(); ++k) {
    const auto& cell = cells[k];
    if (coefficients[k].empty()) continue;
    polynomial_index[k] = descriptors.size();
    const long double budget = 2 * tol / (cell.region ? cell.b / rscale : 1);
    // A per-polynomial L1 discard budget bounds error on [-1,1]^2.
    for (int m = 0; m <= order; ++m) {
      std::vector<double> p(n * n);
      for (int j = 0; j < n * n; ++j)
        p[j] = static_cast<double>(coefficients[k][m * n * n + j]);
      std::vector<unsigned> lengths(n, n);
      using Entry = std::pair<double, unsigned>;
      std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> tails;
      for (int i = 0; i < n; ++i)
        tails.push({std::abs(p[i * n + n - 1]), unsigned(i)});
      long double discarded = 0;
      while (!tails.empty()) {
        const auto entry = tails.top();
        tails.pop();
        if (discarded + entry.first > budget) break;
        discarded += entry.first;
        auto& length = lengths[entry.second];
        --length;
        if (length)
          tails.push(
              {std::abs(p[entry.second * n + length - 1]), entry.second});
      }
      while (!lengths.empty() && lengths.back() == 0) lengths.pop_back();
      unsigned count = 0;
      for (auto length : lengths) count += length;
      if (data.size() % block_size + count > block_size)
        data.resize((data.size() / block_size + 1) * block_size, 0.);
      descriptors.push_back({unsigned(data.size()), lengths});
      for (unsigned i = 0; i < lengths.size(); ++i)
        for (unsigned j = 0; j < lengths[i]; ++j) data.push_back(p[i * n + j]);
    }
  }
  const unsigned blocks = (data.size() + block_size - 1) / block_size;
  const std::string prefix = argv[1];
  std::ofstream h(prefix + ".h"), cc(prefix + ".cc");
  h << "// Independently generated Chebyshev data; LGPL-3.0-or-later.\n"
    << "// Regenerate with tools/generate_truncated_coulomb_table.cc.\n"
    << "#ifndef LIBINT2_TRUNCATED_COULOMB_TABLE_H\n#define "
       "LIBINT2_TRUNCATED_COULOMB_TABLE_H\n"
    << "namespace libint2 { namespace detail {\n"
    << "static const int truncated_coulomb_max_m=" << order << ";\n"
    << "static const int truncated_coulomb_degree=" << degree << ";\n"
    << "static const unsigned truncated_coulomb_block_size=" << block_size
    << ";\n"
    << "static const long double truncated_coulomb_width="
    << std::setprecision(21) << width << "L;\n"
    << "struct TruncatedCoulombNode {double split,a,b,c,d; int "
       "left,right,axis,poly;};\n"
    << "struct TruncatedCoulombPolynomial {unsigned offset; unsigned char "
       "nx,ny["
    << n << "];};\n"
    << "extern const TruncatedCoulombNode truncated_coulomb_nodes[];\n"
    << "extern const TruncatedCoulombPolynomial "
       "truncated_coulomb_polynomials[];\n"
    << "extern const double* const "
       "truncated_coulomb_coefficients[];\n}}\n#endif\n";
  cc << "// Independent generated Chebyshev data. LGPL-3.0-or-later.\n"
     << "#include <libint2/truncated_coulomb_table.h>\nnamespace libint2 { "
        "namespace detail {\n"
     << "const TruncatedCoulombNode truncated_coulomb_nodes[]={\n"
     << std::setprecision(17);
  for (std::size_t k = 0; k < cells.size(); ++k) {
    const auto& c = cells[k];
    cc << '{' << double(c.split) << ',' << double(c.a) << ',' << double(c.b)
       << ',' << double(c.c) << ',' << double(c.d) << ',' << c.left << ','
       << c.right << ',' << c.axis << ',' << polynomial_index[k] << "},\n";
  }
  cc << "};\nconst TruncatedCoulombPolynomial "
        "truncated_coulomb_polynomials[]={\n";
  for (const auto& d : descriptors) {
    cc << '{' << d.offset << ',' << d.lengths.size() << ",{";
    for (auto length : d.lengths) cc << length << ',';
    cc << "}},\n";
  }
  cc << "};\n";
  for (unsigned k = 0; k < blocks; ++k)
    cc << "extern const double truncated_coulomb_block_" << k << "[];\n";
  cc << "const double* const truncated_coulomb_coefficients[]={\n";
  for (unsigned k = 0; k < blocks; ++k)
    cc << "truncated_coulomb_block_" << k << ",\n";
  cc << "};\n}}\n";
  for (unsigned k = 0; k < blocks; ++k) {
    std::ofstream block(prefix + "_" + std::to_string(k) + ".cc");
    block << "// Independently generated coefficients. LGPL-3.0-or-later.\n"
          << "namespace libint2 { namespace detail {\nextern const double "
             "truncated_coulomb_block_"
          << k << "[]={\n"
          << std::setprecision(17);
    for (std::size_t j = k * block_size;
         j < std::min(data.size(), std::size_t(k + 1) * block_size); ++j)
      block << data[j] << ((j % 4 == 3) ? ",\n" : ",");
    block << "\n};\n}}\n";
  }
  std::cerr << "leaves=" << accepted << " coefficients=" << data.size()
            << " bytes=" << data.size() * 8 << '\n';
}
