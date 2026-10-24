// Independently generated Chebyshev data; LGPL-3.0-or-later.
// Regenerate with tools/generate_truncated_coulomb_table.cc.
#ifndef LIBINT2_TRUNCATED_COULOMB_TABLE_H
#define LIBINT2_TRUNCATED_COULOMB_TABLE_H
namespace libint2 { namespace detail {
static const int truncated_coulomb_max_m=64;
static const int truncated_coulomb_degree=16;
static const unsigned truncated_coulomb_block_size=524288;
static const long double truncated_coulomb_width=20.3479128454902600254L;
struct TruncatedCoulombNode {double split,a,b,c,d; int left,right,axis,poly;};
struct TruncatedCoulombPolynomial {unsigned offset; unsigned char nx,ny[17];};
extern const TruncatedCoulombNode truncated_coulomb_nodes[];
extern const TruncatedCoulombPolynomial truncated_coulomb_polynomials[];
extern const double* const truncated_coulomb_coefficients[];
}}
#endif
