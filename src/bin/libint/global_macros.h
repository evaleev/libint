/*
 *  Copyright (C) 2004-2026 Edward F. Valeev
 *
 *  This file is part of Libint compiler.
 *
 *  Libint compiler is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Libint compiler is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Libint compiler.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _libint2_src_bin_libint_globalmacros_h_
#define _libint2_src_bin_libint_globalmacros_h_

/// Feel free to adjust higher, but not yet tested
#define LIBINT_CARTGAUSS_MAX_AM 32

/// For developers only
#define LIBINT_SUPPORT_ONEBODYINTS 1

/// Set to 1 to enable various safety checks which are normally too expensive to
/// perform
#define CHECK_SAFETY 0

/// DirectedGraph uses an associated container (multimap) to store vertices? If
/// not, use a simple container (list)
#define USE_ASSOCCONTAINER_BASED_DIRECTEDGRAPH 1

/// Controls whether classes derived from GenIntegralSet overload its label()
#define OVERLOAD_GENINTEGRALSET_LABEL 0

/// Controls whether DGVertex's typeid_ is used to quickly screen vertex types
/// and thus avoid dynamic cast
#define PTREQUIV_USE_TYPEID 1

/** Controls whether PtrEquiv avoids using operator== to compare vertices in
   favor of a comparing keys associated with each object -- thus avoiding static
   cast */
#define PTREQUIV_USE_KEY_TO_COMPARE 1

/// If 1 then DGVertex's instid_ is used to as a key, otherwise label() is used
/// as a key
#define PTREQUIV_USE_INSTID 1

/** Controls whether AlgebraicOperator avoids using operator== to compare
   vertices in favor of a comparing keys associated with each object -- thus
   avoiding static cast + comparison code */
#define ALGEBRAICOPERATOR_USE_KEY_TO_COMPARE 1

/** If 1 then AlgebraicOperator compares left and right arguments directly
 * (comparing pointers), else DGVertex::equiv() is used */
#define ALGEBRAICOPERATOR_USE_SHAREDPTR 0

/// Use VectorBraket from braket.h
#define USE_BRAKET_H 1

/// Use integer key to hash integrals, rather than string label -- this should
/// reduce memory consumption
#define USE_INT_KEY_TO_HASH 1

/// Use integer key to compare
#define USE_INT_KEY_TO_COMPARE 1
/// GenIntegralSet use unique integer keys to hash integrals and avoid creating
/// temporaries
#if !USE_INT_KEY_TO_COMPARE
#error "For now USE_INT_KEY_TO_COMPARE must be 1"
#endif

/// If set to 1 then avoid using SubIterator to compute size
#define COMPUTE_SIZE_DIRECTLY 1

/// If set to 0 then complex expressions will be condensed into single-line
/// expressions, which should help linewise vectorization
#define DISABLE_SUBTREES 0

/// Set to 1 to produce GraphViz-formatted DAGs
#define PRINT_DAG_GRAPHVIZ 0

/// Produce massive amounts of debugging info
#define DEBUG 0
#define DEBUG_RESTRUCTURE 0
#define DEBUG_TRAVERSAL 0
#define DEBUG_CONSTRUCTION 0

#endif
