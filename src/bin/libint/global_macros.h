
#ifndef _libint2_src_bin_libint_globalmacros_h_
#define _libint2_src_bin_libint_globalmacros_h_

#define DEBUG 0

/// Controls whether classes derived from GenIntegralSet overload its label()
#define OVERLOAD_GENINTEGRALSET_LABEL 0

/// Controls whether DGVertex's typeid_ is used to quickly screen vertex types and thus avoid dynamic cast
#define PTREQUIV_USE_TYPEID 1

/** Controls whether PtrEquiv avoids using operator== to compare vertices in favor of a comparing keys associated with each object --
    thus avoiding static cast */
#define PTREQUIV_USE_KEY_TO_COMPARE 1

/// If 1 then DGVertex's instid_ is used to as a key, otherwise label() is used as a key
#define PTREQUIV_USE_INSTID 1

/** Controls whether AlgebraicOperator avoids using operator== to compare vertices in favor of a comparing keys associated with each object --
    thus avoiding static cast + comparison code */
#define ALGEBRAICOPERATOR_USE_KEY_TO_COMPARE 1

/** If 1 then AlgebraicOperator compares left and right arguments directly (comparing pointers), else DGVertex::equiv() is used */
#define ALGEBRAICOPERATOR_USE_SAFEPTR 0

#endif
