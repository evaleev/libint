set(LIBRARY_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/lib/libint)

# ====  /doc  ===================================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/doc/progman/progman.tex"
    "${PROJECT_SOURCE_DIR}/doc/progman/refs.bib"
    "${PROJECT_SOURCE_DIR}/doc/Libint_Logo3_alt.pdf"
    "${PROJECT_SOURCE_DIR}/doc/Libint_Logo3_alt.eps"
  DESTINATION "${EXPORT_STAGE_DIR}/doc"
  )

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/doc/progman/"
  DESTINATION "${EXPORT_STAGE_DIR}/doc"
  FILES_MATCHING
    PATTERN "*.cc"
  )

# ====  /tests  =================================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/tests/"
  DESTINATION "${EXPORT_STAGE_DIR}/tests"
  FILES_MATCHING
    PATTERN "*.c"
    PATTERN "*.cc"
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*.F90"
    PATTERN "*.py"
    PATTERN "*.xyz"
    PATTERN "hftest.cmake"
    PATTERN "CMakeLists.txt"
    # PATTERN "ssss.nb"
  )

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/export/fortran/fortran_example.F90"
  DESTINATION "${EXPORT_STAGE_DIR}/tests/fortran"
  )

file(READ "${PROJECT_SOURCE_DIR}/export/fortran/test.cc" _file_contents)
string(REPLACE "tests/unit" "unit" _file_contents "${_file_contents}")
file(WRITE "${EXPORT_STAGE_DIR}/tests/fortran/test.cc" "${_file_contents}")

file(READ "${PROJECT_SOURCE_DIR}/export/fortran/test-eri.cc" _file_contents)
string(REPLACE "tests/unit" "unit" _file_contents "${_file_contents}")
string(REPLACE "tests/eri" "eri" _file_contents "${_file_contents}")
file(WRITE "${EXPORT_STAGE_DIR}/tests/fortran/test-eri.cc" "${_file_contents}")

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/src/bin/test_eri/eri.h"
    "${PROJECT_SOURCE_DIR}/src/bin/test_eri/prep_libint2.h"
  DESTINATION "${EXPORT_STAGE_DIR}/tests/eri"
  )

# ====  /lib/basis  =============================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/lib/basis/"
  DESTINATION "${EXPORT_STAGE_DIR}/lib/basis"
  FILES_MATCHING
    PATTERN "*.g94"
  )

# file(INSTALL) preserves symlinks, and tar -xf fails on them on Windows (either
#   the symlink extracted before the target file or problems with the stars in
#   filenames themselves). To make export tarballs broadly usable, we'll copy
#   the symlinked files into real files, which is what libtool effectively did.
file(
  COPY_FILE
    "${PROJECT_SOURCE_DIR}/lib/basis/6-31gs.g94"
    "${EXPORT_STAGE_DIR}/lib/basis/6-31gs.g94"
  )
file(
  COPY_FILE
    "${PROJECT_SOURCE_DIR}/lib/basis/6-31gss.g94"
    "${EXPORT_STAGE_DIR}/lib/basis/6-31gss.g94"
  )
file(
  COPY_FILE
    "${PROJECT_SOURCE_DIR}/lib/basis/6-311gss.g94"
    "${EXPORT_STAGE_DIR}/lib/basis/6-311gss.g94"
  )

# ====  /fortran  ===============================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/export/fortran/c_to_f.py"
    "${PROJECT_SOURCE_DIR}/export/fortran/make_defs.py"
    "${PROJECT_SOURCE_DIR}/export/fortran/libint_f.F90"
  DESTINATION "${EXPORT_STAGE_DIR}/fortran"
  )

# ====  /python  ================================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/python/"
  DESTINATION "${EXPORT_STAGE_DIR}/python"
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.cc"
    PATTERN "*.py"
    PATTERN "*.py.in"
    PATTERN "CMakeLists.txt"
  )

# ====  /external  ==============================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/external/boost.tar.gz"
  DESTINATION "${EXPORT_STAGE_DIR}/external"
  )

# ====  /cmake  =================================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/cmake/modules/autocmake_safeguards.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/int_orderings.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/options.cmake"
    "${PROJECT_BINARY_DIR}/cmake/modules/int_computed.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/int_userreal.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/int_checkboost.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/FindMultiprecision.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/FindEigen3.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/AddCustomTargetSubproject.cmake"
  DESTINATION "${EXPORT_STAGE_DIR}/cmake/modules"
  )

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/cmake/libint2-config.cmake.in"
  DESTINATION "${EXPORT_STAGE_DIR}/cmake"
  )

# ====  /include  ===============================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/include/libint2.h"
    "${PROJECT_SOURCE_DIR}/include/libint2.hpp"
  DESTINATION "${EXPORT_STAGE_DIR}/include"
  )

# Note that libint2_iface.h, libint2_params.h, libint2_types.h shift around.
#   They're generated along with the integrals .h/.cc library src, then get
#   exported to include/ (along with the integrals .h), then are finally
#   installed (not with the integrals .h) into include/libint2/ . The
#   __COMPILING_LIBINT define and the include/libint2/util/generated/libint2_*.h
#   redirection headers take care of the "build tree"/"export" setup.
# In a cmake+cmake buildsystem, one could probably install these three headers
#   to both locations and forego the define.
file(
  INSTALL
    "${PROJECT_BINARY_DIR}/generated/"
  DESTINATION "${EXPORT_STAGE_DIR}/include"
  FILES_MATCHING
    PATTERN "*.h"
  )

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/src/bin/libint/util_types.h"
  DESTINATION "${EXPORT_STAGE_DIR}/include"
  )

file(
  INSTALL
    "${LIBRARY_SOURCE_DIR}/"
  DESTINATION "${EXPORT_STAGE_DIR}/include"
  FILES_MATCHING
    PATTERN "*.h"
  )

# ====  /include/libint2  =======================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/include/libint2/"
  DESTINATION "${EXPORT_STAGE_DIR}/include/libint2"
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.h.cmake.in"
    PATTERN "basis.h.in"  # TODO basis.h.cmake.in after libtool retires
  )

# ====  /src  ===================================================================

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/src/lib/libint/engine.cpp"
    "${PROJECT_SOURCE_DIR}/src/lib/libint/configuration.cc.cmake.in"
  DESTINATION "${EXPORT_STAGE_DIR}/src"
  )

file(
  INSTALL
    "${PROJECT_BINARY_DIR}/generated/"
  DESTINATION "${EXPORT_STAGE_DIR}/src"
  FILES_MATCHING
    PATTERN "*.cc"
  )

file(GLOB generated_sources_list RELATIVE "${EXPORT_STAGE_DIR}" "${EXPORT_STAGE_DIR}/src/*.cc")
file(WRITE ${EXPORT_STAGE_DIR}/srclist.cmake "set(LIBINT2_LIBRARY_CXX_SRC \"${generated_sources_list}\" )")

# ====  /  ======================================================================

configure_file(
  "${PROJECT_SOURCE_DIR}/export/LICENSE.export"
  "${EXPORT_STAGE_DIR}/LICENSE"
  COPYONLY)

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/INSTALL"
    "${PROJECT_SOURCE_DIR}/INSTALL.md"
    "${PROJECT_SOURCE_DIR}/COPYING"
    "${PROJECT_SOURCE_DIR}/COPYING.LESSER"
    "${PROJECT_SOURCE_DIR}/CITATION"
    "${PROJECT_SOURCE_DIR}/README.md"
  DESTINATION "${EXPORT_STAGE_DIR}"
  )

configure_file(
  "${LIBRARY_SOURCE_DIR}/CMakeLists.txt.export"
  "${EXPORT_STAGE_DIR}/CMakeLists.txt"
  COPYONLY)
