set(LIBRARY_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/lib/libint)

file(
  INSTALL
    "${LIBRARY_SOURCE_DIR}/../../../tests/"
  DESTINATION "${EXPORT_STAGE_DIR}/tests"
  FILES_MATCHING
    PATTERN "*.c"
    PATTERN "*.cc"
    PATTERN "*.h"
    PATTERN "*.hpp"
    PATTERN "*.py"
    PATTERN "*.xyz"
    PATTERN "CMakeLists.txt"
  )

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/src/bin/test_eri/eri.h"
    "${PROJECT_SOURCE_DIR}/src/bin/test_eri/prep_libint2.h"
  DESTINATION "${EXPORT_STAGE_DIR}/tests/eri"
  )

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

configure_file(
  "${LIBRARY_SOURCE_DIR}/CMakeLists.txt.export"
  "${EXPORT_STAGE_DIR}/CMakeLists.txt"
  COPYONLY)


# <<<  Headers  >>>

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/include/libint2.h"
    "${PROJECT_SOURCE_DIR}/include/libint2.hpp"
  DESTINATION "${EXPORT_STAGE_DIR}/include"
  )
file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/include/libint2/"
  DESTINATION "${EXPORT_STAGE_DIR}/include/libint2"
  FILES_MATCHING
    PATTERN "*.h"
    PATTERN "*.h.cmake.in"
  )
file(
  INSTALL
    "${LIBRARY_SOURCE_DIR}/"
  DESTINATION "${EXPORT_STAGE_DIR}/include/libint2"
  FILES_MATCHING
    PATTERN "*.h"
  )


# <<<  Source  >>>

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/src/lib/libint/configuration.cc.cmake.in"
  DESTINATION "${EXPORT_STAGE_DIR}/src"
  )

file(GLOB generated_sources_list RELATIVE "${EXPORT_STAGE_DIR}" "${EXPORT_STAGE_DIR}/src/*.cc")
file(WRITE ${EXPORT_STAGE_DIR}/srclist.cmake "set(LIBINT2_LIBRARY_CXX_SRC \"${generated_sources_list}\" )")
