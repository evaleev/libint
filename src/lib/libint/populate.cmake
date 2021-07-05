configure_file("${PROJECT_SOURCE_DIR}/export/LICENSE.export" "${EXPORT_STAGE_DIR}/LICENSE" COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/export/INSTALL.export"
        "${EXPORT_STAGE_DIR}/INSTALL" COPYONLY)

file(INSTALL "${PROJECT_SOURCE_DIR}/COPYING"
        "${PROJECT_SOURCE_DIR}/COPYING.LESSER"
        "${PROJECT_SOURCE_DIR}/CITATION"
        "${PROJECT_SOURCE_DIR}/README.md"
        DESTINATION "${EXPORT_STAGE_DIR}")

file(INSTALL "${PROJECT_SOURCE_DIR}/doc/progman/progman.tex"
        "${PROJECT_SOURCE_DIR}/doc/progman/refs.bib"
        "${PROJECT_SOURCE_DIR}/doc/Libint_Logo3_alt.pdf"
        "${PROJECT_SOURCE_DIR}/doc/Libint_Logo3_alt.eps"
        DESTINATION "${EXPORT_STAGE_DIR}/doc")

file(INSTALL "${PROJECT_SOURCE_DIR}/doc/progman/" DESTINATION "${EXPORT_STAGE_DIR}/doc"
        FILES_MATCHING PATTERN "*.cc")

file(INSTALL "${PROJECT_SOURCE_DIR}/tests/"
        DESTINATION "${EXPORT_STAGE_DIR}/tests"
        FILES_MATCHING PATTERN "*.c"
        PATTERN "*.cc"
        PATTERN "*.h"
        PATTERN "*.hpp"
        PATTERN "*.py"
        PATTERN "*.xyz"
        PATTERN "CMakeLists.txt")

file(INSTALL "${PROJECT_SOURCE_DIR}/src/bin/test_eri/eri.h"
             "${PROJECT_SOURCE_DIR}/src/bin/test_eri/prep_libint2.h"
        DESTINATION "${EXPORT_STAGE_DIR}/tests/eri")

# Q: why is this different from current export location?
file(INSTALL "${PROJECT_SOURCE_DIR}/lib/basis/"
        DESTINATION "${EXPORT_STAGE_DIR}/lib/basis"
        FILES_MATCHING PATTERN "*.g94")

file(INSTALL "${PROJECT_SOURCE_DIR}/external/boost.tar.gz"
        DESTINATION "${EXPORT_STAGE_DIR}/external")

file(INSTALL "${PROJECT_SOURCE_DIR}/cmake/modules/autocmake_safeguards.cmake"
        "${PROJECT_SOURCE_DIR}/cmake/modules/int_orderings.cmake"
        "${PROJECT_SOURCE_DIR}/cmake/modules/options.cmake"
        "${PROJECT_SOURCE_DIR}/cmake/modules/FindMPFR.cmake"
        "${PROJECT_SOURCE_DIR}/cmake/modules/AddCustomTargetSubproject.cmake"
        "${PROJECT_BINARY_DIR}/cmake/modules/int_computed.cmake"
        DESTINATION "${EXPORT_STAGE_DIR}/cmake/modules")

file(INSTALL "${PROJECT_SOURCE_DIR}/cmake/Libint2Config.cmake.in"
        "${PROJECT_SOURCE_DIR}/cmake/libint2.pc.cmake.in"
        DESTINATION "${EXPORT_STAGE_DIR}/cmake")

configure_file("${PROJECT_SOURCE_DIR}/src/lib/libint/CMakeLists.txt.export" "${EXPORT_STAGE_DIR}/CMakeLists.txt" COPYONLY)

file(INSTALL "${PROJECT_SOURCE_DIR}/include/libint2.h"
        "${PROJECT_SOURCE_DIR}/include/libint2.hpp"
        DESTINATION "${EXPORT_STAGE_DIR}/include")

file(INSTALL "${PROJECT_SOURCE_DIR}/include/libint2/"
        DESTINATION "${EXPORT_STAGE_DIR}/include/libint2"
        FILES_MATCHING PATTERN "*.h")

file(INSTALL "${PROJECT_SOURCE_DIR}/src/bin/libint/util_types.h"
        DESTINATION "${EXPORT_STAGE_DIR}/include/libint2")

file(INSTALL "${PROJECT_SOURCE_DIR}/src/lib/libint/"
        DESTINATION "${EXPORT_STAGE_DIR}/include/libint2"
        FILES_MATCHING PATTERN "*.h")


