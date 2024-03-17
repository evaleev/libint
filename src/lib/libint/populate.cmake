set(LIBRARY_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/lib/libint)

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/cmake/modules/int_orderings.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/options.cmake"
    "${PROJECT_BINARY_DIR}/cmake/modules/int_computed.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/modules/int_userreal.cmake"
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

file(
  INSTALL
    "${PROJECT_SOURCE_DIR}/src/lib/libint/configuration.cc.cmake.in"
  DESTINATION "${EXPORT_STAGE_DIR}/src"
  )

