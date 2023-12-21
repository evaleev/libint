list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/modules)

include(int_computed)  # for macros.tex and features

add_custom_target(libint-library-generate DEPENDS "CMakeLists.txt.export") # TODO "${EXPORT_STAGE_DIR}/src/libint2_params.h")

# <<<  Add Metadata To The Library Source  >>>

configure_file(
  "features.cmake.in"
  "${EXPORT_STAGE_DIR}/features"
  @ONLY)

add_custom_command(
  OUTPUT
    ${EXPORT_STAGE_DIR}/CMakeLists.txt
  DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake"
    libint-library-generate
  COMMAND ${CMAKE_COMMAND}
    "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
    "-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}"
    "-DEXPORT_STAGE_DIR=${EXPORT_STAGE_DIR}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake
  COMMENT "Populating Libint2 library"
  )
add_custom_target(libint-library-populate DEPENDS ${EXPORT_STAGE_DIR}/CMakeLists.txt)


# <<<  Export The Library Source  >>>

add_custom_command(
  OUTPUT
    "${EXPORT_STAGE_DIR}.tgz"
  COMMAND ${CMAKE_COMMAND} -E tar "cfz" "${EXPORT_STAGE_DIR}.tgz" "${EXPORT_STAGE_DIR}"
  DEPENDS
    libint-library-generate
    libint-library-populate
  WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/.."
  COMMENT "Exporting tarball of Libint2 library source"
  )
add_custom_target(libint-library-export DEPENDS "${EXPORT_STAGE_DIR}.tgz")


    # <<<  Build Library  >>>

    include(ExternalProject)
    ExternalProject_Add(
      library
      DEPENDS
        libint-library-export
      BUILD_ALWAYS 1
      # download generalized to generating source by unpacking source tarball
      # * DOWNLOAD_DIR <SOURCE_DIR> is tidier but fails on Windows
      URL "${EXPORT_STAGE_DIR}.tgz"
      LOG_DOWNLOAD 1
      CMAKE_CACHE_ARGS
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
        -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
        -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
      )

