list(APPEND CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/cmake/modules)

include(GNUInstallDirs)
include(int_computed)  # for macros.tex and features

# <<<  Generate Library Source  >>>

file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/generated")
add_custom_command(
  OUTPUT
    "${PROJECT_BINARY_DIR}/generated/libint2_params.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.cc"
  COMMAND $<TARGET_FILE:build_libint>
  DEPENDS
    build_libint
  WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/generated"
  COMMENT "Generating Libint2 library source"
  )
add_custom_target(libint-library-generate DEPENDS "${PROJECT_BINARY_DIR}/generated/libint2_params.h")


# <<<  Add Metadata To The Library Source  >>>

# the library source to be tarballed is mainly assembled through populate.cmake.
#   These following configure_file's need a few extra substitution variables that
#   aren't passed to populate.

configure_file(
  "features.cmake.in"
  "${EXPORT_STAGE_DIR}/features"
  @ONLY)

configure_file(
  "${PROJECT_SOURCE_DIR}/doc/progman/macros.tex.in"
  "${EXPORT_STAGE_DIR}/doc/macros.tex"
  @ONLY)

# position basis.h for local buildtree.
# * gets re-configured (CMakeLists.txt.export) for library export with its own
#   CMAKE_INSTALL_PREFIX
# * gets re-re-positioned (end of this file) for case of continuous generator +
#   library build (through EP) with outer CMAKE_INSTALL_PREFIX
set(LIBINT_DATADIR_ABSOLUTE ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION})
configure_file(
  "${PROJECT_SOURCE_DIR}/include/libint2/basis.h.in"  # TODO basis.h.cmake.in after libtool retires
  "include/libint2/basis.h"
  @ONLY)

# populate.cmake copies whole directory trees, so glob them for its dependencies;
#   CONFIGURE_DEPENDS re-globs at build time so that added files are picked up too.
#   Dotfiles (e.g. .DS_Store) are skipped to avoid spurious re-exports.
file(GLOB_RECURSE _populate_inputs CONFIGURE_DEPENDS
  "${PROJECT_SOURCE_DIR}/export/*"
  "${PROJECT_SOURCE_DIR}/include/libint2/*"
  "${PROJECT_SOURCE_DIR}/doc/progman/*"
  "${PROJECT_SOURCE_DIR}/cmake/modules/*.cmake"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
  )
list(FILTER _populate_inputs EXCLUDE REGEX "/\\.[^/]*$")
list(APPEND _populate_inputs
  "${PROJECT_SOURCE_DIR}/include/libint2.h"
  "${PROJECT_SOURCE_DIR}/include/libint2.hpp"
  "${PROJECT_SOURCE_DIR}/src/bin/libint/util_types.h"
  "${PROJECT_SOURCE_DIR}/src/bin/test_eri/eri.h"
  "${PROJECT_SOURCE_DIR}/src/bin/test_eri/prep_libint2.h"
  "${PROJECT_SOURCE_DIR}/external/boost.tar.gz"
  "${PROJECT_SOURCE_DIR}/doc/Libint_Logo3_alt.pdf"
  "${PROJECT_SOURCE_DIR}/doc/Libint_Logo3_alt.eps"
  "${PROJECT_SOURCE_DIR}/INSTALL.md"
  "${PROJECT_SOURCE_DIR}/COPYING"
  "${PROJECT_SOURCE_DIR}/COPYING.LESSER"
  "${PROJECT_SOURCE_DIR}/CITATION"
  "${PROJECT_SOURCE_DIR}/README.md"
  "${PROJECT_BINARY_DIR}/cmake/modules/int_computed.cmake"
  )

# N.B. populate.cmake leaves unchanged files untouched, so the tarball keys off a stamp
set(POPULATE_STAMP "${CMAKE_CURRENT_BINARY_DIR}/libint-library-populate.stamp")
add_custom_command(
  OUTPUT
    "${POPULATE_STAMP}"
    ${EXPORT_STAGE_DIR}/CMakeLists.txt
    ${EXPORT_STAGE_DIR}/include/libint2_params.h
    ${EXPORT_STAGE_DIR}/include/HRRPart0bra0ket0pp.h
    ${EXPORT_STAGE_DIR}/src/HRRPart0bra0ket0pp.cc
  COMMAND ${CMAKE_COMMAND}
          "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}"
          "-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}"
          "-DEXPORT_STAGE_DIR=${EXPORT_STAGE_DIR}"
          -P ${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake
  COMMAND ${CMAKE_COMMAND} -E touch "${POPULATE_STAMP}"
  DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/populate.cmake"
    "${PROJECT_SOURCE_DIR}/cmake/libint2-config.cmake.in"
    ${_populate_inputs}
    "${PROJECT_BINARY_DIR}/generated/libint2_params.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.h"
    "${PROJECT_BINARY_DIR}/generated/HRRPart0bra0ket0pp.cc"
    libint-library-generate
  COMMENT "Populating Libint2 library"
  )
add_custom_target(libint-library-populate DEPENDS "${POPULATE_STAMP}")


# <<<  Export The Library Source  >>>

set(EXPORT_TARBALL "${EXPORT_STAGE_DIR}.${LIBINT_EXPORT_COMPRESSOR_EXT}")
add_custom_command(
  OUTPUT
    "${EXPORT_TARBALL}"
  COMMAND ${CMAKE_COMMAND} -E tar ${LIBINT_EXPORT_COMPRESSOR_CMD} "${EXPORT_TARBALL}" "${EXPORT_STAGE_DIR}"
  DEPENDS
    libint-library-generate
    libint-library-populate
    "${POPULATE_STAMP}"
  WORKING_DIRECTORY "${EXPORT_STAGE_DIR}/.."
  COMMENT "Exporting tarball of Libint2 library source: ${EXPORT_TARBALL}"
  )
if (LIBINT2_BUILD_LIBRARY_AS_SUBPROJECT)
    # rerun cmake if building as subproject to unpack and consume the library as subproject
    add_custom_target(libint-library-export DEPENDS "${EXPORT_TARBALL}"
      COMMAND ${CMAKE_COMMAND} -S "${CMAKE_SOURCE_DIR}" -B "${CMAKE_BINARY_DIR}")
else()
    add_custom_target(libint-library-export DEPENDS "${EXPORT_TARBALL}")
endif()
add_custom_target(export DEPENDS libint-library-export)


# <<<  Build Library  >>>

include(CheckLanguage)
check_language(C)
if(CMAKE_C_COMPILER)
    enable_language(C)
else()
    message(FATAL_ERROR "No C language support; invoke CMake with -D CMAKE_C_COMPILER=/path/to/C/compiler to resolve")
endif()

if (LIBINT2_BUILD_LIBRARY_AS_SUBPROJECT)
    include(FetchContent)

    if (EXISTS "${EXPORT_TARBALL}")
        message("Unpacking ${EXPORT_TARBALL} and setting up as a subproject")
        # the hash makes FetchContent re-unpack a re-exported tarball at the same URL
        file(SHA256 "${EXPORT_TARBALL}" _export_tarball_sha256)
        FetchContent_Declare(
          libint-library
          URL "${EXPORT_TARBALL}"
          URL_HASH SHA256=${_export_tarball_sha256}
          )
        FetchContent_MakeAvailable(libint-library)
    else()
        message(WARNING "LIBINT2_BUILD_LIBRARY_AS_SUBPROJECT=ON but library has not been exported yet; build target 'libint-library-export' first, this will make library build/test targets ('check', 'install') available. Don't run 'check' with CTEST_PARALLEL_LEVEL as this will build the library multiple times.")
    endif()

else()
    include(ExternalProject)

    # the library is installed into a staging area and from there into CMAKE_INSTALL_PREFIX (see
    # below); an absolute install dir would make the library install bypass the staging area, the
    # install prefix, and DESTDIR
    foreach(_dir INCLUDEDIR LIBDIR DATADIR)
      if (IS_ABSOLUTE "${CMAKE_INSTALL_${_dir}}")
        message(FATAL_ERROR "CMAKE_INSTALL_${_dir}=${CMAKE_INSTALL_${_dir}} is an absolute path; this is not supported when the library is built as an ExternalProject. Use relative install directories, set LIBINT2_BUILD_LIBRARY_AS_SUBPROJECT=ON, or build the exported library (target 'export') separately.")
      endif()
    endforeach()

    set(LIBRARY_INSTALL_STAGE_DIR "${CMAKE_CURRENT_BINARY_DIR}/library-install-stage")
    list(APPEND library_CMAKE_ARGS
        -G${CMAKE_GENERATOR}
        -DCMAKE_INSTALL_PREFIX=${LIBRARY_INSTALL_STAGE_DIR}
        -DCMAKE_INSTALL_INCLUDEDIR=${CMAKE_INSTALL_INCLUDEDIR}
        -DCMAKE_INSTALL_LIBDIR=${CMAKE_INSTALL_LIBDIR}
        -DCMAKE_INSTALL_DATADIR=${CMAKE_INSTALL_DATADIR}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DBUILD_TESTING=${BUILD_TESTING}
        -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        -DLIBINT2_BUILD_SHARED_AND_STATIC_LIBS=${LIBINT2_BUILD_SHARED_AND_STATIC_LIBS}
        -DLIBINT2_REQUIRE_CXX_API=${LIBINT2_REQUIRE_CXX_API}
        -DLIBINT2_REQUIRE_CXX_API_COMPILED=${LIBINT2_REQUIRE_CXX_API_COMPILED}
        -DLIBINT2_ENABLE_FORTRAN=${LIBINT2_ENABLE_FORTRAN}
        -DLIBINT2_ALIGN_SIZE=${LIBINT_ALIGN_SIZE}
        -DLIBINT2_REALTYPE=${LIBINT2_REALTYPE}
        -DLIBINT2_USER_DEFINED_REAL_INCLUDES=${LIBINT2_USER_DEFINED_REAL_INCLUDES}
        -DLIBINT2_ENABLE_MPFR=${LIBINT2_ENABLE_MPFR}
        -DLIBINT2_ENABLE_PYTHON=${LIBINT2_ENABLE_PYTHON}
        -DLIBINT2_PREFIX_PYTHON_INSTALL=${LIBINT2_PREFIX_PYTHON_INSTALL}
      )
    if (Multiprecision_ROOT)
        list(APPEND library_CMAKE_ARGS "-DMultiprecision_ROOT=${Multiprecision_ROOT}")
    endif()
    if (Boost_DIR)
        list(APPEND library_CMAKE_ARGS "-DBoost_DIR=${Boost_DIR}")
    endif()
    if (BOOST_ROOT)
        list(APPEND library_CMAKE_ARGS "-DBOOST_ROOT=${BOOST_ROOT}")
    endif()
    if (CMAKE_DISABLE_FIND_PACKAGE_Boost)
        list(APPEND library_CMAKE_ARGS "-DCMAKE_DISABLE_FIND_PACKAGE_Boost=${CMAKE_DISABLE_FIND_PACKAGE_Boost}")
    endif()
    if (Eigen3_DIR)
        list(APPEND library_CMAKE_ARGS "-DEigen3_DIR=${Eigen3_DIR}")
    endif()
    if (Eigen3_ROOT)
        list(APPEND library_CMAKE_ARGS "-DEigen3_ROOT=${Eigen3_ROOT}")
    endif()
    if (LIBINT2_LOCAL_Eigen3_INSTALL)
        list(APPEND library_CMAKE_ARGS "-DLIBINT2_LOCAL_Eigen3_INSTALL=${LIBINT2_LOCAL_Eigen3_INSTALL}")
    endif()
    if (Python_EXECUTABLE)
        list(APPEND library_CMAKE_ARGS "-DPython_EXECUTABLE=${Python_EXECUTABLE}")
    endif()
    if (LIBINT2_ENABLE_FORTRAN AND CMAKE_Fortran_COMPILER)
        # Note: not handling CMAKE_Fortran_COMPILER_ARG1
        list(APPEND library_CMAKE_ARGS "-DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}")
    endif()
    # N.B. the library is installed into LIBRARY_INSTALL_STAGE_DIR and then copied
    #   from there, so the module dir must be passed relative to the prefix: an
    #   absolute one (e.g. Fedora's %{_fmoddir}) would be installed into directly
    #   at build time, bypassing the stage. Accept absolute paths below
    #   CMAKE_INSTALL_PREFIX by making them relative, reject any others.
    if (LIBINT2_ENABLE_FORTRAN AND LIBINT2_INSTALL_FMODDIR)
        file(TO_CMAKE_PATH "${LIBINT2_INSTALL_FMODDIR}" _l2_fmoddir)
        if (IS_ABSOLUTE "${_l2_fmoddir}")
            file(TO_CMAKE_PATH "${CMAKE_INSTALL_PREFIX}" _l2_prefix)
            file(RELATIVE_PATH _l2_fmoddir_rel "${_l2_prefix}" "${_l2_fmoddir}")
            if (_l2_fmoddir_rel MATCHES "^\\.\\.(/|$)" OR IS_ABSOLUTE "${_l2_fmoddir_rel}")
                message(FATAL_ERROR "LIBINT2_INSTALL_FMODDIR=${LIBINT2_INSTALL_FMODDIR} is not below CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}; use a path relative to, or below, CMAKE_INSTALL_PREFIX")
            endif()
            set(_l2_fmoddir "${_l2_fmoddir_rel}")
            unset(_l2_fmoddir_rel)
            unset(_l2_prefix)
        endif()
        list(APPEND library_CMAKE_ARGS "-DLIBINT2_INSTALL_FMODDIR=${_l2_fmoddir}")
    endif()

    ExternalProject_Add(
      library
      DEPENDS
        libint-library-export
      BUILD_ALWAYS 1
      # download generalized to generating source by unpacking source tarball
      # * DOWNLOAD_DIR <SOURCE_DIR> is tidier but fails on Windows
      URL "${EXPORT_TARBALL}"
      LOG_DOWNLOAD 1
      CMAKE_ARGS "${library_CMAKE_ARGS}"
      CMAKE_CACHE_ARGS
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
        -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
        -DCMAKE_PREFIX_PATH:STRING=${CMAKE_PREFIX_PATH}
      TEST_COMMAND ""
      USES_TERMINAL_CONFIGURE TRUE
      USES_TERMINAL_BUILD TRUE
      USES_TERMINAL_INSTALL TRUE
      STEP_TARGETS
        build
        install
      )
    ExternalProject_Get_property(library BINARY_DIR)
    # the download (i.e. unpack) step otherwise only reruns if the URL changes,
    #   so a re-exported tarball would not be picked up by an existing build tree
    ExternalProject_Add_StepDependencies(library download "${EXPORT_TARBALL}")
    # start each library install from an empty stage, so that the files left
    #   behind by a differently-configured earlier build (e.g. shared libraries
    #   after switching to BUILD_SHARED_LIBS=OFF) do not get installed
    ExternalProject_Add_Step(
      library clear-stage
      COMMAND ${CMAKE_COMMAND} -E rm -rf "${LIBRARY_INSTALL_STAGE_DIR}"
      COMMENT "Clearing Libint2 library install stage"
      DEPENDEES build
      DEPENDERS install
      ALWAYS TRUE
      )

    # instead of running tests by default, replicate check-libint2 target here
    #   and use it to invoke library's check-libint2 target
    if (BUILD_TESTING)
        include(AddCustomTargetSubproject)
        add_custom_target_subproject(
          libint2 check  # target check-libint2
          USES_TERMINAL
          COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/"
          WORKING_DIRECTORY ${BINARY_DIR}
          DEPENDS
            "library-build"
          )

        if (LIBINT2_ENABLE_PYTHON)
            add_custom_target_subproject(
              libint2 check-python
              USES_TERMINAL
              COMMAND ${CMAKE_COMMAND} --build ${BINARY_DIR} --target libint2-python-test
              WORKING_DIRECTORY ${BINARY_DIR}/python
            )
        endif()
    endif()

    # install library from ExternalProject staging area into CMAKE_INSTALL_PREFIX
    # N.B. LIBINT2_INSTALL_FMODDIR defaults to CMAKE_INSTALL_INCLUDEDIR but a packager
    #   may point it outside the three dirs below (e.g. lib64/gfortran/modules while
    #   CMAKE_INSTALL_LIBDIR=lib), in which case the staged modules need their own pass
    set(_stage_dirs ${CMAKE_INSTALL_INCLUDEDIR} ${CMAKE_INSTALL_LIBDIR} ${CMAKE_INSTALL_DATADIR})
    if (DEFINED _l2_fmoddir)
        list(APPEND _stage_dirs ${_l2_fmoddir})
    endif()
    list(REMOVE_DUPLICATES _stage_dirs)
    foreach(_dir ${_stage_dirs})
      install(
        DIRECTORY "${LIBRARY_INSTALL_STAGE_DIR}/${_dir}/"
        DESTINATION "${_dir}"
        USE_SOURCE_PERMISSIONS
        )
    endforeach()

    # generate basis.h at install time so LIBINT_DATADIR_ABSOLUTE reflects any --prefix override.
    # Same approach as in export/CMakeLists.txt.export (see the comment there): configure_file into
    # the build tree, then file(INSTALL) with explicit FILE_PERMISSIONS, so that the installed
    # header does not inherit basis.h.in's mode, $ENV{DESTDIR} is honored, and basis.h is recorded
    # in install_manifest.txt.
    install(CODE "
      set(LIBINT_VERSION [==[${LIBINT_VERSION}]==])
      set(LIBINT_DATADIR_ABSOLUTE [==[${CMAKE_INSTALL_DATADIR}/libint/${LIBINT_VERSION}]==])
      set(LIBINT_DATADIR_ABSOLUTE \"\${CMAKE_INSTALL_PREFIX}/\${LIBINT_DATADIR_ABSOLUTE}\")
      set(_libint2_includedir [==[${CMAKE_INSTALL_INCLUDEDIR}]==])
      set(_libint2_includedir \"\${CMAKE_INSTALL_PREFIX}/\${_libint2_includedir}\")
      string(MD5 _libint2_prefix_hash \"\${CMAKE_INSTALL_PREFIX}\")
      set(_libint2_basis_h [==[${CMAKE_CURRENT_BINARY_DIR}/basis.h.install]==])
      set(_libint2_basis_h \"\${_libint2_basis_h}/\${_libint2_prefix_hash}/basis.h\")
      configure_file([==[${PROJECT_SOURCE_DIR}/include/libint2/basis.h.in]==] \"\${_libint2_basis_h}\" @ONLY)
      file(TOUCH \"\${_libint2_basis_h}\")
      file(INSTALL \"\${_libint2_basis_h}\"
        DESTINATION \"\${_libint2_includedir}/libint2\"
        FILE_PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
    ")
endif()
