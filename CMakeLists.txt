cmake_minimum_required(VERSION 3.19)  # string(json
cmake_policy(SET CMP0074 NEW)
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
    # Avoid warning about DOWNLOAD_EXTRACT_TIMESTAMP in CMake 3.24:
    cmake_policy(SET CMP0135 NEW)
endif()

############################# Version and Metadata #############################

# can't use PROJECT_SOURCE_DIR etc. before project() call
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
include(DynamicVersion)
dynamic_version(
  PROJECT_PREFIX Libint2Compiler_
  GIT_ARCHIVAL_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/.git_archival.txt
  OUTPUT_COMMIT LibintRepository_COMMIT
  OUTPUT_VERSION LibintRepository_VERSION 
  OUTPUT_DESCRIBE LibintRepository_DESCRIBE
  OUTPUT_DISTANCE LibintRepository_DISTANCE
  )

project(
  Libint2Compiler
  VERSION ${LibintRepository_VERSION}
  DESCRIPTION
    "Libint: A library for the evaluation of molecular integrals of many-body operators over Gaussian functions"
  HOMEPAGE_URL "http://libint.valeyev.net"
  LANGUAGES CXX
  )
  # * http://libint.valeyev.net/ redirects to https://github.com/evaleev/libint

set(${PROJECT_NAME}_AUTHORS "Edward F. Valeev")
set(${PROJECT_NAME}_LICENSE "GPL-3.0 for generator; LGPL-3.0 for generated")

# along with project(... VERSION) above scanned by dynamic_version() from `git
#   describe`, these are the authoritative version source (formerly in configure.ac)
set(LIBINT_BUILDID "post999")
set(LIBINT_SOVERSION "2:3:0")

include(int_versions)

set(L2 Libint2)  # Namespace
set(pnv libint2) # projectnameversion

################################### Overview ###################################

# CMake build overview:
#
#    >>> ls
#    cmake/  COPYING  src/  tests/  ...
#    >>> cmake -S. -Bbuild -GNinja -DCMAKE_INSTALL_PREFIX=/path/to/install-libint ...
#    ...
#    -- Generating done
#    -- Build files have been written to: /current/dir/build
#    >>> cmake --build build --target install

################################## Main Project #################################

# STRICTLY TEMPORARY FOR DEMONSTRATION PURPOSES
configure_file(src/lib/libint/configuration.cc configuration.cc @ONLY)
