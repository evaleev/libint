if (CMAKE_HOST_WIN32)
    cmake_minimum_required(VERSION 3.21)
else()
    cmake_minimum_required(VERSION 3.16)
endif()
cmake_policy(SET CMP0079 NEW)

project(libint2-python)

if (NOT TARGET Python::Module)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

find_package(pybind11 2.6.0 CONFIG)
if (NOT TARGET pybind11::pybind11)
    message(STATUS "Suitable pybind11 could not be located, building pybind11 instead.")

    # [LAB May 2023] For a long time this was fixed at v2.4.3 from a Valeev group clone.
    #  * recently noticed that C++17 is not getting processed right for this version (at least with the Ubuntu GHA setup)
    #  * so bumping the backup tag to v2.6.0 (c.2021). the v2.4.3 could probably be revived if needed.
    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v2.6.0  # sync below!
      # GIT_REPOSITORY      https://github.com/ValeevGroup/pybind11.git
      # GIT_TAG             80d452484c5409444b0ec19383faa84bb7a4d351  # v2.4.3
      # FIND_PACKAGE_ARGS 2.4.3 CONFIG  # CMake 3.24 integrates find_package() call above into FC_Declare
      )

    FetchContent_MakeAvailable(pybind11)
    set(pybind11_VERSION "2.6.0")  # getting the fetched version isn't reliable from ${CMAKE_PROJECT_VERSION}), so sync explicitly from tag above
endif()

if (pybind11_VERSION VERSION_LESS_EQUAL 2.5.0)
    # remove after minimum pb11 advances beyond v2.5.0 and let cxx_std_17 below do the work
    set(CMAKE_CXX_STANDARD 17)
endif()

if (pybind11_VERSION VERSION_GREATER_EQUAL 2.11.0)
    # ipo/lto triggering changed https://github.com/pybind/pybind11/pull/4643
    if ((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") AND (CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC"))
        # clang-cl

        # this variable needs to be defined, not the property
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
    endif()
endif()

pybind11_add_module(
  libint2-python MODULE
  #EXCLUDE_FROM_ALL
  src/libint2/libint2.cc
  src/libint2/engine.cc
  )

target_compile_features(libint2-python PRIVATE "cxx_std_17")

if (pybind11_VERSION VERSION_GREATER_EQUAL 2.7.0)
    # suppress error in L2
    # * https://github.com/pybind/pybind11/pull/2919
    target_compile_definitions(libint2-python PRIVATE NDEBUG=1)
endif()

target_compile_options(libint2-python
  PRIVATE
    # too many warnings on Windows
    # $<$<CXX_COMPILER_ID:MSVC>:/W4>
    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
  )

#find_package(Eigen3 3.3 REQUIRED)

set_target_properties(
  libint2-python
  PROPERTIES
  #PREFIX ""
  OUTPUT_NAME libint2
  )

if (TARGET libint2_obj)
  set(libint2_python_target libint2_obj)

    if(MSVC)
        target_compile_definitions(libint2-python PUBLIC _USE_MATH_DEFINES)
        target_compile_options(libint2-python PUBLIC "/EHsc")
    endif()

  target_link_libraries(libint2-python PRIVATE libint2_obj)
  target_link_libraries(libint2-python PRIVATE Boost::boost)
else()
  find_package(Libint2 REQUIRED)
  target_link_libraries(libint2-python PRIVATE Libint2::cxx)
endif()

# if (Eigen3::Eigen)
#   target_link_libraries(libint2-python INTERFACE Eigen3::Eigen)
# else()
#   include_directories(${EIGEN3_INCLUDE_DIR})
# endif()

configure_file(setup.py.in ${PROJECT_BINARY_DIR}/setup.py)

add_custom_target(
  libint2-python-test
  DEPENDS libint2-python
  COMMAND ${Python_EXECUTABLE} -m setup test
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

add_custom_target(
  libint2-python-wheel
  DEPENDS libint2-python
  COMMAND ${Python_EXECUTABLE} -m setup bdist_wheel
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

enable_testing()

# add the executable
add_test(
  NAME libint2-python
  COMMAND ${Python_EXECUTABLE} -m setup test
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
