cmake_minimum_required(VERSION 3.12)

project(libint2-python)

set(CMAKE_CXX_STANDARD 17)
add_compile_options(-Wall)

include(FetchContent)

FetchContent_Declare(
  pybind11
  GIT_REPOSITORY      https://github.com/ValeevGroup/pybind11.git
  GIT_TAG             80d452484c5409444b0ec19383faa84bb7a4d351 # v2.4.3
  )
FetchContent_MakeAvailable(pybind11)

#find_package(Eigen3 3.3 REQUIRED)

pybind11_add_module(
  libint2-python MODULE
  #EXCLUDE_FROM_ALL
  src/libint2/libint2.cc
  src/libint2/engine.cc
  )

set_target_properties(
  libint2-python
  PROPERTIES
  #PREFIX ""
  OUTPUT_NAME libint2
  )

if (TARGET libint2_obj)
  set(libint2_python_target libint2_obj)
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
  COMMAND ${PYTHON_EXECUTABLE} -m setup test
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

add_custom_target(
  libint2-python-wheel
  DEPENDS libint2-python
  COMMAND ${PYTHON_EXECUTABLE} -m setup bdist_wheel
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )

enable_testing()

# add the executable
add_test(
  NAME libint2-python
  COMMAND ${PYTHON_EXECUTABLE} -m setup test
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
  )
