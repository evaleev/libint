if (CMAKE_HOST_WIN32)
    cmake_minimum_required(VERSION 3.21)
else()
    cmake_minimum_required(VERSION 3.16)
endif()
cmake_policy(SET CMP0079 NEW)

project(libint2-python)

#  <<<  Build  >>>

if (NOT TARGET Python::Module)
    find_package(Python COMPONENTS Interpreter Development REQUIRED)
endif()

find_package(pybind11 2.4.3 CONFIG)
if (NOT TARGET pybind11::pybind11)
    message(STATUS "Suitable pybind11 could not be located, building pybind11 instead.")

    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v2.4.3  # sync below!
      # GIT_REPOSITORY      https://github.com/ValeevGroup/pybind11.git
      # GIT_TAG             80d452484c5409444b0ec19383faa84bb7a4d351  # v2.4.3
      # FIND_PACKAGE_ARGS 2.4.3 CONFIG  # CMake 3.24 integrates find_package() call above into FC_Declare
      )

    FetchContent_MakeAvailable(pybind11)
    set(pybind11_VERSION "2.4.3")  # getting the fetched version isn't reliable from ${CMAKE_PROJECT_VERSION}), so sync explicitly from tag above
endif()

if (pybind11_VERSION VERSION_LESS_EQUAL 2.5.0)
    # remove after minimum pb11 advances beyond v2.5.0 and let cxx_std_17 below do the work
    set(CMAKE_CXX_STANDARD 17)
endif()

pybind11_add_module(
  libint2-python
  MODULE
  # EXCLUDE_FROM_ALL  # why?
    src/libint2/libint2.cc
    src/libint2/engine.cc
  )

target_compile_features(libint2-python PRIVATE "cxx_std_17")

if (pybind11_VERSION VERSION_GREATER_EQUAL 2.7.0)
    # suppress error in L2
    # * https://github.com/pybind/pybind11/pull/2919
    target_compile_definitions(libint2-python PRIVATE NDEBUG=1)
endif()

# * too many warnings on Windows
# target_compile_options(libint2-python
#   PRIVATE
#     $<$<CXX_COMPILER_ID:MSVC>:/W4>
#     $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall>
#   )

set_target_properties(
  libint2-python
  PROPERTIES
    OUTPUT_NAME libint2
  )

# future note: this whole block could likely be replaced by the below, but perhaps explicit better for now
# * `if(TARGET Libint2::cxx)\n  target_link_libraries(libint2-python PRIVATE Libint2::cxx)\nelse()\n  find_package(Libint2 ...`
if (TARGET int-obj)
  #set(libint2_python_target int-obj)

    target_compile_definitions(libint2-python PRIVATE __COMPILING_LIBINT2)

    if(MSVC)
        target_compile_definitions(libint2-python PUBLIC _USE_MATH_DEFINES)
        target_compile_options(libint2-python PUBLIC "/EHsc")
        set_target_properties(
          libint2-python
          PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS ON
          )
    endif()

    target_include_directories(
      libint2-python
      PRIVATE
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/../src>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/../include>
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/../include/libint2>
      )

    target_link_libraries(
      libint2-python
      PRIVATE
        int-obj
        Eigen3::Eigen
      )
    if (TARGET Boost::headers)
        target_link_libraries(libint2-python PRIVATE Boost::headers)
    endif()

else()
    find_package(Libint2 COMPONENTS CXX_ho REQUIRED)
    target_link_libraries(libint2-python PRIVATE Libint2::cxx)
endif()

configure_file(setup.py.in ${PROJECT_BINARY_DIR}/setup.py)

#  <<<  Install  >>>

if (LIBINT2_PREFIX_PYTHON_INSTALL)
    # * Note that this block is *Linux-style* install to `CMAKE_INSTALL_PREFIX`.
    # * For *Python-style* install to `Python_EXECUTABLE`'s site-packages, use
    #   wheel target below from python/ directory.
    #     > make libint2-python-wheel
    #     > pip3 install dist/libint2-*.whl

    execute_process(
      COMMAND ${Python_EXECUTABLE} -c
        "import os, sys, sysconfig as s; spdir = s.get_path('platlib'); print(spdir.replace(os.path.commonpath([sys.prefix, spdir]), '').lstrip(os.path.sep));"
      OUTPUT_VARIABLE LIBINT2_INSTALL_PYMODDIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    message(STATUS "Showing option LIBINT2_INSTALL_PYMODDIR: ${LIBINT2_INSTALL_PYMODDIR} (non-user-configurable until wanted)")

    install(
      FILES $<TARGET_FILE:libint2-python>
      COMPONENT ${L2}_Python
      DESTINATION ${LIBINT2_INSTALL_PYMODDIR}
      )
endif()

#  <<<  Python Tooling  >>>

add_custom_target(
  libint2-python-check
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
