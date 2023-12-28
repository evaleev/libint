add_executable(
  libint-compiler-eri-test-generate
  EXCLUDE_FROM_ALL
  generate_eri_code.cc
  )
target_link_libraries(libint-compiler-eri-test-generate libint-libcompiler)
target_include_directories(
  libint-compiler-eri-test-generate
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src/bin/
  )
set_target_properties(
  libint-compiler-eri-test-generate
  PROPERTIES
    INSTALL_RPATH_USE_LINK_PATH TRUE
  )
