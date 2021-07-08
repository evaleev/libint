add_executable(libint-generate-eri-test EXCLUDE_FROM_ALL generate_eri_code.cc)
target_link_libraries(libint-generate-eri-test libint-libcompiler)
target_include_directories(libint-generate-eri-test PRIVATE ${PROJECT_SOURCE_DIR}/src/bin/)

set_target_properties(libint-generate-eri-test
    PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)

