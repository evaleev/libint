add_subdirectory(libint)

include(AddCustomTargetSubproject)
add_custom_target_subproject(
  libint2compiler check  # target check-libint2compiler
  USES_TERMINAL
  COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/compiler"
  )
add_dependencies(check-libint2compiler build_libint)

if (BUILD_TESTING)
  # the generated code must not depend on the order in which the generator
  # emits the shell sets and the set-level recurrence relations, hence on the
  # number of processes that generate it (LIBINT2_GENERATOR_JOBS); the library
  # generation target is defined in src/lib/libint
  add_test(
    NAME libint2/compiler/generate
    COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --config $<CONFIG> --target libint-library-generate
    )
  set_tests_properties(
    libint2/compiler/generate
    PROPERTIES
      FIXTURES_SETUP LIBINT2_COMPILER_GENERATED
    )
  add_test(
    NAME libint2/compiler/codegen_order
    COMMAND ${CMAKE_COMMAND}
      -DBUILD_LIBINT=$<TARGET_FILE:build_libint>
      -DREFERENCE_DIR=${PROJECT_BINARY_DIR}/generated
      -DWORK_DIR=${PROJECT_BINARY_DIR}/generated-codegen-order
      -P ${CMAKE_CURRENT_SOURCE_DIR}/check_codegen_order.cmake
    )
  set_tests_properties(
    libint2/compiler/codegen_order
    PROPERTIES
      FIXTURES_REQUIRED LIBINT2_COMPILER_GENERATED
    )

  add_subdirectory(test_eri)
  # borrows main libint/test.cc from libtool layout
  add_subdirectory(profile)
endif()
