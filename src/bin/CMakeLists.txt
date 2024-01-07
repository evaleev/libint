add_subdirectory(libint)

include(AddCustomTargetSubproject)
add_custom_target_subproject(
  libint2compiler check  # target check-libint2compiler
  USES_TERMINAL
  COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/compiler"
  )
add_dependencies(check-libint2compiler build_libint)

if (BUILD_TESTING)
    add_subdirectory(test_eri)
    # borrows main libint/test.cc from libtool layout
    add_subdirectory(profile)
endif()
