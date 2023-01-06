include(AddCustomTargetSubproject)
add_custom_target_subproject(libint2compiler check USES_TERMINAL COMMAND ${CMAKE_CTEST_COMMAND} -V -R "libint2/compiler")

add_subdirectory(core)
add_subdirectory(eri)

