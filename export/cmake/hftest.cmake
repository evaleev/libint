
macro(runtest)

    set(OUTPUT_FILE_NAME "${CMAKE_BINARY_DIR}/${testName}.out")

    set(CHECK_CMD "${pythonExec}" "${srcDir}/tests/hartree-fock/${fileName}-validate.py")
    if (${fileName} STREQUAL "hartree-fock++")
        set(CHECK_ARGS "${srcDir}/MakeVars.features")
    endif()

    execute_process(COMMAND
            ${CMAKE_BINARY_DIR}/${testName} ${testArgs}
            OUTPUT_FILE "${OUTPUT_FILE_NAME}"
            RESULT_VARIABLE TEST_RESULT)

    if(TEST_RESULT)
        message(STATUS "\nOUTPUT of " ${testName} " with args " ${testArgs})
        execute_process(COMMAND
                cat
                ${OUTPUT_FILE_NAME}
                RESULT_VARIABLE
                CAT_RESULT
                )
        message(FATAL_ERROR "Error running ${CMAKE_BINARY_DIR}/${testName}")
    endif(TEST_RESULT)

    execute_process(COMMAND
            ${CHECK_CMD} ${CHECK_ARGS} ${OUTPUT_FILE_NAME}
            RESULT_VARIABLE CHECK_RESULT)

    if(CHECK_RESULT)
        message(STATUS "\nOUTPUT of " ${testName})
        execute_process(COMMAND
                cat
                ${OUTPUT_FILE_NAME}
                RESULT_VARIABLE
                CAT_RESULT
                )
        message(FATAL_ERROR "Error running ${CHECK_CMD} with args " ${CHECK_ARGS} ${OUTPUT_FILE_NAME})
    endif(CHECK_RESULT)

endmacro(runtest)

runtest()
