
macro(runtest)

    set(OUTPUT_FILE_NAME "${PROJECT_BINARY_DIR}/tests/${testName}.out")

    set(CHECK_CMD "${pythonExec}" "${PROJECT_SOURCE_DIR}/tests/hartree-fock/${testName}-validate.py")
    if (testName STREQUAL "hartree-fock++")
        set(CHECK_ARGS "${PROJECT_SOURCE_DIR}/features")
    endif()

    execute_process(COMMAND
            ${PROJECT_BINARY_DIR}/tests/${execName} ${testArgs}
            OUTPUT_FILE "${OUTPUT_FILE_NAME}"
            RESULT_VARIABLE TEST_RESULT)

    if(TEST_RESULT)
        message(STATUS "\nOUTPUT of " ${PROJECT_BINARY_DIR}/tests/${execName} " with args " ${testArgs})
        execute_process(COMMAND
                cat
                ${OUTPUT_FILE_NAME}
                RESULT_VARIABLE
                CAT_RESULT
                )
        message(FATAL_ERROR "Error running ${PROJECT_BINARY_DIR}/tests/${execName}")
    endif(TEST_RESULT)

    execute_process(COMMAND
            ${CHECK_CMD} ${CHECK_ARGS} ${OUTPUT_FILE_NAME}
            RESULT_VARIABLE CHECK_RESULT)

    if(CHECK_RESULT)
        message(STATUS "\nOUTPUT of " ${CHECK_CMD})
        execute_process(COMMAND
                cat
                ${OUTPUT_FILE_NAME}
                RESULT_VARIABLE
                CAT_RESULT
                )
        message(FATAL_ERROR "Error running ${CHECK_CMD} with args " ${CHECK_ARGS} " " ${OUTPUT_FILE_NAME})
    endif(CHECK_RESULT)

endmacro(runtest)

runtest()
