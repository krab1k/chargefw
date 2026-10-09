function(run_cli label expected_exit)
    execute_process(
            COMMAND "${CHARGEFW_CLI}" ${ARGN}
            RESULT_VARIABLE result
            OUTPUT_VARIABLE output
            ERROR_VARIABLE error
    )
    if(NOT "${result}" STREQUAL "${expected_exit}")
        message(FATAL_ERROR "${label}: expected exit ${expected_exit}, got ${result}\n${output}${error}")
    endif()
endfunction()
