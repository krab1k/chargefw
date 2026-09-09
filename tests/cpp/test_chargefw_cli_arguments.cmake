include("${CMAKE_CURRENT_LIST_DIR}/prepare_moved_install.cmake")

function(expect_argument_error label)
    execute_process(
            COMMAND "${CHARGEFW_CLI}" ${ARGN}
            RESULT_VARIABLE result
            ERROR_VARIABLE error
    )
    if(NOT result EQUAL 2)
        message(FATAL_ERROR "${label} exit status was ${result}: ${error}")
    endif()
endfunction()

expect_argument_error(missing_required_argument calculate "${CHARGEFW_INPUT}")
expect_argument_error(invalid_threads calculate --threads not-a-count "${CHARGEFW_INPUT}" output)
expect_argument_error(unknown_option calculate --unknown-option "${CHARGEFW_INPUT}" output)
expect_argument_error(unknown_method methods not-a-method)

execute_process(
        COMMAND "${CHARGEFW_CLI}" --help
        RESULT_VARIABLE help_result
        ERROR_VARIABLE help_error
)
if(NOT help_result EQUAL 0)
    message(FATAL_ERROR "help command failed with exit status ${help_result}: ${help_error}")
endif()
