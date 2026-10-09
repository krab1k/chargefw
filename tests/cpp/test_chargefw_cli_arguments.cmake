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

function(expect_success label)
    execute_process(
            COMMAND "${CHARGEFW_CLI}" ${ARGN}
            RESULT_VARIABLE result
            ERROR_VARIABLE error
    )
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "${label} exit status was ${result}: ${error}")
    endif()
endfunction()

expect_argument_error(missing_subcommand)
expect_argument_error(missing_required_argument calculate "${CHARGEFW_INPUT}")
expect_argument_error(invalid_threads calculate --threads not-a-count "${CHARGEFW_INPUT}" output)
expect_argument_error(unknown_option calculate --unknown-option "${CHARGEFW_INPUT}" output)
expect_argument_error(unknown_method methods not-a-method)
expect_argument_error(unknown_parameter_set parameters not-a-parameter-set)
expect_argument_error(unknown_parameter_method parameters --method not-a-method)
expect_argument_error(parameter_filter_and_detail parameters --method qeq QEq_original)
expect_success(method_option_preserves_positionals applicability --method-option peoe.iters=8
               "${CHARGEFW_INPUT}")

expect_success(help --help)
expect_success(calculate_help calculate --help)
expect_success(inspect inspect "${CHARGEFW_INPUT}")
expect_success(applicable applicability --method formal "${CHARGEFW_INPUT}")
expect_success(inapplicable applicability --method smpqeq "${CHARGEFW_INPUT}")
expect_success(methods methods)
expect_success(method_details methods eem)
expect_success(parameters parameters)
expect_success(parameter_filter parameters --method qeq)
expect_success(parameter_details parameters QEq_original)
