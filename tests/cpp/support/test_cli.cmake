include_guard(GLOBAL)

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

# Sets <out_var> to an empty, nonexistent output directory named chargefw_cli_<name>.
function(fresh_output_directory out_var name)
    set(directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${name}")
    file(REMOVE_RECURSE "${directory}")
    set(${out_var} "${directory}" PARENT_SCOPE)
endfunction()

# Fails unless the JSON value at the member path given after <expected> equals <expected>.
function(expect_json label json expected)
    string(JSON actual GET "${json}" ${ARGN})
    if(NOT actual STREQUAL "${expected}")
        message(FATAL_ERROR "${label}: expected '${expected}' at '${ARGN}', got '${actual}'")
    endif()
endfunction()

function(expect_files_exist label)
    foreach(path IN LISTS ARGN)
        if(NOT EXISTS "${path}")
            message(FATAL_ERROR "${label}: expected file was not created: ${path}")
        endif()
    endforeach()
endfunction()

function(expect_files_absent label)
    foreach(path IN LISTS ARGN)
        if(EXISTS "${path}")
            message(FATAL_ERROR "${label}: unexpected file was created: ${path}")
        endif()
    endforeach()
endfunction()

# Runs the command given after <label> and fails with its error output unless it succeeds.
function(run_checked label)
    execute_process(COMMAND ${ARGN} RESULT_VARIABLE result ERROR_VARIABLE error)
    if(NOT result EQUAL 0)
        message(FATAL_ERROR "${label} failed: ${error}")
    endif()
endfunction()
