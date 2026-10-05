set(input_path "${CHARGEFW_INPUT}")
set(output_directory "${CHARGEFW_TEST_DIR}/chargefw_cli_fixed_charge_groups")
set(output_json "${output_directory}/fixed_charge_groups.chargefw.json")
execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen
                --fixed-charge-group MG --fixed-charge-group CA "${input_path}"
        RESULT_VARIABLE applicability_result
        OUTPUT_VARIABLE applicability_output
        ERROR_VARIABLE applicability_error
)
if(NOT applicability_result EQUAL 0 OR
   NOT applicability_output MATCHES "runnable plans: [1-9]" OR
   NOT applicability_output MATCHES "plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen")
    message(FATAL_ERROR "fixed-charge applicability did not produce SQEqp: ${applicability_output}${applicability_error}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen "${input_path}"
        RESULT_VARIABLE unselected_result
        OUTPUT_VARIABLE unselected_output
        ERROR_VARIABLE unselected_error
)
if(NOT unselected_result EQUAL 0 OR
   NOT unselected_output MATCHES "runnable plans: 0" OR
   NOT unselected_output MATCHES "has no atom parameter matching")
    message(FATAL_ERROR "unselected ions should remain active and lack parameters: ${unselected_output}${unselected_error}")
endif()

file(REMOVE_RECURSE "${output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen --execution full
                --fixed-charge-group MG --fixed-charge-group CA "${input_path}" "${output_directory}"
        RESULT_VARIABLE calculate_result
        OUTPUT_VARIABLE calculate_output
        ERROR_VARIABLE calculate_error
)
if(NOT calculate_result EQUAL 0)
    message(FATAL_ERROR "fixed-charge calculation failed: ${calculate_error}")
endif()
file(READ "${output_json}" result_json)
string(JSON status GET "${result_json}" results 0 status)
string(JSON charges GET "${result_json}" results 0 assignments 0 charges)
string(JSON charge_count LENGTH "${charges}")
if(NOT status STREQUAL "success" OR NOT charge_count EQUAL 5)
    message(FATAL_ERROR "unexpected fixed-charge result status or charge count")
endif()
string(JSON mg_charge GET "${charges}" 3)
string(JSON ca_charge GET "${charges}" 4)
if(NOT mg_charge EQUAL 2 OR NOT ca_charge EQUAL 2)
    message(FATAL_ERROR "fixed source charges were not restored in source atom order: ${charges}")
endif()
string(JSON total_charge GET "${result_json}" results 0 assignments 0 total_charge)
if(NOT total_charge EQUAL 4)
    message(FATAL_ERROR "expected active water total plus two fixed ions to equal 4, got ${total_charge}")
endif()
string(JSON fixed_groups GET "${result_json}" calculation_provenance effective fixed_charge_groups)
string(JSON provenance GET "${fixed_groups}" charge_provenance)
string(JSON source_count LENGTH "${fixed_groups}" sources)
if(NOT provenance STREQUAL "chargefw:fixed-charge-ions:v1" OR NOT source_count EQUAL 2)
    message(FATAL_ERROR "missing effective fixed-charge provenance: ${fixed_groups}")
endif()
string(JSON mg_index GET "${fixed_groups}" sources 0 atom_index)
string(JSON mg_value GET "${fixed_groups}" sources 0 charge)
string(JSON ca_index GET "${fixed_groups}" sources 1 atom_index)
string(JSON ca_value GET "${fixed_groups}" sources 1 charge)
if(NOT mg_index EQUAL 3 OR NOT mg_value EQUAL 2 OR NOT ca_index EQUAL 4 OR NOT ca_value EQUAL 2)
    message(FATAL_ERROR "fixed sources differ from input order or preset values: ${fixed_groups}")
endif()
string(JSON active_total GET "${fixed_groups}" charge_totals 0 active_total_charge)
if(NOT active_total EQUAL 0)
    message(FATAL_ERROR "expected zero active formal-charge total for water, got ${active_total}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen --execution full
                --fixed-charge-group MG --fixed-charge-group MG --fixed-charge-group CA
                "${input_path}" "${output_directory}"
        RESULT_VARIABLE repeated_result
        ERROR_VARIABLE repeated_error
)
if(NOT repeated_result EQUAL 0)
    message(FATAL_ERROR "repeated component IDs should be idempotent: ${repeated_error}")
endif()
file(READ "${output_json}" repeated_json)
string(JSON repeated_sources LENGTH "${repeated_json}"
       calculation_provenance effective fixed_charge_groups sources)
if(NOT repeated_sources EQUAL 2)
    message(FATAL_ERROR "repeated component IDs duplicated resolved sources")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen --fixed-charge-group FE "${input_path}"
        RESULT_VARIABLE absent_result
        OUTPUT_VARIABLE absent_output
        ERROR_VARIABLE absent_error
)
if(NOT absent_result EQUAL 0 OR NOT absent_output STREQUAL "${unselected_output}")
    message(FATAL_ERROR "an absent component should match unselected assessment behavior: ${absent_output}${absent_error}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --fixed-charge-group UNKNOWN "${input_path}"
        RESULT_VARIABLE unknown_result
        ERROR_VARIABLE unknown_error
)
if(NOT unknown_result EQUAL 2 OR NOT unknown_error MATCHES "unknown fixed-charge component ID: UNKNOWN")
    message(FATAL_ERROR "unknown component ID was not reported clearly: ${unknown_error}")
endif()
execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability "${input_path}" --fixed-charge-group
        RESULT_VARIABLE missing_value_result
        ERROR_VARIABLE missing_value_error
)
if(NOT missing_value_result EQUAL 2 OR NOT missing_value_error MATCHES "--fixed-charge-group")
    message(FATAL_ERROR "missing component ID was not reported by argument parsing: ${missing_value_error}")
endif()

file(REMOVE_RECURSE "${output_directory}")
