set(input_path "${CHARGEFW_INPUT}")
set(output_directory "${CHARGEFW_TEST_DIR}/chargefw_cli_fixed_ions")
set(output_json "${output_directory}/fixed_ions.chargefw.json")
execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen
                --fixed-ions MG --fixed-ions CA "${input_path}"
        RESULT_VARIABLE applicability_result
        ERROR_VARIABLE applicability_error
)
if(NOT applicability_result EQUAL 0)
    message(FATAL_ERROR "fixed-charge applicability failed: ${applicability_error}")
endif()

file(REMOVE_RECURSE "${output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method sqeqp
                --parameter-set SQEqp_Schindler2021_CCD_gen --execution full
                --fixed-ions MG --fixed-ions MG --fixed-ions CA
                "${input_path}" "${output_directory}"
        RESULT_VARIABLE calculate_result
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
string(JSON fixed_groups GET "${result_json}" calculation_provenance effective fixed_ions)
string(JSON group_count LENGTH "${fixed_groups}")
if(NOT group_count EQUAL 2)
    message(FATAL_ERROR "missing effective fixed-charge provenance: ${fixed_groups}")
endif()
string(JSON mg_id GET "${fixed_groups}" 0 component_id)
string(JSON mg_index GET "${fixed_groups}" 0 instances 0 atom_index)
string(JSON mg_value GET "${fixed_groups}" 0 charge)
string(JSON ca_id GET "${fixed_groups}" 1 component_id)
string(JSON ca_index GET "${fixed_groups}" 1 instances 0 atom_index)
string(JSON ca_value GET "${fixed_groups}" 1 charge)
if(NOT mg_id STREQUAL "MG" OR NOT ca_id STREQUAL "CA" OR
   NOT mg_index EQUAL 3 OR NOT mg_value EQUAL 2 OR NOT ca_index EQUAL 4 OR NOT ca_value EQUAL 2)
    message(FATAL_ERROR "fixed sources differ from input order or preset values: ${fixed_groups}")
endif()
foreach(removed_field charge_provenance charge_totals original_total_charge active_total_charge)
    string(JSON removed_type ERROR_VARIABLE missing_field TYPE "${fixed_groups}" 0 "${removed_field}")
    if(NOT missing_field)
        message(FATAL_ERROR "unexpected fixed-ion audit field: ${removed_field}")
    endif()
endforeach()

execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --fixed-ions UNKNOWN "${input_path}"
        RESULT_VARIABLE unknown_result
        ERROR_VARIABLE unknown_error
)
if(NOT unknown_result EQUAL 2)
    message(FATAL_ERROR "unknown component ID was not rejected: ${unknown_error}")
endif()
execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability "${input_path}" --fixed-ions
        RESULT_VARIABLE missing_value_result
        ERROR_VARIABLE missing_value_error
)
if(NOT missing_value_result EQUAL 2)
    message(FATAL_ERROR "missing component ID was not reported by argument parsing: ${missing_value_error}")
endif()

file(REMOVE_RECURSE "${output_directory}")
