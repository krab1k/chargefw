set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_outputs")
set(output_prefix "${output_directory}/water.chargefw")
file(REMOVE_RECURSE "${output_directory}")

set(json_failure_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_json_failure")
set(json_failure_output_prefix "${json_failure_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${json_failure_output_directory}")
file(MAKE_DIRECTORY "${json_failure_output_prefix}.json")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method eem "${CHARGEFW_INPUT}"
                "${json_failure_output_directory}"
        RESULT_VARIABLE json_failure_result
        ERROR_VARIABLE json_failure_error
)
if(NOT json_failure_result EQUAL 2)
    message(FATAL_ERROR "JSON output failure was not reported correctly: ${json_failure_error}")
endif()
file(REMOVE_RECURSE "${json_failure_output_directory}")

set(export_failure_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_export_failure")
set(export_failure_output_prefix "${export_failure_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${export_failure_output_directory}")
file(MAKE_DIRECTORY "${export_failure_output_prefix}.mol2")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method eem --output-mol2 "${CHARGEFW_INPUT}"
                "${export_failure_output_directory}"
        RESULT_VARIABLE export_failure_result
        ERROR_VARIABLE export_failure_error
)
if(NOT export_failure_result EQUAL 6 OR
   NOT EXISTS "${export_failure_output_prefix}.json")
    message(FATAL_ERROR "molecular export failure was not reported correctly: ${export_failure_error}")
endif()
file(READ "${export_failure_output_prefix}.json" export_failure_json)
string(JSON export_failure_status GET "${export_failure_json}" status)
if(NOT export_failure_status STREQUAL "success")
    message(FATAL_ERROR "molecular export failure changed the JSON calculation status")
endif()
file(REMOVE_RECURSE "${export_failure_output_directory}")

foreach(mode IN ITEMS full cutoff cover)
    set(mode_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${mode}")
    set(mode_output_prefix "${mode_output_directory}/water.chargefw")
    file(REMOVE_RECURSE "${mode_output_directory}")

    set(execution_arguments --method eem --execution ${mode})
    if(NOT mode STREQUAL "full")
        list(APPEND execution_arguments --radius 8)
    endif()
    if(mode STREQUAL "cutoff")
        list(APPEND execution_arguments --progress)
    endif()

    execute_process(
            COMMAND "${CHARGEFW_CLI}" calculate ${execution_arguments} "${CHARGEFW_INPUT}"
                    "${mode_output_directory}"
            RESULT_VARIABLE mode_result
            ERROR_VARIABLE mode_error
    )
    if(NOT mode_result EQUAL 0)
        message(FATAL_ERROR "${mode} CLI calculation failed: ${mode_error}")
    endif()
    file(READ "${mode_output_prefix}.json" mode_json)
    string(JSON mode_status GET "${mode_json}" results 0 status)
    string(JSON effective_mode GET "${mode_json}" calculation_provenance effective execution mode)
    if(NOT mode_status STREQUAL "success" OR NOT effective_mode STREQUAL "${mode}")
        message(FATAL_ERROR "Unexpected ${mode} CLI result")
    endif()
    file(REMOVE_RECURSE "${mode_output_directory}")
endforeach()

function(expect_invalid_policy label)
    set(policy_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_invalid_${label}")
    file(REMOVE_RECURSE "${policy_output_directory}")
    execute_process(
            COMMAND "${CHARGEFW_CLI}" calculate ${ARGN} "${CHARGEFW_INPUT}" "${policy_output_directory}"
            RESULT_VARIABLE policy_result
            ERROR_VARIABLE policy_error
    )
    if(NOT policy_result EQUAL 2)
        message(FATAL_ERROR "${label} policy exit status was ${policy_result}: ${policy_error}")
    endif()
    if(EXISTS "${policy_output_directory}/water.chargefw.json")
        message(FATAL_ERROR "${label} policy unexpectedly wrote a result document")
    endif()
    file(REMOVE_RECURSE "${policy_output_directory}")
endfunction()

expect_invalid_policy(missing_radius --execution cutoff)
expect_invalid_policy(short_radius --execution cover --radius 7)
expect_invalid_policy(parameter_without_method --parameter-set QEq_original)

set(malformed_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_malformed")
file(REMOVE_RECURSE "${malformed_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate "${CHARGEFW_MALFORMED_INPUT}" "${malformed_output_directory}"
        RESULT_VARIABLE malformed_result
        ERROR_VARIABLE malformed_error
)
if(NOT malformed_result EQUAL 2 OR
   EXISTS "${malformed_output_directory}/malformed_then_water.chargefw.json")
    message(FATAL_ERROR "malformed input did not fail before writing output: ${malformed_error}")
endif()
file(REMOVE_RECURSE "${malformed_output_directory}")

set(mixed_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_mixed")
set(mixed_output_prefix "${mixed_output_directory}/mixed_v2000_v3000.chargefw")
file(REMOVE_RECURSE "${mixed_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method eem --execution full "${CHARGEFW_MIXED_INPUT}"
                "${mixed_output_directory}"
        RESULT_VARIABLE mixed_result
        ERROR_VARIABLE mixed_error
)
if(NOT mixed_result EQUAL 0)
    message(FATAL_ERROR "mixed V2000/V3000 CLI calculation failed: ${mixed_error}")
endif()
file(READ "${mixed_output_prefix}.json" mixed_json)
string(JSON mixed_record_count LENGTH "${mixed_json}" results)
string(JSON mixed_first_id GET "${mixed_json}" results 0 input record_id)
string(JSON mixed_second_id GET "${mixed_json}" results 1 input record_id)
if(NOT mixed_record_count EQUAL 2 OR NOT mixed_first_id STREQUAL "v2000" OR
   NOT mixed_second_id STREQUAL "v3000")
    message(FATAL_ERROR "mixed V2000/V3000 records were not preserved")
endif()
file(REMOVE_RECURSE "${mixed_output_directory}")

set(multiconformer_input "${CMAKE_CURRENT_BINARY_DIR}/multi_conformer.json")
set(multiconformer_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_multiconformer")
set(multiconformer_output_prefix "${multiconformer_output_directory}/multi_conformer.chargefw")
file(WRITE "${multiconformer_input}" [=[
{
  "schema_version": "1.0",
  "molecules": [{
    "id": "oxygen",
    "atoms": [{"atomic_number": 8, "formal_charge": 0}],
    "conformers": [
      {"id": "first", "coordinates": [[1.25, 2.5, 3.75]]},
      {"id": "second", "coordinates": [[9.5, 8.5, 7.5]]}
    ]
  }]
}
]=])
file(REMOVE_RECURSE "${multiconformer_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method formal --output-mol2 --output-mmcif
                "${multiconformer_input}"
                "${multiconformer_output_directory}"
        RESULT_VARIABLE multiconformer_result
        ERROR_VARIABLE multiconformer_error
)
if(NOT multiconformer_result EQUAL 0)
    message(FATAL_ERROR "multi-conformer JSON CLI calculation failed: ${multiconformer_error}")
endif()
foreach(extension IN ITEMS json mol2 cif)
    if(NOT EXISTS "${multiconformer_output_prefix}.${extension}")
        message(FATAL_ERROR "multi-conformer JSON CLI output was not created: ${extension}")
    endif()
endforeach()

file(REMOVE "${multiconformer_input}")
file(REMOVE_RECURSE "${multiconformer_output_directory}")

set(coordinate_free_input "${CMAKE_CURRENT_BINARY_DIR}/coordinate_free.json")
set(coordinate_free_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_coordinate_free")
set(coordinate_free_output_prefix "${coordinate_free_output_directory}/coordinate_free.chargefw")
file(WRITE "${coordinate_free_input}" [=[
{"schema_version":"1.0","molecules":[{"atoms":[{"atomic_number":1,"formal_charge":0}]}]}
]=])
file(REMOVE_RECURSE "${coordinate_free_output_directory}")
file(MAKE_DIRECTORY "${coordinate_free_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method formal --output-mol2
                "${coordinate_free_input}"
                "${coordinate_free_output_directory}"
        RESULT_VARIABLE coordinate_free_result
        ERROR_VARIABLE coordinate_free_error
)
if(NOT coordinate_free_result EQUAL 6)
    message(FATAL_ERROR "coordinate-free MOL2 failure was not reported as an export error: ${coordinate_free_error}")
endif()
if(NOT EXISTS "${coordinate_free_output_prefix}.json")
    message(FATAL_ERROR "failed MOL2 generation did not write JSON")
endif()
file(REMOVE "${coordinate_free_input}")
file(REMOVE_RECURSE "${coordinate_free_output_directory}")

set(range_input "${CMAKE_CURRENT_BINARY_DIR}/out_of_range_charge.json")
set(range_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_out_of_range_charge")
set(range_output_prefix "${range_output_directory}/out_of_range_charge.chargefw")
file(WRITE "${range_input}" [=[
{"schema_version":"1.0","molecules":[{"atoms":[{"atomic_number":1,"formal_charge":6}],"conformers":[{"coordinates":[[0,0,0]]}]}]}
]=])
file(REMOVE_RECURSE "${range_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method formal --output-mol2 --output-mmcif
                "${range_input}"
                "${range_output_directory}"
        RESULT_VARIABLE range_result
        ERROR_VARIABLE range_error
)
if(NOT range_result EQUAL 6)
    message(FATAL_ERROR "mmCIF range failure was not reported as an export error: ${range_error}")
endif()
if(NOT EXISTS "${range_output_prefix}.json" OR
   NOT EXISTS "${range_output_prefix}.mol2")
    message(FATAL_ERROR "second export failure discarded an earlier output")
endif()
file(REMOVE "${range_input}")
file(REMOVE_RECURSE "${range_output_directory}")

set(no_plan_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_no_plan")
set(no_plan_output_prefix "${no_plan_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${no_plan_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method smpqeq "${CHARGEFW_INPUT}"
                "${no_plan_output_directory}"
        RESULT_VARIABLE no_plan_result
        ERROR_VARIABLE no_plan_error
)
if(NOT no_plan_result EQUAL 3)
    message(FATAL_ERROR "no-plan CLI exit status was ${no_plan_result}: ${no_plan_error}")
endif()
if(NOT EXISTS "${no_plan_output_prefix}.json")
    message(FATAL_ERROR "no-plan CLI JSON output was not written to the output directory")
endif()
foreach(extension IN ITEMS mol2 cif)
    if(EXISTS "${no_plan_output_prefix}.${extension}")
        message(FATAL_ERROR "no-plan CLI wrote unexpected molecular output: ${extension}")
    endif()
endforeach()
file(READ "${no_plan_output_prefix}.json" no_plan_json)
string(JSON no_plan_status GET "${no_plan_json}" status)
string(JSON no_plan_record_status GET "${no_plan_json}" results 0 status)
string(JSON no_plan_diagnostic GET "${no_plan_json}" diagnostics 0 code)
string(JSON no_plan_cause GET "${no_plan_json}" results 0 diagnostics 1 code)
if(NOT no_plan_status STREQUAL "no_executable_plan" OR
   NOT no_plan_record_status STREQUAL "no_executable_plan" OR
   NOT no_plan_diagnostic STREQUAL "no_executable_plan" OR
   no_plan_cause STREQUAL "no_executable_plan")
    message(FATAL_ERROR "no-plan CLI JSON does not report a structured no-plan result")
endif()
file(REMOVE_RECURSE "${no_plan_output_directory}")

set(warning_input_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_input_warning")
set(warning_input_output_prefix "${warning_input_output_directory}/aromatic.chargefw")
file(REMOVE_RECURSE "${warning_input_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method eem "${CHARGEFW_MOL2_INPUT}"
                "${warning_input_output_directory}"
        RESULT_VARIABLE warning_input_result
        ERROR_VARIABLE warning_input_error
)
if(NOT warning_input_result EQUAL 0)
    message(FATAL_ERROR "calculation with an import warning failed: ${warning_input_error}")
endif()
file(READ "${warning_input_output_prefix}.json" warning_input_json)
string(JSON warning_input_code GET "${warning_input_json}" results 0 diagnostics 0 code)
if(NOT warning_input_code STREQUAL "partial_charges_ignored")
    message(FATAL_ERROR "successful import warning was not retained in JSON")
endif()
file(REMOVE_RECURSE "${warning_input_output_directory}")

set(invalid_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_invalid")
set(invalid_output_prefix "${invalid_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${invalid_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --threads 18446744073709551615 "${CHARGEFW_INPUT}"
                "${invalid_output_directory}"
        RESULT_VARIABLE invalid_result
        ERROR_VARIABLE invalid_error
)
if(NOT invalid_result EQUAL 2)
    message(FATAL_ERROR "invalid-request CLI exit status was ${invalid_result}: ${invalid_error}")
endif()
if(NOT EXISTS "${invalid_output_prefix}.json")
    message(FATAL_ERROR "invalid-request CLI JSON output was not written to the output directory")
endif()
file(READ "${invalid_output_prefix}.json" invalid_json)
string(JSON invalid_status GET "${invalid_json}" status)
string(JSON invalid_diagnostic GET "${invalid_json}" diagnostics 0 code)
if(NOT invalid_status STREQUAL "invalid_input_or_request" OR
   NOT invalid_diagnostic STREQUAL "invalid_input_or_request")
    message(FATAL_ERROR "invalid-request CLI JSON does not report a structured invalid result")
endif()
file(REMOVE_RECURSE "${invalid_output_directory}")

execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate "${CHARGEFW_INPUT}" "${output_directory}"
        RESULT_VARIABLE result
        ERROR_VARIABLE error
)

if(NOT result EQUAL 0)
    message(FATAL_ERROR "chargefw failed with exit code ${result}: ${error}")
endif()

if(NOT IS_DIRECTORY "${output_directory}")
    message(FATAL_ERROR "Output directory was not created: ${output_directory}")
endif()

if(NOT EXISTS "${output_prefix}.json")
    message(FATAL_ERROR "Default JSON output was not created")
endif()
foreach(extension IN ITEMS mol2 cif)
    if(EXISTS "${output_prefix}.${extension}")
        message(FATAL_ERROR "Unrequested molecular output was created: ${extension}")
    endif()
endforeach()

file(READ "${output_prefix}.json" json_output)
string(JSON schema_version GET "${json_output}" schema_version)
if(NOT schema_version STREQUAL "1.0")
    message(FATAL_ERROR "Unexpected schema version: ${schema_version}")
endif()
string(JSON generator_version GET "${json_output}" generator version)
if(NOT generator_version STREQUAL "${CHARGEFW_VERSION}")
    message(FATAL_ERROR "Unexpected generator version: ${generator_version}")
endif()
string(JSON status GET "${json_output}" results 0 status)
if(NOT status STREQUAL "success")
    message(FATAL_ERROR "Expected successful calculation, got ${status}")
endif()
string(JSON requested_execution GET "${json_output}" calculation_provenance requested execution kind)
if(NOT requested_execution STREQUAL "auto")
    message(FATAL_ERROR "Unexpected requested execution kind: ${requested_execution}")
endif()
string(JSON execution_mode GET "${json_output}" calculation_provenance effective execution mode)
if(NOT execution_mode STREQUAL "full")
    message(FATAL_ERROR "Unexpected effective execution mode: ${execution_mode}")
endif()
string(JSON selected_method GET "${json_output}" calculation_provenance effective method id)
if(selected_method STREQUAL "")
    message(FATAL_ERROR "Selected method provenance is empty")
endif()
string(JSON permissive_types GET "${json_output}" calculation_provenance requested classification permissive_types)
if(permissive_types)
    message(FATAL_ERROR "Expected strict type classification by default")
endif()
string(JSON cutoff_atom_threshold GET "${json_output}" calculation_provenance requested resource_policy cutoff_atom_threshold)
if(NOT cutoff_atom_threshold STREQUAL "20000")
    message(FATAL_ERROR "Unexpected default cutoff atom threshold: ${cutoff_atom_threshold}")
endif()
string(JSON cover_atom_threshold GET "${json_output}" calculation_provenance requested resource_policy cover_atom_threshold)
if(NOT cover_atom_threshold STREQUAL "80000")
    message(FATAL_ERROR "Unexpected default cover atom threshold: ${cover_atom_threshold}")
endif()

set(explicit_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_explicit_outputs")
set(explicit_output_prefix "${explicit_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${explicit_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --output-mol2 --output-mmcif "${CHARGEFW_INPUT}"
                "${explicit_output_directory}"
        RESULT_VARIABLE explicit_result
        ERROR_VARIABLE explicit_error
)
if(NOT explicit_result EQUAL 0)
    message(FATAL_ERROR "explicit molecular output failed: ${explicit_error}")
endif()
foreach(extension IN ITEMS json mol2 cif)
    if(NOT EXISTS "${explicit_output_prefix}.${extension}")
        message(FATAL_ERROR "Requested output file was not created: ${extension}")
    endif()
endforeach()
file(REMOVE_RECURSE "${output_directory}")
file(REMOVE_RECURSE "${explicit_output_directory}")
