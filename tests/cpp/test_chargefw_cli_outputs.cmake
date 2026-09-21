set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_outputs")
set(output_prefix "${output_directory}/water.chargefw")
include("${CMAKE_CURRENT_LIST_DIR}/prepare_moved_install.cmake")
file(REMOVE_RECURSE "${output_directory}")

set(json_failure_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_json_failure")
set(json_failure_output_prefix "${json_failure_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${json_failure_output_directory}")
file(MAKE_DIRECTORY "${json_failure_output_prefix}.json")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method eem "${CHARGEFW_INPUT}"
                "${json_failure_output_directory}"
        RESULT_VARIABLE json_failure_result
        OUTPUT_VARIABLE json_failure_output
        ERROR_VARIABLE json_failure_error
)
file(GLOB json_failure_temporary_files "${json_failure_output_prefix}.json.tmp.*")
if(NOT json_failure_result EQUAL 2 OR
   NOT json_failure_error MATCHES "Unable to publish output file: ${json_failure_output_prefix}.json" OR
   json_failure_output MATCHES "Wrote" OR json_failure_temporary_files)
    message(FATAL_ERROR "atomic JSON publication failure was not reported correctly: ${json_failure_error}")
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
        OUTPUT_VARIABLE export_failure_output
        ERROR_VARIABLE export_failure_error
)
file(GLOB export_failure_temporary_files "${export_failure_output_prefix}.mol2.tmp.*")
if(NOT export_failure_result EQUAL 6 OR
   NOT export_failure_error MATCHES "Export error: Unable to publish output file: ${export_failure_output_prefix}.mol2" OR
   NOT EXISTS "${export_failure_output_prefix}.json" OR
   NOT export_failure_output MATCHES "Wrote ${export_failure_output_prefix}.json" OR
   export_failure_temporary_files)
    message(FATAL_ERROR "atomic molecular export failure was not reported correctly: ${export_failure_error}")
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
    if(mode STREQUAL "cutoff" AND NOT mode_error MATCHES "Calculating eem")
        message(FATAL_ERROR "cutoff CLI progress output was not emitted: ${mode_error}")
    endif()
    if(mode STREQUAL "cutoff")
        string(ASCII 27 erase_escape)
        set(clear_line "\r${erase_escape}[2K")
        string(REPLACE "${clear_line}" "" progress_without_clear_sequences "${mode_error}")
        if(progress_without_clear_sequences MATCHES "\r")
            message(FATAL_ERROR "progress output did not erase a previous rendered line: ${mode_error}")
        endif()
        string(FIND "${mode_error}" "${clear_line}\n" final_clear_position)
        if(NOT mode_error MATCHES "Fragments" OR NOT mode_error MATCHES "Targets" OR
           final_clear_position EQUAL -1)
            message(FATAL_ERROR "progress output did not clear its final rendered line: ${mode_error}")
        endif()
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

file(READ "${multiconformer_output_prefix}.cif" multiconformer_cif)
string(FIND "${multiconformer_cif}" "1.25 2.5 3.75" cif_first_position)
string(FIND "${multiconformer_cif}" "9.5 8.5 7.5" cif_second_position)
if(cif_first_position EQUAL -1 OR cif_second_position EQUAL -1)
    message(FATAL_ERROR "generated mmCIF did not retain every conformer")
endif()
file(READ "${multiconformer_output_prefix}.mol2" multiconformer_mol2)
string(FIND "${multiconformer_mol2}" "1.25 2.5 3.75" mol2_first_position)
string(FIND "${multiconformer_mol2}" "9.5 8.5 7.5" mol2_second_position)
if(mol2_first_position EQUAL -1 OR mol2_second_position EQUAL -1)
    message(FATAL_ERROR "generated MOL2 did not retain every conformer")
endif()
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
file(WRITE "${coordinate_free_output_prefix}.mol2" "existing export\n")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method formal --output-mol2
                "${coordinate_free_input}"
                "${coordinate_free_output_directory}"
        RESULT_VARIABLE coordinate_free_result
        ERROR_VARIABLE coordinate_free_error
)
if(NOT coordinate_free_result EQUAL 6 OR
   NOT coordinate_free_error MATCHES "Export error: MOL2 output requires coordinates")
    message(FATAL_ERROR "coordinate-free MOL2 failure was not reported as an export error: ${coordinate_free_error}")
endif()
file(READ "${coordinate_free_output_prefix}.mol2" coordinate_free_mol2)
file(GLOB coordinate_free_temporary_files "${coordinate_free_output_prefix}.mol2.tmp.*")
if(NOT coordinate_free_mol2 STREQUAL "existing export\n" OR
   NOT EXISTS "${coordinate_free_output_prefix}.json" OR coordinate_free_temporary_files)
    message(FATAL_ERROR "failed MOL2 generation changed a published file or left a temporary file")
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
if(NOT range_result EQUAL 6 OR
   NOT range_error MATCHES "Export error: mmCIF charge is outside the dictionary range")
    message(FATAL_ERROR "mmCIF range failure was not reported as an export error: ${range_error}")
endif()
file(GLOB range_temporary_files "${range_output_prefix}.*.tmp.*")
if(NOT EXISTS "${range_output_prefix}.json" OR
   NOT EXISTS "${range_output_prefix}.mol2" OR
   EXISTS "${range_output_prefix}.cif" OR range_temporary_files)
    message(FATAL_ERROR "second export failure discarded an earlier output or left a temporary file")
endif()
file(REMOVE "${range_input}")
file(REMOVE_RECURSE "${range_output_directory}")

set(deterministic_result "")
foreach(run IN ITEMS 1 2)
    set(deterministic_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_deterministic_${run}")
    set(deterministic_output_prefix "${deterministic_output_directory}/water.chargefw")
    file(REMOVE_RECURSE "${deterministic_output_directory}")
    execute_process(
            COMMAND "${CHARGEFW_CLI}" calculate --method eem --execution full "${CHARGEFW_INPUT}"
                    "${deterministic_output_directory}"
            RESULT_VARIABLE deterministic_run_result
            ERROR_VARIABLE deterministic_error
    )
    if(NOT deterministic_run_result EQUAL 0)
        message(FATAL_ERROR "deterministic CLI calculation failed: ${deterministic_error}")
    endif()
    file(READ "${deterministic_output_prefix}.json" deterministic_json)
    string(JSON deterministic_method GET "${deterministic_json}" calculation_provenance effective method id)
    string(JSON deterministic_mode GET "${deterministic_json}" calculation_provenance effective execution mode)
    string(JSON deterministic_charge_0 GET "${deterministic_json}" results 0 assignments 0 charges 0)
    string(JSON deterministic_charge_1 GET "${deterministic_json}" results 0 assignments 0 charges 1)
    string(JSON deterministic_charge_2 GET "${deterministic_json}" results 0 assignments 0 charges 2)
    set(current_deterministic_result
        "${deterministic_method}|${deterministic_mode}|${deterministic_charge_0}|${deterministic_charge_1}|${deterministic_charge_2}")
    if(deterministic_result STREQUAL "")
        set(deterministic_result "${current_deterministic_result}")
    elseif(NOT deterministic_result STREQUAL current_deterministic_result)
        message(FATAL_ERROR "CLI stable result fields differ between identical runs")
    endif()
    file(REMOVE_RECURSE "${deterministic_output_directory}")
endforeach()

set(no_plan_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_no_plan")
set(no_plan_output_prefix "${no_plan_output_directory}/water.chargefw")
file(REMOVE_RECURSE "${no_plan_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method smpqeq "${CHARGEFW_INPUT}"
                "${no_plan_output_directory}"
        RESULT_VARIABLE no_plan_result
        OUTPUT_VARIABLE no_plan_output
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
if(NOT no_plan_error MATCHES "Error: method 'smpqeq'")
    message(FATAL_ERROR "no-plan CLI did not print its applicability cause: ${no_plan_error}")
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
if(NOT warning_input_result EQUAL 0 OR NOT warning_input_error MATCHES "Warning: MOL2 partial charges were ignored")
    message(FATAL_ERROR "successful import warning was not reported: ${warning_input_error}")
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

set(warning_output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_warning")
file(REMOVE_RECURSE "${warning_output_directory}")
execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate --method eem --execution full --cutoff-atom-threshold 0
                --cover-atom-threshold 0
                "${CHARGEFW_INPUT}" "${warning_output_directory}"
        RESULT_VARIABLE warning_result
        ERROR_VARIABLE warning_error
)
if(NOT warning_result EQUAL 0 OR NOT warning_error MATCHES "Warning: ")
    message(FATAL_ERROR "Explicit-full warning was not emitted: ${warning_error}")
endif()
file(REMOVE_RECURSE "${warning_output_directory}")

execute_process(
        COMMAND "${CHARGEFW_CLI}" calculate "${CHARGEFW_INPUT}" "${output_directory}"
        RESULT_VARIABLE result
        OUTPUT_VARIABLE output
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
file(READ "${explicit_output_prefix}.mol2" mol2_output)
if(NOT mol2_output MATCHES "@<TRIPOS>MOLECULE" OR
   NOT mol2_output MATCHES "USER_CHARGES")
    message(FATAL_ERROR "Generated MOL2 output is incomplete")
endif()

file(REMOVE_RECURSE "${output_directory}")
file(REMOVE_RECURSE "${explicit_output_directory}")

execute_process(
        COMMAND "${CHARGEFW_CLI}" inspect "${CHARGEFW_INPUT}"
        RESULT_VARIABLE inspect_result
        OUTPUT_VARIABLE inspect_output
)
if(NOT inspect_result EQUAL 0 OR NOT inspect_output MATCHES "records: 1" OR
   NOT inspect_output MATCHES "atoms=3")
    message(FATAL_ERROR "inspect output is incomplete: ${inspect_output}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --method formal "${CHARGEFW_INPUT}"
        RESULT_VARIABLE applicability_result
        OUTPUT_VARIABLE applicability_output
)
if(NOT applicability_result EQUAL 0 OR NOT applicability_output MATCHES "plan method=formal")
    message(FATAL_ERROR "applicability output is incomplete: ${applicability_output}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" applicability --method smpqeq "${CHARGEFW_INPUT}"
        RESULT_VARIABLE rejected_applicability_result
        OUTPUT_VARIABLE rejected_applicability_output
)
if(NOT rejected_applicability_result EQUAL 0 OR
   NOT rejected_applicability_output MATCHES "rejected method=smpqeq")
    message(FATAL_ERROR "rejected applicability output has incorrect method identity: ${rejected_applicability_output}")
endif()

execute_process(COMMAND "${CHARGEFW_CLI}" methods OUTPUT_VARIABLE methods_output)
if(NOT methods_output MATCHES "formal.*Formal atomic charges" OR
   methods_output MATCHES "default=")
    message(FATAL_ERROR "method summary output is incorrect: ${methods_output}")
endif()

execute_process(COMMAND "${CHARGEFW_CLI}" methods eem OUTPUT_VARIABLE eem_method_output)
string(FIND "${eem_method_output}" "time complexity: O(n^3)" eem_time_complexity)
string(FIND "${eem_method_output}" "memory complexity: O(n^2)" eem_memory_complexity)
if(NOT eem_method_output MATCHES "full name: Electronegativity Equalization Method" OR
   NOT eem_method_output MATCHES "requires coordinates: yes" OR
   eem_time_complexity EQUAL -1 OR eem_memory_complexity EQUAL -1 OR
   NOT eem_method_output MATCHES "supports cutoff: yes" OR
   NOT eem_method_output MATCHES "supports cover: yes")
    message(FATAL_ERROR "EEM method details are incomplete: ${eem_method_output}")
endif()

execute_process(COMMAND "${CHARGEFW_CLI}" methods denr OUTPUT_VARIABLE denr_method_output)
if(NOT denr_method_output MATCHES "notes: This implementation fixes initial charges to zero" OR
   NOT denr_method_output MATCHES "minimum>=0" OR
   NOT denr_method_output MATCHES "minimum>0")
    message(FATAL_ERROR "method option details are incomplete: ${denr_method_output}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" parameters
        OUTPUT_VARIABLE parameters_output
)
if(NOT parameters_output MATCHES "QEq_original.*Rappe 1991")
    message(FATAL_ERROR "parameter summary output is incorrect: ${parameters_output}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" parameters --method qeq
        OUTPUT_VARIABLE filtered_parameters_output
)
if(NOT filtered_parameters_output MATCHES "QEq_original" OR
   filtered_parameters_output MATCHES "EEM_Baek1991")
    message(FATAL_ERROR "parameter filtering is incorrect: ${filtered_parameters_output}")
endif()

execute_process(
        COMMAND "${CHARGEFW_CLI}" parameters QEq_original
        OUTPUT_VARIABLE parameter_details_output
)
if(NOT parameter_details_output MATCHES "id: QEq_original" OR
   NOT parameter_details_output MATCHES "method: qeq" OR
   NOT parameter_details_output MATCHES "name: Rappe 1991" OR
   NOT parameter_details_output MATCHES "publication: 10.1021/j100161a070" OR
   NOT parameter_details_output MATCHES "notes: Derived from experimental atomic IPs and EAs" OR
   NOT parameter_details_output MATCHES "priority: 0")
    message(FATAL_ERROR "parameter-set details are incomplete: ${parameter_details_output}")
endif()
