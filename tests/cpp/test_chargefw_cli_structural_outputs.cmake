function(run_structural_input extension contents input_stem expected_selection expected_bonds)
    set(input_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${input_stem}.${extension}")
    set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_${input_stem}_outputs")
    set(output_prefix "${output_directory}/chargefw_cli_${input_stem}.chargefw")
    file(WRITE "${input_path}" "${contents}")
    file(REMOVE_RECURSE "${output_directory}")

    execute_process(
            COMMAND "${CHARGEFW_CLI}" calculate ${ARGN} "${input_path}" "${output_directory}"
            RESULT_VARIABLE result
            OUTPUT_VARIABLE output
            ERROR_VARIABLE error
    )

    if(NOT result EQUAL 0)
        message(FATAL_ERROR "chargefw failed for .${extension} with exit code ${result}: ${error}")
    endif()

    if(NOT EXISTS "${output_prefix}.json")
        message(FATAL_ERROR "JSON output was not created for .${extension}")
    endif()
    if(NOT EXISTS "${output_prefix}.cif")
        message(FATAL_ERROR "mmCIF output was not created for .${extension}")
    endif()
    if(EXISTS "${output_prefix}.sdf" OR EXISTS "${output_prefix}.mol2")
        message(FATAL_ERROR "Structural input must not produce SDF or MOL2 output")
    endif()

    file(READ "${output_prefix}.json" json_output)
    string(JSON status GET "${json_output}" results 0 status)
    if(NOT status STREQUAL "success")
        message(FATAL_ERROR "Expected successful .${extension} calculation, got ${status}")
    endif()
    string(JSON structural_selection GET "${json_output}" calculation_provenance requested structural_input selection)
    string(JSON structural_bonds GET "${json_output}" calculation_provenance requested structural_input bonds)
    if(NOT structural_selection STREQUAL "${expected_selection}" OR
       NOT structural_bonds STREQUAL "${expected_bonds}")
        message(FATAL_ERROR "Unexpected structural provenance: selection=${structural_selection}, bonds=${structural_bonds}")
    endif()

    file(REMOVE "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

function(expect_non_structural_options_rejected)
    set(input_path "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_nonstructural.json")
    set(output_directory "${CMAKE_CURRENT_BINARY_DIR}/chargefw_cli_nonstructural_outputs")
    file(WRITE "${input_path}" "{\"schema_version\": \"1.0\", \"molecules\": [{\"atoms\": [{\"atomic_number\": 8, \"formal_charge\": 0}]}]}\n")
    file(REMOVE_RECURSE "${output_directory}")

    execute_process(
            COMMAND "${CHARGEFW_CLI}" calculate --structural-bonds templates "${input_path}" "${output_directory}"
            RESULT_VARIABLE result
            ERROR_VARIABLE error
    )

    if(result EQUAL 0 OR NOT error MATCHES "Structural input options are only supported for PDB and mmCIF input")
        message(FATAL_ERROR "Expected structural options to be rejected for JSON input: ${error}")
    endif()

    file(REMOVE "${input_path}")
    file(REMOVE_RECURSE "${output_directory}")
endfunction()

run_structural_input(
        pdb
        "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  \nATOM      2  H1  HOH A   1       0.957   0.000   0.000  1.00 20.00           H  \nATOM      3  H2  HOH A   1      -0.239   0.927   0.000  1.00 20.00           H  \nEND\n"
        structural_pdb
        polymers-and-ligands
        templates
        --structural-selection polymers-and-ligands
        --structural-bonds templates
)

run_structural_input(
        cif
        "data_structural_cif\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\nHETATM 1 O O . HOH A 1 ? 0.000 0.000 0.000 1.00 20.00 0 1 HOH A O 1\nHETATM 2 H H1 . HOH A 1 ? 0.957 0.000 0.000 1.00 20.00 0 1 HOH A H1 1\nHETATM 3 H H2 . HOH A 1 ? -0.239 0.927 0.000 1.00 20.00 0 1 HOH A H2 1\n#\n"
        structural_cif
        all
        hybrid
        --structural-selection all
        --structural-bonds hybrid
)

run_structural_input(
        cif
        "data_structural_default_bonds\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\nHETATM 1 O O . HOH A 1 ? 0.000 0.000 0.000 1.00 20.00 0 1 HOH A O 1\nHETATM 2 H H1 . HOH A 1 ? 0.957 0.000 0.000 1.00 20.00 0 1 HOH A H1 1\nHETATM 3 H H2 . HOH A 1 ? -0.239 0.927 0.000 1.00 20.00 0 1 HOH A H2 1\n#\n"
        structural_default_bonds
        all
        hybrid
)

expect_non_structural_options_rejected()
