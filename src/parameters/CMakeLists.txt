target_sources(chargefw_core
        PRIVATE
        atom_parameters.cpp
        bond_parameters.cpp
        common_parameters.cpp
        parameter_set.cpp
        parameter_helpers.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/parameters/parameter_set_metadata.h
)