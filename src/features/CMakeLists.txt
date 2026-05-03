target_sources(chargefw_core
        PRIVATE
        topology_features.cpp
        conformer_features.cpp
        prepared_molecule.cpp
        prepared_molecule_collection.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/features/topology_features.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/features/conformer_features.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/features/prepared_molecule.h
        ${PROJECT_SOURCE_DIR}/include/chargefw/features/prepared_molecule_collection.h
)