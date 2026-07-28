target_sources(chargefw_core
        PRIVATE
        molecule_record.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/adapters/molecule_record.h
)
