target_sources(chargefw_core
        PRIVATE
        calculation.cpp
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
        FILES
        ${PROJECT_SOURCE_DIR}/include/chargefw/calculation/calculation.h
)
