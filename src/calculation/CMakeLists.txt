target_sources(chargefw_core
        PRIVATE
        calculation.cpp
        cutoff_execution.cpp
        execution_policy.cpp
        parallel_for.h
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
         FILES
         ${PROJECT_SOURCE_DIR}/include/chargefw/calculation/calculation.h
         ${PROJECT_SOURCE_DIR}/include/chargefw/calculation/execution_policy.h
)
