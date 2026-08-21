target_sources(chargefw_core
        PRIVATE
        calculation.cpp
        cover_execution.cpp
        cutoff_execution.cpp
        execution_policy.cpp
        reduced_execution.cpp
        cover_execution.h
        parallel_for.h
        reduced_execution.h
        PUBLIC
        FILE_SET HEADERS
        BASE_DIRS
        ${PROJECT_SOURCE_DIR}/include
         FILES
         ${PROJECT_SOURCE_DIR}/include/chargefw/calculation/calculation.h
         ${PROJECT_SOURCE_DIR}/include/chargefw/calculation/execution_policy.h
)
