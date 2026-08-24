chargefw_add_test(test_calculation
        test_calculation.cpp
)

chargefw_add_test(test_observer
        test_observer.cpp
)

chargefw_add_test(test_execution_policy
        test_execution_policy.cpp
)

chargefw_add_test(test_calculation_targets
        test_calculation_targets.cpp
)
target_include_directories(test_calculation_targets PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(test_calculation_targets PRIVATE
        TBB::tbb
)

chargefw_add_test(test_cutoff_execution
        test_cutoff_execution.cpp
)

chargefw_add_test(test_cover_execution
        test_cover_execution.cpp
)
target_include_directories(test_cover_execution PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
