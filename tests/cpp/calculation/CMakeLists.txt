chargefw_add_test(test_calculation
        test_calculation.cpp
)

chargefw_add_test(test_observer
        test_observer.cpp
)

chargefw_add_test(test_execution_policy
        test_execution_policy.cpp
)

chargefw_add_test(test_target_execution
        test_target_execution.cpp
)
target_include_directories(test_target_execution PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)
target_link_libraries(test_target_execution PRIVATE
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
