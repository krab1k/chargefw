cmake_minimum_required(VERSION 3.28)

project(
  chargefw
  VERSION 0.1
  DESCRIPTION "Empirical partial atomic charge calculation framework"
  LANGUAGES CXX
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

file(RELATIVE_PATH CHARGEFW_PARAMETER_DIRECTORY_FROM_LIBRARY
  "${CMAKE_INSTALL_FULL_LIBDIR}"
  "${CMAKE_INSTALL_FULL_DATADIR}/chargefw/parameters"
)

option(CHARGEFW_ENABLE_CCACHE "Use ccache if available" ON)
option(CHARGEFW_BUILD_TESTS "Build tests" ON)
option(CHARGEFW_BUILD_CLI "Build the chargefw command-line application" ON)
option(CHARGEFW_ENABLE_WARNINGS "Enable project warnings" ON)
option(CHARGEFW_ENABLE_ASAN "Enable AddressSanitizer" OFF)
option(CHARGEFW_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)
option(CHARGEFW_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
option(CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS "Optimize for the build host CPU" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

include(cmake/dependencies.cmake)
include(cmake/developer_tools.cmake)
include(cmake/optimization.cmake)
include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)

chargefw_setup_dependencies()
chargefw_setup_developer_tools()

add_library(chargefw_core SHARED)
add_library(chargefw::core ALIAS chargefw_core)
set_target_properties(chargefw_core PROPERTIES EXPORT_NAME core)

target_compile_features(chargefw_core
  PUBLIC
    cxx_std_23
)

target_include_directories(chargefw_core
  PRIVATE
    src

  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/generated>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(chargefw_core SYSTEM PRIVATE
  $<TARGET_PROPERTY:TBB::tbb,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(chargefw_core
  PUBLIC
    nlohmann_json::nlohmann_json
    gemmi::gemmi_cpp
  PRIVATE
    $<BUILD_INTERFACE:Eigen3::Eigen>
    $<BUILD_INTERFACE:nanoflann::nanoflann>
    $<BUILD_INTERFACE:$<TARGET_FILE:TBB::tbb>>
    $<BUILD_INTERFACE:${CMAKE_DL_LIBS}>
)

if(TARGET tbb)
  add_dependencies(chargefw_core tbb)
endif()

target_compile_definitions(chargefw_core
  PRIVATE
    EIGEN_DONT_PARALLELIZE
)

chargefw_enable_optimizations(chargefw_core)

configure_file(
  ${PROJECT_SOURCE_DIR}/include/chargefw/config.h.in
  ${PROJECT_BINARY_DIR}/generated/chargefw/config.h
  @ONLY
)

if(CHARGEFW_ENABLE_WARNINGS)
  chargefw_enable_warnings(chargefw_core)
endif()

if(CHARGEFW_ENABLE_ASAN)
  chargefw_enable_sanitizer(chargefw_core address)
endif()

if(CHARGEFW_ENABLE_UBSAN)
  chargefw_enable_sanitizer(chargefw_core undefined)
endif()

add_subdirectory(src)

if(CHARGEFW_BUILD_CLI)
  add_subdirectory(apps/chargefw)
endif()

if(CHARGEFW_BUILD_TESTS)
  chargefw_setup_test_dependencies()
  enable_testing()
  add_subdirectory(tests)
endif()

include(cmake/install_rules.cmake)
