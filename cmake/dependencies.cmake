include(FetchContent)

function(chargefw_setup_dependencies)
    find_package(CLI11 2.7.2 CONFIG QUIET)
    if(NOT TARGET CLI11::CLI11)
        FetchContent_Declare(
                cli11
                SYSTEM
                URL https://github.com/CLIUtils/CLI11/archive/refs/tags/v2.7.2.tar.gz
        )
        FetchContent_MakeAvailable(cli11)
    endif()

    find_package(nlohmann_json 3.12 CONFIG QUIET)
    if(NOT TARGET nlohmann_json::nlohmann_json)
        FetchContent_Declare(
                nlohmann_json
                SYSTEM
                URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
        )
        FetchContent_MakeAvailable(nlohmann_json)
    endif()

    find_package(Eigen3 5.0 CONFIG QUIET)
    if(NOT TARGET Eigen3::Eigen)
        FetchContent_Declare(
                eigen
                SYSTEM
                URL https://gitlab.com/libeigen/eigen/-/archive/5.0.1/eigen-5.0.1.tar.gz
        )
        FetchContent_MakeAvailable(eigen)
    endif()

    find_package(nanoflann 1.12 CONFIG QUIET)
    if(NOT TARGET nanoflann::nanoflann)
        set(NANOFLANN_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
        set(NANOFLANN_BUILD_TESTS OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(
                nanoflann
                SYSTEM
                URL https://github.com/jlblancoc/nanoflann/archive/refs/tags/1.12.1.tar.gz
        )
        FetchContent_MakeAvailable(nanoflann)
    endif()

    find_package(TBB 2023.1 CONFIG QUIET)
    if(NOT TARGET TBB::tbb)
        set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
        set(TBB_TEST OFF CACHE BOOL "" FORCE)
        set(TBB_STRICT OFF CACHE BOOL "" FORCE)
        set(CHARGEFW_TBB_FETCHED ON CACHE INTERNAL "oneTBB was fetched by ChargeFW")
        FetchContent_Declare(
                onetbb
                SYSTEM
                URL https://github.com/uxlfoundation/oneTBB/archive/refs/tags/v2023.1.0.tar.gz
        )
        FetchContent_MakeAvailable(onetbb)
    endif()

    find_package(gemmi 0.7.4 CONFIG QUIET)
    if(NOT TARGET gemmi::gemmi_cpp)
        set(BUILD_GEMMI_PROGRAM OFF CACHE BOOL "" FORCE)
        set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)
        set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(GEMMI_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(CHARGEFW_GEMMI_FETCHED ON CACHE INTERNAL "Gemmi was fetched by ChargeFW")
        FetchContent_Declare(
                gemmi
                SYSTEM
                URL https://github.com/project-gemmi/gemmi/archive/refs/tags/v0.7.4.tar.gz
        )
        FetchContent_MakeAvailable(gemmi)
    endif()
endfunction()

function(chargefw_setup_test_dependencies)
    find_package(snitch 1.3.2 CONFIG QUIET)
    if(NOT TARGET snitch::snitch)
        FetchContent_Declare(
                snitch
                SYSTEM
                URL https://github.com/cschreib/snitch/archive/refs/tags/v1.3.2.tar.gz
        )
        FetchContent_MakeAvailable(snitch)
    endif()
endfunction()
