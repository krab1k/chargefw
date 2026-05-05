include(FetchContent)

function(chargefw_setup_dependencies)
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
endfunction()