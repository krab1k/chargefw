include(FetchContent)

function(chargefw_setup_dependencies)
  FetchContent_Declare(
    nlohmann_json
    SYSTEM
    URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
  )

  FetchContent_MakeAvailable(nlohmann_json)
endfunction()
