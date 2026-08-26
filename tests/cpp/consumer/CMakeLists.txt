cmake_minimum_required(VERSION 3.28)

project(chargefw_consumer LANGUAGES CXX)

find_package(chargefw CONFIG REQUIRED)

add_executable(chargefw_consumer main.cpp)
target_link_libraries(chargefw_consumer PRIVATE chargefw::core)
