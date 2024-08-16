cmake_minimum_required(VERSION 3.20...3.24)

project(CycleLengthVerifier VERSION 1.0 LANGUAGES CXX)

find_package(Boost 1.85 REQUIRED COMPONENTS dynamic_bitset multiprecision)

add_executable(CycleLengthVerifier "main.cpp" "TransitionMatrix.cpp")

target_link_libraries(CycleLengthVerifier PRIVATE Boost::dynamic_bitset Boost::multiprecision)

target_compile_features(CycleLengthVerifier PUBLIC cxx_std_20)
set_target_properties(CycleLengthVerifier PROPERTIES CXX_EXTENSIONS OFF)