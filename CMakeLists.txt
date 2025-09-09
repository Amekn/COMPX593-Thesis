cmake_minimum_required(VERSION 4.0)
project(OntTools VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(VariantGenerator src/VariantGenerator.cpp)
add_executable(UmiFilter src/UmiFilter.cpp)
add_executable(LengthFilter src/LengthFilter.cpp)
add_executable(UmiExtractor.cpp src/UmiExtractor.cpp)
