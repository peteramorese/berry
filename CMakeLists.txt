cmake_minimum_required(VERSION 3.16)

project(berry)

find_package (Eigen3 3.3.8 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)

if(NOT DEFINED BRY_BUILD_EXECUTABLES)
    option(BRY_BUILD_EXECUTABLES "Build executables (set to OFF for building just the library target)" ON)
endif()

if(NOT DEFINED BRY_BUILD_YAML)
    option(BRY_BUILD_YAML "Build the yaml-cpp dependency" ON)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Dependencies
#if(BRY_BUILD_YAML)
#    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/dependencies/yaml-cpp)
#endif()

set(BRY_DEPENDENCY_INCLUDE_DIRS 
    #${CMAKE_CURRENT_SOURCE_DIR}/dependencies/yaml-cpp/include
    ${EIGEN3_INCLUDE_DIRS}
)
set(BRY_LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/berry)

#add_subdirectory(berry)

#if(BRY_BUILD_EXECUTABLES)
    add_subdirectory(src)
#endif()