cmake_minimum_required(VERSION 3.16)

project(berry)

find_package (Eigen3 3.3.8 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)

# TODO: Optionfy this
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

if(NOT DEFINED BRY_BUILD_EXECUTABLES)
    option(BRY_BUILD_EXECUTABLES "Build executables (ON by default, set to OFF for building just the library target)" ON)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

set(LMN_BUILD_EXECUTABLES OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modules/lemon)

set(BRY_DEPENDENCY_INCLUDE_DIRS 
    ${EIGEN3_INCLUDE_DIRS}
    ${LMN_INCLUDE_DIRS}
)

set(BRY_INCLUDE_DIRS 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/berry
    ${CMAKE_CURRENT_SOURCE_DIR}/berry/impl
    ${BRY_DEPENDENCY_INCLUDE_DIRS}
    CACHE INTERNAL ""
)

if(BRY_BUILD_EXECUTABLES)
   add_subdirectory(src)
endif()