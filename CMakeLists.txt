cmake_minimum_required(VERSION 2.4 FATAL_ERROR)

project(HTree_TO_DT)

set(CMAKE_BUILD_TYPE debug)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{DATATREE_HOME})

find_package(ROOT REQUIRED COMPONENTS RIO Tree)
find_package(DataTree REQUIRED)

set(CMAKE_CXX_FLAGS "-Wall")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_PREFIX_PATH $ENV{DATATREE})

#include(${ROOT_USE_FILE})

include_directories(
    ${CMAKE_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${DataTree_INCLUDE_DIR}/src
    )
add_definitions(${ROOT_CXX_FLAGS})

add_executable(HTree_to_DT HADES_constants.h main.cpp)
target_link_libraries(HTree_to_DT
    ${ROOT_LIBRARIES}
		DataTree
    )
