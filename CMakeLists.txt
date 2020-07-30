# To cmake project, from project root do:
# 1. cd build
# 2. cmake ..
# 3. cmake --build .
# 4. ctest -VV

# Basic project information
cmake_minimum_required(VERSION 3.6)
set(PROJECT_NAME_STR RangeTree)
project(${PROJECT_NAME_STR} C CXX)
set (RangeTree_VERSION_MAJOR 0)
set (RangeTree_VERSION_MINOR 1)

find_package(Threads REQUIRED)

if(CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-Wall -ansi -Wno-deprecated)
endif()

# Flags
set(CMAKE_CXX_STANDARD 11)

#-------------------
# set common include folder for module
#-------------------
set(COMMON_INCLUDES ./)
set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

#-------------------
# Module source
#-------------------
#set(PROJECT_LIB_NAME ${PROJECT_NAME_STR}_impl)
#
#include_directories(${COMMON_INCLUDES})
#file(GLOB SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
#add_library(${PROJECT_LIB_NAME} ${SRC_FILES})
#
add_subdirectory(${EXT_PROJECTS_DIR}/gtest)

#-------------------
# Test
#-------------------
add_definitions(${MSVC_COMPILER_DEFS})
enable_testing()
set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_test)
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES})

file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)
add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)

#target_link_libraries(${PROJECT_TEST_NAME} ${PROJECT_LIB_NAME})

target_link_libraries(${PROJECT_TEST_NAME}
        ${GTEST_LIBS_DIR}/libgtest.a
        ${GTEST_LIBS_DIR}/libgtest_main.a
        )

target_link_libraries(${PROJECT_TEST_NAME} ${CMAKE_THREAD_LIBS_INIT})

add_test(RangeTreeTests ${PROJECT_TEST_NAME})
