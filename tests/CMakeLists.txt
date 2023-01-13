cmake_minimum_required(VERSION 3.14)

set(TEST_PROJECT_NAME ${LIB_TARGET_NAME}_tests)
set(TEST_TARGET_NAME ${TEST_PROJECT_NAME})

project(${TEST_PROJECT_NAME})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.12.1
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(
    ${TEST_PROJECT_NAME}
    remez_test.cpp
)

set(GTEST_INCLUDE_DIR "${CMAKE_BINARY_DIR}/_deps/googletest-src/googletest/include")

target_include_directories(
    ${TEST_PROJECT_NAME}
    SYSTEM BEFORE PRIVATE
    "${GTEST_INCLUDE_DIR}"
    "${LIB_DEV_PUBLIC_HEADER_DIR}"
)

target_link_libraries(
    ${TEST_PROJECT_NAME}
    ${LIB_TARGET_NAME}
    gtest_main
)

add_test(NAME ${TEST_PROJECT_NAME}_tests1 COMMAND ${TEST_PROJECT_NAME})
