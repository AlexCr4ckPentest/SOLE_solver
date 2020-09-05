cmake_minimum_required(VERSION 3.1)

project(SOLE-solver)

set(CMAKE_CXX_STANDART 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -Wpedantic -O3 -Os")

file(GLOB SRC_FILES "src/*.cpp")
file(GLOB INC_FILES "include/*.hpp")

add_executable(${PROJECT_NAME} ${SRC_FILES} ${INC_FILES})
target_link_libraries(${PROJECT_NAME})