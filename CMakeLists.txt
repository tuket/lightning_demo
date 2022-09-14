cmake_minimum_required(VERSION 3.8)
project(lightning_demo)
set (CMAKE_CXX_STANDARD 20)

add_subdirectory(libs/sokol)

set(SRCS src/main.cpp)
add_executable(lightning_demo ${SRCS})
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lightning_demo)

target_link_libraries(lightning_demo sokol)

set_property(TARGET lightning_demo PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

source_group("" FILES ${SRCS})