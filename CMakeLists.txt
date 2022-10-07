cmake_minimum_required(VERSION 3.8)
project(lightning_demo)
set (CMAKE_CXX_STANDARD 20)

add_subdirectory(libs/sokol)
add_subdirectory(libs/glm)
add_subdirectory(libs/imgui)

if(WIN32)
    set(SHDC "${CMAKE_SOURCE_DIR}/bin/sokol-tools/bin/win32/sokol-shdc.exe")
elseif(APPLE)
    set(SHDC "${CMAKE_SOURCE_DIR}/bin/sokol-tools/bin/osx/sokol-shdc")
else()
    set(SHDC "${CMAKE_SOURCE_DIR}/bin/sokol-tools/bin/linux/sokol-shdc")
endif()

file(GLOB_RECURSE glsl_files "${CMAKE_SOURCE_DIR}/src/*.glsl")
foreach(glsl_file ${glsl_files})
    string(REGEX REPLACE "[.]glsl$" ".shader.h" shader_header_file ${glsl_file})
    #message(${shader_header_file})
    add_custom_command(
        DEPENDS ${glsl_file}
        OUTPUT ${shader_header_file}
        COMMAND ${SHDC} --slang=glsl330 --input=\"${glsl_file}\" --output=\"${shader_header_file}\"
    )
    list(APPEND shader_header_files ${shader_header_file})
endforeach()
add_custom_target(shaders_target DEPENDS ${shader_header_files})

set(SRCS src/main.cpp)
add_executable(lightning_demo ${SRCS})
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lightning_demo)

target_link_libraries(lightning_demo sokol glm imgui)
add_dependencies(lightning_demo shaders_target)

set_property(TARGET lightning_demo PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

source_group("" FILES ${SRCS})