cmake_minimum_required(VERSION 3.8)
project(lightning_demo)
set (CMAKE_CXX_STANDARD 20)

add_subdirectory(libs/sokol)
add_subdirectory(libs/glm)
add_subdirectory(libs/imgui)

if(CMAKE_HOST_WIN32)
    set(SHDC "${CMAKE_SOURCE_DIR}/bin/sokol-tools/bin/win32/sokol-shdc.exe")
elseif(CMAKE_HOST_APPLE)
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
        COMMAND ${SHDC} --slang=glsl330:glsl300es:glsl100 --input=\"${glsl_file}\" --output=\"${shader_header_file}\"
    )
    list(APPEND shader_header_files ${shader_header_file})
endforeach()
add_custom_target(shaders_target DEPENDS ${shader_header_files})

set(SRCS
    src/main.cpp
    src/lightning.hpp src/lightning.cpp
)
add_executable(lightning_demo ${SRCS})
set_property(DIRECTORY ${PROJECT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lightning_demo)
if(EMSCRIPTEN)
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -s WASM=1 -s USE_WEBGL2=1 -s MIN_WEBGL_VERSION=2 -s MAX_WEBGL_VERSION=2 -Os \
        -s NO_FILESYSTEM=1 \
        --closure 1 \
        -s ELIMINATE\_DUPLICATE\_FUNCTIONS=1 \
        -fno-rtti \
        -fno-exceptions")
    set_target_properties(lightning_demo PROPERTIES LINK_FLAGS
		"-s DEMANGLE_SUPPORT=1 \
        --shell-file ${CMAKE_SOURCE_DIR}/src/minimal_shell.html"
    )
else()
    set_property(TARGET lightning_demo PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
endif()

target_link_libraries(lightning_demo sokol glm imgui)
add_dependencies(lightning_demo shaders_target)


source_group("" FILES ${SRCS})