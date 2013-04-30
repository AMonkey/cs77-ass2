cmake_minimum_required(VERSION 2.8.8)

option(USE_CLANG "build application with clang" ON)
option(USE_DEBUG "build with debug symbols"     ON)

# project
project(igl)

# flags output directories
if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(BIN_DIR_SUFFIX mk)
    set(MAKEFILES ON)
    set(CMAKE_VERBOSE_MAKEFILE OFF)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    if(USE_CLANG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
    endif(USE_CLANG)
    if(USE_DEBUG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    endif(USE_DEBUG)
    # may need to wrap / modify next two lines to depend on USE_CLANG
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGL_GLEXT_PROTOTYPES")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -lGL -lGLU -lglut")
elseif(CMAKE_GENERATOR STREQUAL "Xcode")
    set(BIN_DIR_SUFFIX xcode)
    set(XCODE ON)
    set_property(TARGET PROPERTY XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
elseif(CMAKE_GENERATOR STREQUAL "Visual Studio 10 Win64")
    set(BIN_DIR_SUFFIX vs10)
    add_definitions(-D_SCL_SECURE_NO_WARNINGS)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(CMAKE_GENERATOR STREQUAL "Unix Makefiles")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/${BIN_DIR_SUFFIX})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/${BIN_DIR_SUFFIX})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/${BIN_DIR_SUFFIX})

#include directories
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/src/ext)

#opengl and glut
find_package(GLUT)
find_package(OpenGL)

add_library(igl
    src/common/common.h
    src/common/debug.h
    src/common/debug.cpp
    src/common/iterators.h
    src/common/json.h
    src/common/json.cpp
    src/common/std.h
    src/common/std_utils.h
    src/common/stream.h
    
    src/vmath/frame.h
    src/vmath/geom.cpp
    src/vmath/geom.h
    src/vmath/interpolate.cpp
    src/vmath/interpolate.h
    src/vmath/mat.h
    src/vmath/random.h
    src/vmath/range.h
    src/vmath/ray.h
    src/vmath/stdmath.h
    src/vmath/transform.h
    src/vmath/vec.h
    src/vmath/vmath.h
    
    src/igl/camera.h
    src/igl/camera.cpp
    src/igl/draw.cpp
    src/igl/draw.h
    src/igl/gizmo.cpp
    src/igl/gizmo.h
    src/igl/gl_utils.cpp
    src/igl/gl_utils.h
    src/igl/gls.h
    src/igl/igl.h
    src/igl/image.h
    src/igl/image.cpp
    src/igl/light.h
    src/igl/light.cpp
    src/igl/material.h
    src/igl/material.cpp
    src/igl/node.h
    src/igl/node.cpp
    src/igl/primitive.h
    src/igl/primitive.cpp
    src/igl/scene.h
    src/igl/scene.cpp
    src/igl/serialize.h
    src/igl/serialize.cpp
    src/igl/shade.h
    src/igl/shade.cpp
    src/igl/shape.h
    src/igl/shape.cpp
    src/igl/tesselate.cpp
    src/igl/tesselate.h
    
    src/ext/lodepng/lodepng.h
    src/ext/lodepng/lodepng.cpp
)

target_link_libraries(igl ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${GLUT_LIBRARIES})

add_executable(view src/apps/view.cpp)
target_link_libraries(view igl)

if(CMAKE_GENERATOR STREQUAL "Xcode")
    set_property(TARGET igl view PROPERTY XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD c++0x)
    set_property(TARGET igl view PROPERTY XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY libc++)
endif(CMAKE_GENERATOR STREQUAL "Xcode")
