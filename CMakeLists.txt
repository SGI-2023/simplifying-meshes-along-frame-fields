cmake_minimum_required(VERSION 3.16)
project(some-project)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_GLFW       "Use GLFW"           ON)

include(libigl)
igl_include(copyleft comiso)

include_directories(${CMAKE_SOURCE_DIR}/src)
# Add your project files
include_directories("include/")
if(USE_SOLUTION)
  file(GLOB SRCFILES solution/*.cpp)
else()
  file(GLOB SRCFILES src/*.cpp)
endif()

add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::glfw)

add_executable(main "main.cpp")
add_executable(view_mesh "view_mesh.cpp")
add_executable(cross_field "cross_field.cpp")
add_executable(frame_field_tutorial "frame_field_tutorial.cpp")
add_executable(alignment_check "alignment_check.cpp")

target_link_libraries(main core igl::core igl::opengl igl::glfw)
target_link_libraries(view_mesh core igl::core igl::opengl igl::glfw)
target_link_libraries(cross_field core igl::core igl::opengl igl::glfw)
target_link_libraries(frame_field_tutorial core igl::core igl::glfw igl_copyleft::comiso)
target_link_libraries(alignment_check core igl::core igl::opengl igl::glfw)
