cmake_minimum_required(VERSION 3.16)
project(some-project)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_GLFW       "Use GLFW"           ON)

include(libigl)

include_directories(${CMAKE_SOURCE_DIR}/src)
# Add your project files
include_directories("include/")
if(USE_SOLUTION)
  file(GLOB SRCFILES solution/*.cpp)
else()
  file(GLOB SRCFILES src/*.cpp include/*.h)
endif()

add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::glfw)

add_executable(main "main.cpp")

target_link_libraries(main core igl::core igl::opengl igl::glfw)
