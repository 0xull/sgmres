cmake_minimum_required(VERSION 3.20)
project(sgmres LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(SGMRES_NATIVE "Tune code generation for the build machine" ON)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(sgmres_options INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(sgmres_options INTERFACE -Wall -Wextra -Wpedantic -Wshadow)
  if(SGMRES_NATIVE)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-march=native HAVE_MARCH_NATIVE)
    if(HAVE_MARCH_NATIVE)
      target_compile_options(sgmres_options INTERFACE -march=native)
    endif()
  endif()
elseif(MSVC)
  target_compile_options(sgmres_options INTERFACE /W4)
endif()

add_library(sgmres src/problems.cpp)
target_include_directories(sgmres PUBLIC include)
target_link_libraries(sgmres PUBLIC Eigen3::Eigen PRIVATE sgmres_options)

add_executable(inspect apps/inspect.cpp)
target_link_libraries(inspect PRIVATE sgmres sgmres_options)

include(CTest)
if(BUILD_TESTING)
  add_executable(problems_test tests/problems_test.cpp)
  target_link_libraries(problems_test PRIVATE sgmres sgmres_options)
  add_test(NAME problems COMMAND problems_test)
endif()
