cmake_minimum_required(VERSION 3.10)
project(GeneVariants)

set(CMAKE_CXX_STANDARD 17)

# include spoa locally
add_subdirectory(external/spoa)

include_directories(external/spoa/include)

add_executable(project project.cpp)

target_link_libraries(project spoa)
