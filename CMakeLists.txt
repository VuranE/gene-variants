cmake_minimum_required(VERSION 3.10)
project(GeneVariants)

set(CMAKE_CXX_STANDARD 17)

# Add SPOA as a subdirectory (build it as part of your project)
add_subdirectory(external/spoa)

# Include SPOA headers
include_directories(external/spoa/include)

# Add your executable
add_executable(project project.cpp)

# Link your executable with SPOA library
target_link_libraries(project spoa)
