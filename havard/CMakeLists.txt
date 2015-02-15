project(proj1)
cmake_minimum_required(VERSION 2.8)

enable_language(CXX)

# MPI stuff
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
list(APPEND CMAKE_CXX_FLAGS ${MPI_C_FLAGS})

# openMP stuff
find_package(OpenMP REQUIRED)
list(APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})

# -Wpendantic gives warnings form the mpi library, so not added
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

add_executable(serialAndOpenMP serialAndOpenMP.cpp)
#target_link_libraries(serial ${MPI_CXX_LIBRARIES})

add_executable(mpi mpi.cpp)
target_link_libraries(mpi ${MPI_CXX_LIBRARIES})


# Use example in examples/common/ folder !!!!!!!!!!!!!

#---------------------------------------------------------------------------
# Some regular mistakes in CMake :
# http://voices.canonical.com/jussi.pakkanen/2013/03/26/a-list-of-common-cmake-antipatterns/

# Some useful variables
# http://www.cmake.org/Wiki/CMake_Useful_Variables

# How options() work:
# http://stackoverflow.com/questions/5998186/cmake-adding-command-line-options
