project(poisson-mpi)

cmake_minimum_required(VERSION 2.6)

# We use the C language
enable_language(C)
enable_language(Fortran)

# CMake Options to enable/disable MPI and OpenMP
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
endif()

if(MPI_FOUND)
  add_definitions(-DHAVE_MPI=1)
  set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  set(DEPLIBS_C ${DEPLIBS_C} ${MPI_C_LIBRARIES})
  set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()
endif()

# Send the include directories to the compiler
include_directories(${INCLUDES})

# Add a library with common sources
add_library(common fst.f)

# The C based poisson application
add_executable(poisson-mpi poisson-mpi.c)

# Link to libraries
target_link_libraries(poisson-mpi common ${DEPLIBS_C})

# We want to link this as a C program
set_target_properties(poisson-mpi PROPERTIES LINKER_LANGUAGE C)

 ########################################### 
#project(poisson-mpi) 

#cmake_minimum_required(VERSION 2.8)

#enable_language(CXX)
#enable_language(Fortran)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pedantic")

## MPI - stuff
#find_package(MPI REQUIRED)
#include_directories(${MPI_CXX_INCLUDE_PATH} ${MPI_Fortran_INCLUDE_PATH})

## OpenMP - stuff
#find_package(OpenMP)
#add_definitions(-DHAVE_OPENMP=1)
#set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

## Fortran stuff
## Add a library with common sources
## Add a library with common sources
#add_library(common fst.f)
#add_library(common fst.f)

#add_definitions(${COMMON_DEFINITIONS})
##---------------------------------------------------------------------
#add_executable(poisson-mpi poisson-mpi.cpp)
#target_link_libraries(poisson-mpi common ${MPI_CXX_LIBRARIES})



