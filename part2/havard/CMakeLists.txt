project(poisson)

cmake_minimum_required(VERSION 2.8)

# We use the C language
enable_language(C)
enable_language(CXX)
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
  set(DEPLIBS_CXX ${DEPLIBS_CXX} ${MPI_CXX_LIBRARIES})
  set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
  list(APPEND CMAKE_CXX_FLAGS ${MPI_C_FLAGS})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    list(APPEND CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
  endif()
endif()

# Add some extra compilerflags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")

# Check which compiler is used
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    list(APPEND COMP_DEFINITIONS -DCOMP_GNU)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    list(APPEND COMP_DEFINITIONS -DCOMP_INTEL)
endif()

add_definitions(${COMP_DEFINITIONS})

# Send the include directories to the compiler
include_directories(${INCLUDES})

# Add a library with common sources
add_library(common fst.f)

# Temporary test program
add_executable(temp temp.cpp)
target_link_libraries(temp ${DEPLIBS_CXX})

# The C++ based poisson application
add_executable(poisson poisson.cpp)

# Link to libraries
target_link_libraries(poisson common ${DEPLIBS_CXX})

# We want to link this as a C++ program
set_target_properties(poisson PROPERTIES LINKER_LANGUAGE CXX)
