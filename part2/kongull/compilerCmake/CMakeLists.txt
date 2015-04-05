project(poisson-mpi)

cmake_minimum_required(VERSION 2.8)

#add_subdirectory(../common common)

add_definitions(${COMMON_DEFINITIONS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_C_FLAGS}")
include_directories(${COMMON_INCLUDE_PATHS})


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

# Check which compiler is used
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
		# using G
		list(APPEND COMP_DEFINITIONS -DCOMP_GNU)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
		# using Intel C
		list(APPEND COMP_DEFINITIONS -DCOMP_INTEL)
endif()

add_definitions(${COMP_DEFINITIONS})

# Send the include directories to the compiler
include_directories(${INCLUDES})

# Add a library with common sources
add_library(newcommon fst.f)

# The C based poisson application
add_executable(poisson-mpi poisson-mpi.c)

# Link to libraries
target_link_libraries(poisson-mpi newcommon ${DEPLIBS_C} ${COMMON_LIBRARIES})

# We want to link this as a C program
set_target_properties(poisson-mpi PROPERTIES LINKER_LANGUAGE C)

