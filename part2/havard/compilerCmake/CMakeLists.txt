project(name)
cmake_minimum_required(VERSION 2.8)

enable_language(CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -pedantic")

#---------------------------------------------------------------------
# test stuff!
option(ENABLE_TE "Enable test?" ON)
#set(ENABLE_TE OFF)
if(ENABLE_TE)
    list(APPEND COMMON_DEFINITIONS -DHAVE_TE)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    list(APPEND COMMON_DEFINITIONS -DCOMP_GNU)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    list(APPEND COMMON_DEFINITIONS -DCOMP_INTEL)
endif()

add_definitions(${COMMON_DEFINITIONS})
#---------------------------------------------------------------------
add_executable(ex6 ex6.cpp)
#target_link_libraries(name ${MPI_CXX_LIBRARIES})



