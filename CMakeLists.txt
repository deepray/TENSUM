cmake_minimum_required(VERSION 2.8.9)
project(tensum_solver)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall -O3")

if(CMAKE_HOST_SYSTEM_NAME STREQUAL Linux)
    set(LDFLAGS "-ldl")
else()
    set(LDFLAGS "")
endif()
#set(CMAKE_BUILD_TYPE Release)

# Find MPI
find_package(MPI REQUIRED)

# Mesh partitioner
file(GLOB SOURCES_MESH "grid_gen/*.cc"
                       "grid_gen/*.h"
                       )
add_executable(grid_part ${SOURCES_MESH})

# Solver
file(GLOB SOURCES_SOLVER "src/*.h"
                         "src/*.cc"
                         "extern/fparser/fparser.cc"
        )
add_executable(tensum ${SOURCES_SOLVER})
target_include_directories(tensum PRIVATE ${MPI_CXX_INCLUDE_PATH} extern/fparser)
target_link_libraries(tensum ${LDFLAGS} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})
