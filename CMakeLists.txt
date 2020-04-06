cmake_minimum_required(VERSION 3.5)

project(QubitArraySimulator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_library(QuEST SHARED IMPORTED)
set_target_properties(QuEST PROPERTIES
	IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/../QuEST/libQuEST.so"
	INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/../QuEST/include"
)

set(SIMULATOR_SRC "main.cpp" "qubitarray.cpp" "randqalg.cpp")
add_executable(QubitArraySimulator ${SIMULATOR_SRC})
target_link_libraries(QubitArraySimulator QuEST)
