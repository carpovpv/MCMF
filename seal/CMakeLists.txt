
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        seal)

set (HEADERS 
	seal.h
	ks_seal.h
        ../fields.h
)

set (SOURCE
	align.cpp
	ks_seal.cpp
	seal.cpp
        ../fields.cpp
        )

project (${PROJECT})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake/modules)
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_REENTERANT")

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${HEADERS} ${SOURCE})
target_link_libraries(seal ${LIBS} nlopt pthread)

install (TARGETS seal DESTINATION ${CMAKE_BINARY_DIR}/bin)

