
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        occ)

set (HEADERS 
	Net.h
	vector_io.h
	types.h
)

set (SOURCE
	occ.cpp
	NetLayer.cpp
        Net.cpp
	occ.cpp
        )

project (${PROJECT})

find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${HEADERS} ${SOURCE})

install (TARGETS occ DESTINATION ${CMAKE_BINARY_DIR}/bin)


