
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        dsc2set)

set (SOURCE
	dsc2set.cpp
        )

project (${PROJECT})

find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT}  ${SOURCE})

install (TARGETS dsc2set DESTINATION ${CMAKE_BINARY_DIR}/bin)


