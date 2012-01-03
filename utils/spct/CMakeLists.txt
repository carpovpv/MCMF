
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        spct)

set (SOURCE
	spct.cpp
        )

project (${PROJECT})

find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${SOURCE})

install (TARGETS spct DESTINATION ${CMAKE_BINARY_DIR}/bin)


