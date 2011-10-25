
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        pseal)

set (SOURCE
	pseal.cpp
        )

project (${PROJECT})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../cmake/modules)
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -D_REENTERANT")

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${SOURCE})
target_link_libraries(pseal ${LIBS} nlopt pthread)

install (TARGETS pseal DESTINATION ${CMAKE_BINARY_DIR}/bin)

