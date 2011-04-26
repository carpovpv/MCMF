
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        prognosis)

set (HEADERS 
	../Net.h
	../vector_io.h

)

set (SOURCE
	prognosis.cpp
	../NetLayer.cpp
        ../Net.cpp
        )

project (${PROJECT})

find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${HEADERS} ${SOURCE})
target_link_libraries(prognosis ${LIBS} )


install (TARGETS prognosis DESTINATION ${CMAKE_BINARY_DIR}/bin)


