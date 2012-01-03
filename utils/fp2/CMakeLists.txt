
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        fp2)

set (SOURCE
	fp2.cpp
        )

project (${PROJECT})

find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${SOURCE})
target_link_libraries(fp2 ${LIBS}  )

install (TARGETS fp2 DESTINATION ${CMAKE_BINARY_DIR}/bin)


