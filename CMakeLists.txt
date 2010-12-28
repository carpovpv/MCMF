
cmake_minimum_required (VERSION 2.6)
SET(CMAKE_BUILD_TYPE "Release")
set (PROJECT 
        MCMF)

project (${PROJECT})

set (HEADERS 
	kernel.h
	cmfa.h
	svm.h
	descfact.h
        machine.h
	kernels/electro.h
	kernels/gauss.h
	kernels/hydropho.h
	kernels/steric.h
	kernels/tanimoto.h
	seal/ks_seal.h
	seal/seal.h
        descrs/fp2.h
        machines/oneclasssvm.h
        machines/svr.h
)

set (SOURCE
	cmfa.cpp
	descfact.cpp
	mcmf.cpp
	svm.cpp
        machine.cpp
	kernels/electro.cpp
	kernels/gauss.cpp
	kernels/hydropho.cpp
	kernels/steric.cpp
	kernels/tanimoto.cpp
	seal/ks_seal.cpp
	seal/seal.cpp
	descrs/fp2.cpp
        machines/oneclasssvm.cpp
        machines/svr.cpp
        )

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(OpenBabel2 REQUIRED)
find_package(Boost REQUIRED)

include_directories(${OPENBABEL2_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(${PROJECT_SOURCE_DIR}/seal)

set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${HEADERS} ${SOURCE} )
target_link_libraries(MCMF ${LIBS} nlopt )

install (TARGETS MCMF DESTINATION ${CMAKE_BINARY_DIR}/bin)

