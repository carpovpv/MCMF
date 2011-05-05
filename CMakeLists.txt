
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
        kernels/linear.h
        kernels/sterick.h
        kernels/hydrophov.h
        kernels/abraham.h
	seal/ks_seal.h
	seal/seal.h
        descrs/fp2s.h
        descrs/spectrophores.h
        machines/oneclasssvm.h
        machines/svr.h
        fields.h
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
        kernels/linear.cpp
        kernels/sterick.cpp
        kernels/hydrophov.cpp
        kernels/abraham.cpp
	seal/ks_seal.cpp
	seal/seal.cpp
        descrs/fp2s.cpp
        descrs/spectrophores.cpp
        machines/oneclasssvm.cpp
        machines/svr.cpp
        fields.cpp
        )

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(OpenBabel2 REQUIRED)
find_package(Boost REQUIRED)

include_directories(${OPENBABEL2_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

find_package( Qt4 REQUIRED)
include (${QT_USE_FILE})

INCLUDE_DIRECTORIES(${QT_QTSQL_INCLUDE_DIR})

add_subdirectory(${PROJECT_SOURCE_DIR}/seal)
add_subdirectory(${PROJECT_SOURCE_DIR}/seal/thread)
add_subdirectory(${PROJECT_SOURCE_DIR}/nnet)
add_subdirectory(${PROJECT_SOURCE_DIR}/nnet/prognosis)
add_subdirectory(${PROJECT_SOURCE_DIR}/nnet/occprog)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/auc)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/dsc2set)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/ot)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/spc)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/spct)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/tnt)
add_subdirectory(${PROJECT_SOURCE_DIR}/utils/fp2)


set(LIBS ${OPENBABEL2_LIBRARIES})

add_executable ( ${PROJECT} ${HEADERS} ${SOURCE} )
target_link_libraries(MCMF ${LIBS} nlopt ${QT_LIBRARIES} ${QT_QTSQL_LIBRARY} )

install (TARGETS MCMF DESTINATION ${CMAKE_BINARY_DIR}/bin)

