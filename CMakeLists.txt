
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
	seal/seal.h
        descrs/fp2s.h
        descrs/spectrophores.h
        descrs/mnadescr.h
        machines/oneclasssvm.h
        machines/svr.h
        fields.h
        parser.h
)



set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

find_package(OpenBabel2 REQUIRED)
find_package(Boost REQUIRED)

FIND_PACKAGE( Flex  REQUIRED)
FIND_PACKAGE( Bison REQUIRED)

include_directories(${OPENBABEL2_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(${PROJECT_SOURCE_DIR}/seal)
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
#Integration CMAKE with flex/bison thanks to Andy Cedilnik Kitware Inc.

ADD_CUSTOM_TARGET(miparser echo "Creating grammar.tab.cpp")
#
# Create custom command for flex/lex (note the outputs)
#
ADD_CUSTOM_COMMAND(
   SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/cmdparams.lex"
   COMMAND ${FLEX_EXECUTABLE}
   ARGS -o"${CMAKE_CURRENT_BINARY_DIR}/lex.yy.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmdparams.lex"
   TARGET miparser
   OUTPUTS "${CMAKE_CURRENT_BINARY_DIR}/lex.yy.cpp"
   )

#
# Create custom commands for bison/yacc (note the DEPENDS)
#
ADD_CUSTOM_COMMAND(
  SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/grammar.y"
  COMMAND ${BISON_EXECUTABLE}
  ARGS -y -o"${CMAKE_CURRENT_BINARY_DIR}/grammar.tab.cpp"
             "${CMAKE_CURRENT_SOURCE_DIR}/grammar.y"
  TARGET miparser
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/lex.yy.cpp"
  OUTPUTS "${CMAKE_CURRENT_BINARY_DIR}/grammar.tab.cpp"
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set (SOURCE
        cmfa.cpp
        descfact.cpp
        kernel.cpp
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
        seal/seal.cpp
        descrs/fp2s.cpp
        descrs/spectrophores.cpp
        descrs/mnadescr.cpp
        machines/oneclasssvm.cpp
        machines/svr.cpp
        fields.cpp
        "${CMAKE_CURRENT_BINARY_DIR}/grammar.tab.cpp"
        )
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/grammar.tab.c GENERATED)

add_executable ( ${PROJECT} ${HEADERS} ${SOURCE} )
target_link_libraries(MCMF ${LIBS} nlopt)

install (TARGETS MCMF DESTINATION ${CMAKE_BINARY_DIR}/bin)

