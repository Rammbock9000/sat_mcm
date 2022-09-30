cmake_minimum_required (VERSION 2.8)
project (satscm)

# find ScaLP
find_library(SCALP_LIB
        NAMES libScaLP.so libScaLP.dylib
        DOC "Directory of the ScaLP library")
FIND_PATH(SCALP_H
        ScaLP/Solver.h
        DOC "Path to main ScaLP header file Solver.h"
        )
IF (SCALP_LIB AND SCALP_H)
    message(STATUS "Found ScaLP: ${SCALP_LIB} and ${SCALP_H}")
    add_definitions(-DUSE_SCALP)
ELSE (SCALP_LIB AND SCALP_H)
    message(STATUS "Could not find ScaLP (${SCALP_LIB} and ${SCALP_H}). Please specify path by including it in CMAKE_PREFIX_PATH (-DCMAKE_PREFIX_PATH=...)")
ENDIF (SCALP_LIB AND SCALP_H)

# find Z3
find_library(Z3_LIB
        NAMES libz3.so libz3.a libz3.dylib
        DOC "Directory of the Z3 library")
FIND_PATH(Z3_H
        z3.h
        DOC "Path to main Z3 header file z3.h"
        )
IF (Z3_LIB AND Z3_H)
    message(STATUS "Found Z3: ${Z3_LIB} and ${Z3_H}")
    add_definitions(-DUSE_Z3)
ELSE (Z3_LIB AND Z3_H)
    message(STATUS "Could not find Z3 (${Z3_LIB} and ${Z3_H}). Please specify path by including it in CMAKE_PREFIX_PATH (-DCMAKE_PREFIX_PATH=...)")
ENDIF (Z3_LIB AND Z3_H)

# find CaDiCaL
find_library(CADICAL_LIB
        NAMES libcadical.so libcadical.a libcadical.dylib
        DOC "Directory of the CaDiCaL library")
FIND_PATH(CADICAL_H
        cadical.hpp
        DOC "Path to main CaDiCaL header file cadical.hpp"
        )
IF (CADICAL_LIB AND CADICAL_H)
    message(STATUS "Found CaDiCaL: ${CADICAL_LIB} and ${CADICAL_H}")
    add_definitions(-DUSE_CADICAL)
ELSE (CADICAL_LIB AND CADICAL_H)
    message(STATUS "Could not find CaDiCaL (${CADICAL_LIB} and ${CADICAL_H}). Please specify path by including it in CMAKE_PREFIX_PATH (-DCMAKE_PREFIX_PATH=...)")
ENDIF (CADICAL_LIB AND CADICAL_H)

# set C++ standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# specify sources
file(GLOB SOURCES
        # main
        main.cpp

        # base class
        src/scm.cpp

        # derived classes
        src/scm_cadical.cpp
        src/scm_z3.cpp
     )

# add executable
add_executable(satscm ${SOURCES})

# include directories
include_directories(${CMAKE_SOURCE_DIR}/src)
if(SCALP_H)
    include_directories(${SCALP_H})
ENDIF(SCALP_H)
if(Z3_H)
    include_directories(${Z3_H})
ENDIF(Z3_H)
if(CADICAL_H)
    include_directories(${CADICAL_H})
ENDIF(CADICAL_H)

# link libraries
if(SCALP_LIB)
    target_link_libraries(satscm ${SCALP_LIB})
ENDIF(SCALP_LIB)
if(Z3_LIB)
    target_link_libraries(satscm ${Z3_LIB})
ENDIF(Z3_LIB)
if(CADICAL_LIB)
    target_link_libraries(satscm ${CADICAL_LIB})
ENDIF(CADICAL_LIB)

# install
install(FILES ${HEADERS}
        DESTINATION src
        )

install(TARGETS satscm
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )