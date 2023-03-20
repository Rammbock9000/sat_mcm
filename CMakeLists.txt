cmake_minimum_required (VERSION 2.8)
project (satscm)

# setup for correctly linking to shared libraries
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# find Z3
find_library(Z3_LIB
        NAMES libz3.a libz3.so libz3.dylib
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

# find pthread lib for Glucose Syrup
find_library(PTHREAD_LIB
        NAMES libpthread.so libpthread.a libpthread.dylib
        DOC "pthread lib used by Glucose Syrup"
)

# find Glucose Syrup
find_library(SYRUP_LIB
        NAMES libsyrup.so libsyrup.a libsyrup.dylib
        DOC "Directory of the Glucose Syrup library")
FIND_PATH(SYRUP_H
        NAMES core/Solver.h
        DOC "Path to main Glucose Syrup header file Solver.h"
        )
IF (SYRUP_LIB AND SYRUP_H AND PTHREAD_LIB)
    message(STATUS "Found Glucose Syrup: ${SYRUP_LIB} and ${SYRUP_H} and ${PTHREAD_LIB}")
    add_definitions(-DUSE_SYRUP)
ELSE (SYRUP_LIB AND SYRUP_H AND PTHREAD_LIB)
    message(STATUS "Could not find Glucose Syrup (${SYRUP_LIB} and ${SYRUP_H} and ${PTHREAD_LIB}). Please specify path by including it in CMAKE_PREFIX_PATH (-DCMAKE_PREFIX_PATH=...)")
ENDIF (SYRUP_LIB AND SYRUP_H AND PTHREAD_LIB)

# set C++ standard
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# specify sources
file(GLOB lib_sources
        # base class
        src/mcm.cpp

        # derived classes
        src/mcm_cadical.cpp
        src/mcm_z3.cpp
        src/mcm_syrup.cpp
        )

file(GLOB header_sources
        # base class
        src/mcm.h

        # derived classes
        src/mcm_cadical.h
        src/mcm_z3.h
        src/mcm_syrup.h
        )

# compilation settings
add_compile_options(-O3)

# create library
add_library(satmcm-lib ${lib_sources})
set_target_properties(satmcm-lib PROPERTIES OUTPUT_NAME "satmcm")

# include directories
target_include_directories(satmcm-lib PUBLIC ${CMAKE_SOURCE_DIR}/src)
if(Z3_H)
    target_include_directories(satmcm-lib PUBLIC ${Z3_H})
ENDIF(Z3_H)
if(CADICAL_H)
    target_include_directories(satmcm-lib PUBLIC ${CADICAL_H})
ENDIF(CADICAL_H)
if(SYRUP_H)
    target_include_directories(satmcm-lib PUBLIC ${SYRUP_H})
ENDIF(SYRUP_H)

# link libraries
if(Z3_LIB)
    target_link_libraries(satmcm-lib PUBLIC ${Z3_LIB})
ENDIF(Z3_LIB)
if(CADICAL_LIB)
    target_link_libraries(satmcm-lib PUBLIC ${CADICAL_LIB})
ENDIF(CADICAL_LIB)
if(PTHREAD_LIB)
    target_link_libraries(satmcm-lib PUBLIC ${PTHREAD_LIB})
ENDIF(PTHREAD_LIB)
if(SYRUP_LIB)
    target_link_libraries(satmcm-lib PUBLIC ${SYRUP_LIB})
ENDIF(SYRUP_LIB)

# create executable, set name and link it to library
add_executable(satmcm-bin main.cpp)
set_target_properties(satmcm-bin PROPERTIES OUTPUT_NAME "satmcm")
target_link_libraries(satmcm-bin PUBLIC satmcm-lib)

# define install paths
set(bin_dest "bin/")
set(include_dest "include/")
set(lib_dest "lib/")

# install binary
install(TARGETS satmcm-bin DESTINATION "${bin_dest}")

# install lib
install(TARGETS satmcm-lib DESTINATION "${lib_dest}")

# install header files
install(FILES ${header_sources} DESTINATION "${include_dest}")