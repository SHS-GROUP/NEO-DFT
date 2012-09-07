cmake_minimum_required (VERSION 2.8)
project (CChem)

SET(CMAKE_VERBOSE_MAKEFILE ON)

#enable_language(Fortran)

# if(CMAKE_SIZEOF_VOID_P EQUAL 8)
#     MESSAGE(STATUS "64 bit")
#     set(ARCH 64)
# else(CMAKE_SIZEOF_VOID_P EQUAL 8)
#     set(ARCH 32)
# endif(CMAKE_SIZEOF_VOID_P EQUAL 8)
# add_definitions(-DINTEGER${ARCH})


find_package(Threads)


include (${CMAKE_ROOT}/Modules/FindBLAS.cmake)
if(BLAS_FOUND)
    MESSAGE(STATUS "BLAS libraries: ${BLAS_LIBRARIES}")
    MESSAGE(STATUS "BLAS linker flags: ${BLAS_LINKER_FLAGS}")
    set(LIBRARIES "${LIBRARIES};${BLAS_LIBRARIES}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS_INIT} ${BLAS_LINKER_FLAGS}")
    set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS_INIT} ${BLAS_LINKER_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS_INIT} ${BLAS_LINKER_FLAGS}")
endif(BLAS_FOUND)


INCLUDE (CheckIncludeFiles)
# usage: CHECK_INCLUDE_FILES (<header> <RESULT_VARIABLE> )
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/config/config.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/config/config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/config/)

CONFIGURE_FILE(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/basis/config.hpp.in
  ${CMAKE_CURRENT_SOURCE_DIR}/src/basis/config.hpp)

# if(UNIX)
#     add_definitions(-DUNIX)
# endif(UNIX)

# string(TOUPPER ${CMAKE_SYSTEM_NAME} TARGET_)
# if(ARCH EQUAL 64)
#     set(TARGET_ "${TARGET_}64")
# endif()

# set(TARGET "${TARGET_}" CACHE STRING "TARGET")
# add_definitions(-D${TARGET})
# if("${TARGET}" MATCHES "64")
#     string(REGEX REPLACE "64" "" TARGET_ "${TARGET}")
#     add_definitions(-D${TARGET_})
# endif()
# MESSAGE(STATUS "TARGET: ${TARGET}")


# if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
#     set(LIBRARIES "${LIBRARIES};${CMAKE_THREAD_LIBS_INIT}")
# endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")


# find_library(CCHEM cchem)
# if (CCHEM)
#     GET_FILENAME_COMPONENT(CCHEM_PATH "${CCHEM}" PATH)
#     link_directories("${CCHEM_PATH}")
#     MESSAGE(STATUS "Found libcchem: ${CCHEM_PATH}")
#     set(LIBRARIES "${LIBRARIES};cchem_gamess;cchem;rysq")
#     if(UNIX)
#         set(LIBRARIES "${LIBRARIES};stdc++")
#     endif(UNIX)
#     set(Boost_ADDITIONAL_VERSIONS "1.46" "1.46.0")
#     #set(Boost_NO_SYSTEM_PATHS true)
#     find_package(Boost COMPONENTS thread)
#     #set(Boost_USE_MULTITHREADED OFF)
#     if(Boost_FOUND)
#         MESSAGE(STATUS "Found Boost: ${Boost_LIBRARY_DIRS}")
# 	MESSAGE(STATUS "Found Boost libraries: ${Boost_LIBRARIES}") 
# 	set(LIBRARIES "${LIBRARIES};${Boost_LIBRARIES}")
#     else()
#         MESSAGE(FATAL_ERROR "Boost Thread NOT FOUND")
#     endif()
# endif()


set(CUDA_ROOT "/usr/local/cuda" CACHE PATH "CUDA ROOT")
find_library(CUDART cudart "${CUDA_ROOT}/lib64")
if (CUDART)
    set(LIBRARIES "${LIBRARIES};${CUDART}")
endif()
find_library(CUBLAS cublas "${CUDA_ROOT}/lib64")
if (CUBLAS)
    set(LIBRARIES "${LIBRARIES};${CUBLAS}")
endif()


find_package(HDF5 COMPONENTS CXX)
if (HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_LIBRARIES}")
    set(LIBRARIES "${LIBRARIES};${HDF5_LIBRARIES}")
endif()


option(OPENMP "Enable OpenMP" ON)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
   set(CMAKE_Fortran_FLAGS "-std=legacy ${CMAKE_Fortran_FLAGS}")
   set(GFORTRAN)
   add_definitions(-DUNDERSCORE)
   if(ARCH EQUAL 64)
       set(CMAKE_Fortran_FLAGS "-fdefault-integer-8 ${CMAKE_Fortran_FLAGS}")
   endif(ARCH EQUAL 64)
   if(OPENMP)
       set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
   endif()
endif()


# sources
set(CCHEM_SOURCES)
#file(GLOB CCHEM_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp") 
# list(APPEND SOURCES "${PROJECT_SOURCE_DIR}/source/zunix.c")
# list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/source/global.F")
# list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/source/vector.F")
# list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/source/serial.F")
# list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/source/blas.F")

include_directories(src)
include_directories(src/externals)

add_subdirectory(src)

#list(APPEND SOURCES "${DFT_SOURCES}")

add_library(cchem ${CCHEM_SOURCES})

#add_library(cchem_gamess ${cchem_gamess_sources})

# # serial
# add_executable(gamess.serial.x ${SOURCES} "${PROJECT_SOURCE_DIR}/source/serial.f")


# # parallel
# link_directories("${PROJECT_SOURCE_DIR}/ddi/armci/lib")
# add_executable(gamess.ddi.x EXCLUDE_FROM_ALL ${SOURCES})
# add_subdirectory(ddi EXCLUDE_FROM_ALL)
# add_dependencies(gamess.ddi.x ddi)
# target_link_libraries(gamess.ddi.x ddi)
# target_link_libraries(gamess.ddi.x libarmci.a)

# find_package(MPI)
# if(MPI_FOUND)
#     include_directories(${MPI_INCLUDE_PATH})
#     set_target_properties(gamess.ddi.x PROPERTIES COMPILE_DEFINITIONS MPI)
#     target_link_libraries(gamess.ddi.x ${MPI_LIBRARIES})
# endif()


# # serial
# add_executable(gamess.vanilla.x ${SOURCES} "${PROJECT_SOURCE_DIR}/source/serial.f")
# set_target_properties(gamess.vanilla.x PROPERTIES COMPILE_FLAGS -DVANILLA)

# # common libraries
# target_link_libraries(gamess.serial.x ${LIBRARIES})
# target_link_libraries(gamess.ddi.x ${LIBRARIES})
# target_link_libraries(gamess.vanilla.x ${LIBRARIES})
