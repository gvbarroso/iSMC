# CMake script for iSMC
# Authors:
#   Julien Dutheil
# Created: 16/04/2018

cmake_minimum_required (VERSION 3.10)
project (ismc CXX)

# Compile options
add_compile_options(-std=c++14 -Wall -Weffc++ -Wshadow -Wconversion -O3)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

IF(NOT NO_DEP_CHECK)
  SET(NO_DEP_CHECK FALSE CACHE BOOL
      "Disable dependencies check for building distribution only."
      FORCE)
ENDIF(NOT NO_DEP_CHECK)

IF(NO_DEP_CHECK)
  MESSAGE(STATUS "Dependencies checking disabled. Only distribution can be built.")
ELSE(NO_DEP_CHECK)


# Static linkage?
IF(NOT BUILD_STATIC)
  SET(BUILD_STATIC FALSE CACHE BOOL
      "Enable static linkage."
      FORCE)
ENDIF()
IF(BUILD_STATIC)
  MESSAGE(STATUS "Static linkage requested.")
  #SET(CMAKE_CXX_FLAGS "-static -static-libgcc -static-libstdc++ -static-glibc ${CMAKE_CXX_FLAGS}")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++")
ENDIF()


# Use SimpleHMM?
IF(NOT USE_LIBSIMPLEZIPHMM)
  set(USE_LIBSIMPLEZIPHMM TRUE CACHE BOOL
      "Use LibSimpleZipHMM."
      FORCE)
ENDIF()

# Find dependencies (add install directory to search)
if (CMAKE_INSTALL_PREFIX)
  set (CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
endif (CMAKE_INSTALL_PREFIX)

INCLUDE(GNUInstallDirs)
FIND_PACKAGE(bpp-phyl3 1.0.0 REQUIRED)

SET(THREADS_PREFER_PTHREAD_FLAG ON)
FIND_PACKAGE(Threads REQUIRED)

IF(USE_LIBSIMPLEZIPHMM)
  FIND_PACKAGE(SimpleZipHMM 1.0.0 REQUIRED)
  ADD_DEFINITIONS(-DSIMPLEZIPHMM)
ELSE()
  # As ZipHMM does not use CMake, we need to find the libraries manually:
  MESSAGE("TODO")
ENDIF()

# Boost:
SET(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE( Boost 1.36.0 COMPONENTS iostreams )
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
  SET(LIBS ${LIBS} ${Boost_LIBRARIES})
  MESSAGE("-- Boost libraries found here:")
  MESSAGE("   includes: ${Boost_INCLUDE_DIRS}")
  MESSAGE("   dynamic libraries: ${Boost_LIBRARIES}")
ENDIF()

#here is a useful function:
MACRO(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})

  IF(BUILD_STATIC)
    SET(${lib_name}_STATIC_NAMES lib${lib_name}.a)
    FIND_LIBRARY(${lib_name}_STATIC_LIBRARY NAMES ${${lib_name}_STATIC_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_STATIC_LIBRARY)
      MESSAGE(STATUS "Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
  ELSE()
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES} PATH_SUFFIXES lib${LIB_SUFFIX})
    IF(${lib_name}_LIBRARY)
      MESSAGE(STATUS "Library ${lib_name} found here:")
      MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
      MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    ELSE()
      MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF()
    #add the dependency:
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
  ENDIF()
  
ENDMACRO(IMPROVED_FIND_LIBRARY)

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

# BLAS and other dependencies:

IF(BUILD_STATIC)
  # Find the zlib installation
  IMPROVED_FIND_LIBRARY(LIBS z zlib.h)
  # Find the bz2 installation
  IMPROVED_FIND_LIBRARY(LIBS bz2 bzlib.h)
  # Find blas installation
  IMPROVED_FIND_LIBRARY(LIBS blas cblas.h)
ELSE(BUILD_STATIC)
  FIND_PACKAGE(ZLIB REQUIRED)
  IF(ZLIB_FOUND)
    INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${ZLIB_LIBRARIES})
    MESSAGE(STATUS "Zlib found here:")
    MESSAGE("   includes: ${ZLIB_INCLUDE_DIR}")
    MESSAGE("   libraries: ${ZLIB_LIBRARIES}")
  ENDIF()

  FIND_PACKAGE(BZip2 REQUIRED)
  IF(BZIP2_FOUND)
    INCLUDE_DIRECTORIES(${BZIP2_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${BZIP2_LIBRARIES})
    MESSAGE(STATUS "BZip2 found here:")
    MESSAGE("   includes: ${BZIP2_INCLUDE_DIR}")
    MESSAGE("   libraries: ${BZIP2_LIBRARIES}")
  ENDIF()

  #Does not work with static libraries:
  #SET(BLA_STATIC ${BUILD_STATIC})
  #SET(BLA_VENDOR "ATLAS")
  FIND_PACKAGE(BLAS REQUIRED)
ENDIF(BUILD_STATIC)

# Subdirectories
ADD_SUBDIRECTORY(src)

ENDIF(NO_DEP_CHECK)

# Packager
SET(CPACK_PACKAGE_NAME "ismc")
SET(CPACK_PACKAGE_VENDOR "Gustavo V. Barroso")
SET(CPACK_PACKAGE_VERSION "0.1.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "1")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The iSMC demography inferrence package.")
#Note (Julien): we will have to take care of that when we distribute the software.
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
#SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
#SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
SET(CPACK_SOURCE_IGNORE_FILES 
       "/build/" 
       "/\\\\.git/" 
       "/\\\\.gitignore" 
       ${CPACK_SOURCE_IGNORE_FILES}
       )

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

IF(NOT NO_DEP_CHECK)
IF (UNIX)

ENDIF(UNIX)

ENDIF(NOT NO_DEP_CHECK)

