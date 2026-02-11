# CMake script for iSMC
# Authors:
#   Julien Dutheil
# Created: 16/04/2018

cmake_minimum_required(VERSION 3.10)
project(ismc CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compile options
add_compile_options(
    -Wall
    -Weffc++
    -Wshadow
    -Wconversion
    -O3
)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    # May help on MacOS:
    add_compile_options(
        -DACCELERATE_NEW_LAPACK
        -DACCELERATE_LAPACK_ILP64
        -flax-vector-conversions
    )
    add_link_options(-framework Accelerate)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        RelWithDebInfo
        CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE
    )
endif(NOT CMAKE_BUILD_TYPE)

if(NOT NO_DEP_CHECK)
    set(NO_DEP_CHECK
        FALSE
        CACHE BOOL
        "Disable dependencies check for building distribution only."
        FORCE
    )
endif(NOT NO_DEP_CHECK)

if(NO_DEP_CHECK)
    message(
        STATUS
        "Dependencies checking disabled. Only distribution can be built."
    )
else(NO_DEP_CHECK)
    # Static linkage?
    if(NOT BUILD_STATIC)
        set(BUILD_STATIC FALSE CACHE BOOL "Enable static linkage." FORCE)
    endif()
    if(BUILD_STATIC)
        message(STATUS "Static linkage requested.")
	#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    endif()

    # Use SimpleHMM?
    if(NOT USE_LIBSIMPLEZIPHMM)
        set(USE_LIBSIMPLEZIPHMM TRUE CACHE BOOL "Use LibSimpleZipHMM." FORCE)
    endif()

    # Find dependencies (add install directory to search)
    if(CMAKE_INSTALL_PREFIX)
        set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
    endif(CMAKE_INSTALL_PREFIX)

    include(GNUInstallDirs)
    find_package(bpp-phyl 13 REQUIRED)

    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)

    if(USE_LIBSIMPLEZIPHMM)
        find_package(SimpleZipHMM 1.0.0 REQUIRED)
        add_definitions(-DSIMPLEZIPHMM)
    else()
        # As ZipHMM does not use CMake, we need to find the libraries manually:
        message("TODO")
    endif()

    # Boost:
    set(Boost_USE_STATIC_LIBS ${BUILD_STATIC})
    set(Boost_USE_MULTITHREADED ON)
    find_package(Boost CONFIG COMPONENTS iostreams)
    if(Boost_FOUND)
        include_directories(${Boost_INCLUDE_DIRS})
        set(LIBS ${LIBS} ${Boost_LIBRARIES})
        message("-- Boost libraries found here:")
        message("   includes: ${Boost_INCLUDE_DIRS}")
        message("   libraries: ${Boost_LIBRARIES}")
    endif()

    #here is a useful function:
    macro(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
        #start:
        find_path(${lib_name}_INCLUDE_DIR ${include_to_find})
        include_directories(${${lib_name}_INCLUDE_DIR})

        if(BUILD_STATIC)
            set(${lib_name}_STATIC_NAMES lib${lib_name}.a)
            find_library(
                ${lib_name}_STATIC_LIBRARY
                NAMES ${${lib_name}_STATIC_NAMES}
                PATH_SUFFIXES lib${LIB_SUFFIX}
            )
            if(${lib_name}_STATIC_LIBRARY)
                message(STATUS "Library ${lib_name} found here:")
                message("   includes: ${${lib_name}_INCLUDE_DIR}")
                message("   static libraries: ${${lib_name}_STATIC_LIBRARY}")
            else()
                message(FATAL_ERROR "${lib_name} required but not found.")
            endif()
            #add the dependency:
            set(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_STATIC_LIBRARY})
        else()
            set(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
            find_library(
                ${lib_name}_LIBRARY
                NAMES ${${lib_name}_NAMES}
                PATH_SUFFIXES lib${LIB_SUFFIX}
            )
            if(${lib_name}_LIBRARY)
                message(STATUS "Library ${lib_name} found here:")
                message("   includes: ${${lib_name}_INCLUDE_DIR}")
                message("   dynamic libraries: ${${lib_name}_LIBRARY}")
            else()
                message(FATAL_ERROR "${lib_name} required but not found.")
            endif()
            #add the dependency:
            set(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
        endif()
    endmacro(IMPROVED_FIND_LIBRARY)

    # Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
    # standard install location
    if(CMAKE_INSTALL_PREFIX)
        set(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
    endif(CMAKE_INSTALL_PREFIX)

    # BLAS and other dependencies:

    if(BUILD_STATIC)
        # Find the zlib installation
        improved_find_library(LIBS z zlib.h)
        # Find the bz2 installation
        improved_find_library(LIBS bz2 bzlib.h)
        # Find blas installation
        improved_find_library(LIBS blas cblas.h)
    else(BUILD_STATIC)
        find_package(ZLIB REQUIRED)
        if(ZLIB_FOUND)
            include_directories(${ZLIB_INCLUDE_DIR})
            set(LIBS ${LIBS} ${ZLIB_LIBRARIES})
            message(STATUS "Zlib found here:")
            message("   includes: ${ZLIB_INCLUDE_DIR}")
            message("   libraries: ${ZLIB_LIBRARIES}")
        endif()

        find_package(BZip2 REQUIRED)
        if(BZIP2_FOUND)
            include_directories(${BZIP2_INCLUDE_DIR})
            set(LIBS ${LIBS} ${BZIP2_LIBRARIES})
            message(STATUS "BZip2 found here:")
            message("   includes: ${BZIP2_INCLUDE_DIR}")
            message("   libraries: ${BZIP2_LIBRARIES}")
        endif()

        #Does not work with static libraries:
        #SET(BLA_STATIC ${BUILD_STATIC})
        #SET(BLA_VENDOR "ATLAS")
        find_package(BLAS REQUIRED)
    endif(BUILD_STATIC)

    # Subdirectories
    add_subdirectory(src)
endif(NO_DEP_CHECK)

# Packager
set(CPACK_PACKAGE_NAME "ismc")
set(CPACK_PACKAGE_VENDOR "Gustavo V. Barroso")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The iSMC demography inferrence package.")
#Note (Julien): we will have to take care of that when we distribute the software.
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
#SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
#SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
set(CPACK_SOURCE_GENERATOR "TGZ")
# /!\ This assumes that an external build is used
set(CPACK_SOURCE_IGNORE_FILES
    "/build/"
    "/\\\\.git/"
    "/\\\\.gitignore"
    ${CPACK_SOURCE_IGNORE_FILES}
)

set(CPACK_SOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}"
)
set(CPACK_DEBSOURCE_PACKAGE_FILE_NAME
    "${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig"
)
include(CPack)

#This adds the 'dist' target
add_custom_target(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)

if(NOT NO_DEP_CHECK)
    if(UNIX)
    endif(UNIX)
endif(NOT NO_DEP_CHECK)
