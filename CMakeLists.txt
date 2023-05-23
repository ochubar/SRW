# Minimum cmake version
# Picked 3.12 as this one introduces a proper FindOpenMP
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Project name
project(SRW)

cmake_policy(SET CMP0074 NEW)

# Folder for helper modules (Find*.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cpp/cmake" ${CMAKE_MODULE_PATH})

# Uncomment for VERBOSE Makefiles. Useful for debugging.
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Useful debug for Find commands
#set(CMAKE_FIND_DEBUG_MODE ON)

# Uncomment for install messages
#set(CMAKE_INSTALL_MESSAGE ALWAYS)

# Uncomment for Debug Build
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

# Option for OpenMP build, OFF by default.
option(USE_OPENMP "Activate OpenMP build" OFF)

# Option for C client library to be built or not - See cpp/env/cmake/CMakeLists.txt
option(BUILD_CLIENT_C "Activate C Client library build" OFF)

# Option for Python client library to be built or not - See cpp/env/cmake/CMakeLists.txt
option(BUILD_CLIENT_PYTHON "Activate Python Client library build" OFF)

# Option for Igor Pro client library to be built or not - See cpp/env/cmake/CMakeLists.txt
option(BUILD_CLIENT_IGOR "Activate Igor Pro Client library build" OFF)


# Temporary path for FFTW in case we need to build it
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")


#if(UNIX)
#  # Use ./lib ./bin ./include ./share folder structure
#  include(GNUInstallDirs)
#else()
#  if (WIN32)
#    set(${CMAKE_INSTALL_LIBDIR} "lib")
#    set(${CMAKE_INSTALL_DATADIR} "share")
#    set(${CMAKE_INSTALL_INCLUDEDIR} "include")
#    set(${CMAKE_INSTALL_BINDIR} "bin")
#    message(STATUS "Setting installation destination on Windows to: ${CMAKE_INSTALL_PREFIX}")
#  else()
#    message(FATAL_ERROR "System not UNIX nor WIN32 - not implemented yet")
#  endif()
#endif()


# Force FFTW to look for STATIC libraries that we need
set(FFTW_USE_STATIC_LIBS ON)

# External Libraries
if(USE_OPENMP)
    message(STATUS "Looking for FFTW2")
    add_subdirectory(ext_lib/fftw)
else()
    message(STATUS "Looking for FFTW3")
    add_subdirectory(ext_lib/fftw3)
endif()

# Main library
add_subdirectory(cpp/cmake)
