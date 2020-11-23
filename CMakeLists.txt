# Minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Project name
project(SRW LANGUAGES C)

cmake_policy(SET CMP0074 NEW)

include(GNUInstallDirs)

# Folder for helper modules (Find*.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Uncomment for VERBOSE Makefiles. Useful for debugging.
set(CMAKE_VERBOSE_MAKEFILE ON)

# Option for OpenMP build, OFF by default.
option(USE_OPENMP "Activate OpenMP build" OFF)

# Option for Client libraries to be built or not - See src/clients/CMakeLists.txt
option(BUILD_CLIENTS "Activate Client libraries build" OFF)

# Option for C client library to be built or not - See src/clients/CMakeLists.txt
option(BUILD_CLIENT_C "Activate C Client library build" OFF)

# Option for Python client library to be built or not - See src/clients/CMakeLists.txt
option(BUILD_CLIENT_PYTHON "Activate Python Client library build" OFF)

# Option for Igor Pro client library to be built or not - See src/clients/CMakeLists.txt
option(BUILD_CLIENT_IGOR "Activate Igor Pro Client library build" OFF)

# Temporary path for FFTW in case we need to build it
set(STAGED_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/stage)
message(STATUS "${PROJECT_NAME} staged install: ${STAGED_INSTALL_PREFIX}")

# External Libraries
add_subdirectory(ext_lib)

if(USE_OPENMP)
    set(CORE_DEPENDS fftw_external)
    #LDFLAGS += -lfftw
else()
    set(CORE_DEPENDS fftw3_external)
    #LDFLAGS+= -lfftw3f -lfftw3
endif()

include(ExternalProject)
ExternalProject_Add(${PROJECT_NAME}_core
  DEPENDS
    ${CORE_DEPENDS}
  SOURCE_DIR
    ${CMAKE_CURRENT_LIST_DIR}/src
  CMAKE_ARGS
    -DFFTW_ROOT=${FFTW_ROOT}
    -DBUILD_CLIENTS=${BUILD_CLIENTS}
    -DBUILD_CLIENT_PYTHON=${BUILD_CLIENT_PYTHON}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  CMAKE_CACHE_ARGS
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
  BUILD_ALWAYS
    1
  INSTALL_COMMAND
    ""
  )
