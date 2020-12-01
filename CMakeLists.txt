# Minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Project name
project(SRW)

cmake_policy(SET CMP0074 NEW)

include(GNUInstallDirs)

# Folder for helper modules (Find*.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Uncomment for VERBOSE Makefiles. Useful for debugging.
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_INSTALL_MESSAGE ALWAYS)

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
add_subdirectory(ext_lib/fftw3)

# Main library
add_subdirectory(src)