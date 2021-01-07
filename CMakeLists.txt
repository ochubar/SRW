# Minimum cmake version
# Picked 3.12 as this one introduces a proper FindOpenMP
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Project name
project(SRW)

cmake_policy(SET CMP0074 NEW)

include(GNUInstallDirs)

# Folder for helper modules (Find*.cmake)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Uncomment for VERBOSE Makefiles. Useful for debugging.
set(CMAKE_VERBOSE_MAKEFILE ON)

# Useful debug for Find commands
#set(CMAKE_FIND_DEBUG_MODE ON)

# Uncomment for install messages
#set(CMAKE_INSTALL_MESSAGE ALWAYS)

# Uncomment for Debug Build
#set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_BUILD_TYPE Release)

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

# Force FFTW to look for STATIC libraries that we need
set(FFTW_USE_STATIC_LIBS ON)

# External Libraries
add_subdirectory(ext_lib)

# Main library
add_subdirectory(src)



# TODO:
#- Finish the Python package setup -> Done
#- Test with OpenMP and make adjustments (FFTW2, etc) -> Done
#- Test that
#- Test with Windows
#- Make sure Oleg's way still works with the VisualC++.
