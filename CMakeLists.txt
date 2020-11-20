# Minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# Project name
project(SRW)

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

add_subdirectory(ext_lib)

# OS Dependent configurations
if(UNIX)
    if (APPLE)  # MacOS
        add_definitions(-D__MAC__)
    else() # Linux
        add_definitions(-DLINUX)
    endif()
else() # Windows
    add_definitions(-DWIN32)
endif()

# Add definitions needed to build SRW
add_definitions(-D_GNU_SOURCE -D__USE_XOPEN2K8 -DFFTW_ENABLE_FLOAT)
add_definitions(-D_GM_WITHOUT_BASE -DSRWLIB_STATIC -DNO_TIMER)
add_definitions(-DANSI_DECLARATORS -DTRILIBRARY)

if(USE_OPENMP)
    add_definitions(-D_WITH_OMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -Wno-write-strings")
    set(USE_FFTW3 OFF)
    #LDFLAGS += -lfftw
else()
    add_definitions(-D_FFTW3)
    set(USE_FFTW3 ON)
    #LDFLAGS+= -lfftw3f -lfftw3
endif()

# Main source directory
add_subdirectory(src)