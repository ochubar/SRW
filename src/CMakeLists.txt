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

add_compile_options(-fPIC)

message(STATUS "Compiler is: ${CMAKE_C_COMPILER} | ${CMAKE_CXX_COMPILER}")

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    add_definitions(-D_WITH_OMP)
    message(STATUS "OpenMP libomp location: ${OpenMP_libomp_LIBRARY} | Include: ${OpenMP_CXX_INCLUDE_DIRS}")
    message(STATUS "OpenMP Flags: C -> ${OpenMP_C_FLAGS} | CXX -> ${OpenMP_CXX_FLAGS}")
else()
    add_definitions(-D_FFTW3)
endif()

add_subdirectory(ext)
add_subdirectory(core)
add_subdirectory(lib)

if(BUILD_CLIENTS)
    add_subdirectory(clients)
endif()