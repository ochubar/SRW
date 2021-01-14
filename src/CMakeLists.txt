# OS Dependent configurations
if(UNIX)
    if (APPLE)  # MacOS
        add_definitions(-D__MAC__)
    else() # Linux
        add_definitions(-DLINUX)
    endif()
    add_compile_options(-fPIC)
else() # Windows
    set(WINDOWS_EXPORT_ALL_SYMBOLS ON)
    add_definitions(-DWIN32 -D_WINDOWS -D_USRDLL -D_CRT_SECURE_NO_WARNINGS -DNON_UNIX_STDIO -D__VC__)
endif()


message(STATUS "Compiler is: ${CMAKE_C_COMPILER} | ${CMAKE_CXX_COMPILER}")

if(USE_OPENMP)
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP libomp location: ${OpenMP_libomp_LIBRARY} | Include: ${OpenMP_CXX_INCLUDE_DIRS}")
    message(STATUS "OpenMP Flags: C -> ${OpenMP_C_FLAGS} | CXX -> ${OpenMP_CXX_FLAGS}")
    set(SRW_DEFINITIONS -D_GNU_SOURCE -D__USE_XOPEN2K8 -DFFTW_ENABLE_FLOAT -D_GM_WITHOUT_BASE -DSRWLIB_STATIC -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -D_WITH_OMP)
else()
    set(SRW_DEFINITIONS -D_GNU_SOURCE -D__USE_XOPEN2K8 -DFFTW_ENABLE_FLOAT -D_GM_WITHOUT_BASE -DSRWLIB_STATIC -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -D_FFTW3)
endif()

add_subdirectory(ext)
add_subdirectory(core)
add_subdirectory(lib)

add_subdirectory(clients)