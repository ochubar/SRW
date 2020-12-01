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

if(USE_OPENMP)
    add_definitions(-D_WITH_OMP)
    add_compile_options(-fopenmp -Wno-write-strings)
    set(USE_FFTW3 OFF)
    #LDFLAGS += -lfftw
else()
    add_definitions(-D_FFTW3)
    set(USE_FFTW3 ON)
    #LDFLAGS+= -lfftw3f -lfftw3
endif()

add_subdirectory(ext)
add_subdirectory(core)
add_subdirectory(lib)

if(BUILD_CLIENTS)
    add_subdirectory(clients)
endif()