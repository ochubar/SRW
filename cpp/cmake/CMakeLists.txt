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

    # Setup OpenMP as a linked library for all targets
    link_libraries(OpenMP::OpenMP_CXX)
    # Setup OpenMP include dir for all targets
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
else()
    set(SRW_DEFINITIONS -D_GNU_SOURCE -D__USE_XOPEN2K8 -DFFTW_ENABLE_FLOAT -D_GM_WITHOUT_BASE -DSRWLIB_STATIC -DNO_TIMER -DANSI_DECLARATORS -DTRILIBRARY -D_FFTW3)
endif()

# Setup definitions for all targets
add_definitions(${SRW_DEFINITIONS})

#ext
#auxparse
set(auxparse_source_files
    ../src/ext/auxparse/auxparse.cpp
    ../src/ext/auxparse/auxparse.h
)

add_library(auxparse OBJECT ${auxparse_source_files})
target_include_directories(auxparse PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src/ext/auxparse)
target_compile_definitions(auxparse PRIVATE ${SRW_DEFINITIONS})

#genesis
set(genesis_source_files
    ../src/ext/genesis/genesis_july08/all_com.c
    ../src/ext/genesis/genesis_july08/check.c
    ../src/ext/genesis/genesis_july08/diagno.c
    ../src/ext/genesis/genesis_july08/esource.c
    ../src/ext/genesis/genesis_july08/field.c
    ../src/ext/genesis/genesis_july08/incoherent.c
    ../src/ext/genesis/genesis_july08/initrun.c
    ../src/ext/genesis/genesis_july08/input.c
    ../src/ext/genesis/genesis_july08/loadbeam.c
    ../src/ext/genesis/genesis_july08/loadrad.c
    ../src/ext/genesis/genesis_july08/magfield.c
    ../src/ext/genesis/genesis_july08/main.c
    ../src/ext/genesis/genesis_july08/math.c
    ../src/ext/genesis/genesis_july08/mpi.c
    ../src/ext/genesis/genesis_july08/output.c
    ../src/ext/genesis/genesis_july08/partsim.c
    ../src/ext/genesis/genesis_july08/pushp.c
    ../src/ext/genesis/genesis_july08/rpos.c
    ../src/ext/genesis/genesis_july08/scan.c
    ../src/ext/genesis/genesis_july08/source.c
    ../src/ext/genesis/genesis_july08/stepz.c
    ../src/ext/genesis/genesis_july08/string.c
    ../src/ext/genesis/genesis_july08/tdepend.c
    ../src/ext/genesis/genesis_july08/timerec.c
    ../src/ext/genesis/genesis_july08/track.c
)

add_library(genesis OBJECT ${genesis_source_files})
target_link_libraries(genesis genmath)
target_include_directories(genesis PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src/ext/genesis/genesis_july08/)
target_compile_definitions(genesis PRIVATE ${SRW_DEFINITIONS})

#genmath
set(genmath_source_files
    ../src/ext/genmath/gmfft.cpp
    ../src/ext/genmath/gmfft.h
    ../src/ext/genmath/gmfit.cpp
    ../src/ext/genmath/gmfit.h
    ../src/ext/genmath/gminterp.cpp
    ../src/ext/genmath/gminterp.h
    ../src/ext/genmath/gmmeth.cpp
    ../src/ext/genmath/gmmeth.h
    ../src/ext/genmath/gmtrans.cpp
    ../src/ext/genmath/gmtrans.h
)

add_library(genmath OBJECT ${genmath_source_files})
target_link_libraries(genmath PUBLIC auxparse)
target_include_directories(genmath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src/ext/genmath/ ${CMAKE_CURRENT_SOURCE_DIR}/../src/lib)

if(USE_OPENMP)
    target_include_directories(genmath PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
endif()

target_compile_definitions(genmath PRIVATE ${SRW_DEFINITIONS})

#core
set(core_source_files
    ../src/core/srclcuti.cpp
    ../src/core/srcradint.cpp
    ../src/core/srctrjdt.cpp
    ../src/core/sremitpr.cpp
    ../src/core/srgsnbm.cpp
    ../src/core/srgtrjdt.cpp
    ../src/core/srisosrc.cpp
    ../src/core/srmagcnt.cpp
    ../src/core/srmagfld.cpp
    ../src/core/srmatsta.cpp
    ../src/core/sroptapt.cpp
    ../src/core/sroptcnt.cpp
    ../src/core/sroptdrf.cpp
    ../src/core/sroptel2.cpp
    ../src/core/sroptel3.cpp
    ../src/core/sroptelm.cpp
    ../src/core/sroptfoc.cpp
    ../src/core/sroptgrat.cpp
    ../src/core/sroptgtr.cpp
    ../src/core/sropthck.cpp
    ../src/core/sroptcryst.cpp
    ../src/core/sroptmat.cpp
    ../src/core/sroptpsh.cpp
    ../src/core/sroptshp.cpp
    ../src/core/sroptsmr.cpp
    ../src/core/sroptwgr.cpp
    ../src/core/sroptzp.cpp
    ../src/core/sroptzps.cpp
    ../src/core/srpersto.cpp
    ../src/core/srpowden.cpp
    ../src/core/srprdint.cpp
    ../src/core/srprgind.cpp
    ../src/core/srpropme.cpp
    ../src/core/srptrjdt.cpp
    ../src/core/srradinc.cpp
    ../src/core/srradint.cpp
    ../src/core/srradmnp.cpp
    ../src/core/srradstr.cpp
    ../src/core/srremflp.cpp
    ../src/core/srsase.cpp
    ../src/core/srsend.cpp
    ../src/core/srstowig.cpp
    ../src/core/srsysuti.cpp
    ../src/core/srthckbm.cpp
    ../src/core/srthckbm2.cpp
    ../src/core/srtrjaux.cpp
    ../src/core/srtrjdat.cpp
    ../src/core/srtrjdat3d.cpp
)

add_library(core OBJECT ${core_source_files})
target_link_libraries(core genesis)
target_compile_definitions(core PRIVATE ${SRW_DEFINITIONS})

target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../src/core)
target_include_directories(core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/lib>
)

#lib
set(srw_lib_source_files
    ../src/lib/srwlib.cpp
    ../src/lib/srwlib.h
#    ../src/lib/auxgpu.h
#    ../src/lib/auxgpu.cpp
    ../src/lib/srerror.cpp
    ../src/lib/srerror.h
)

add_library(srw STATIC
    ${srw_lib_source_files}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:genesis>
    $<TARGET_OBJECTS:genmath>
    $<TARGET_OBJECTS:auxparse>
)


# Add definitions needed to build SRW
target_compile_definitions(srw PUBLIC ${SRW_DEFINITIONS})

# Libraries in which SRW lib depends
target_link_libraries(srw core ${FFTW_DOUBLE_LIB} ${FFTW_FLOAT_LIB})

if (UNIX)
    # Math library is only needed by unix. When using Windows, VS does it for you.
    target_link_libraries(srw m)
endif()

target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/lib>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/core>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/ext/auxparse>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/ext/genesis>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/ext/genmath>
)

install(TARGETS srw LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# clients
# c
set(srw_clients_c_source_files
    ../src/clients/c/srwlclient.cpp
)

add_library(srwlclient STATIC ${srw_clients_c_source_files})

target_link_libraries(srwlclient srw)

target_include_directories(srwlclient PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/clients/c/>
)

install(TARGETS srwlclient RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# python
# for testing we will need the python interpreter
# we require python development headers
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)


set(srw_clients_python_source_files
    ../src/clients/python/srwlpy.cpp
)

add_library(srwlpy SHARED ${srw_clients_python_source_files})

target_include_directories(srwlpy PUBLIC ${Python_INCLUDE_DIRS})

target_include_directories(srwlpy PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../src/clients/python>)

# prevent cmake from creating a "lib" prefix
set_target_properties(srwlpy PROPERTIES PREFIX "")

if (UNIX)
    # Math library is only needed by unix. When using Windows, VS does it for you.
    target_link_libraries(srwlpy m ${FFTW_DOUBLE_LIB} ${FFTW_FLOAT_LIB} srw)

    set_target_properties(srwlpy PROPERTIES SUFFIX ".so")

    # The -shared flag is only valid for UNIX systems.
    set_target_properties(srwlpy PROPERTIES LINK_FLAGS "-shared")
endif()

if(APPLE)
    # We need the undefined dynamic_lookup to overcome issues with a statically linked
    # python. More details here: https://gitlab.kitware.com/vtk/vtk/-/issues/17214
    set_target_properties(srwlpy PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -shared")
endif()

if(WIN32)
  target_link_libraries(srwlpy srw ${Python_LIBRARIES})
  set_target_properties(srwlpy PROPERTIES SUFFIX ".pyd")
endif()

add_custom_command(TARGET srwlpy POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:srwlpy> ${CMAKE_SOURCE_DIR}/env/python/srwpy/
)

#install(TARGETS srwlpy LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
