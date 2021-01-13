# for testing we will need the python interpreter
find_package(PythonInterp REQUIRED)

# we require python development headers
find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} EXACT REQUIRED)


set(srw_clients_python_source_files
    srwlpy.cpp
)

add_library(srwlpy SHARED ${srw_clients_python_source_files})

target_include_directories(srwlpy PUBLIC ${PYTHON_INCLUDE_DIRS})

target_include_directories(srwlpy PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>)

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
  target_link_libraries(srwlpy srw ${PYTHON_LIBRARIES})
  set_target_properties(srwlpy PROPERTIES SUFFIX ".pyd")
endif()


install(TARGETS srwlpy LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})