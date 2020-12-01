set(srw_clients_python_source_files
    srwlpy.cpp
)

add_library(srwlpy SHARED ${srw_clients_python_source_files})

# for testing we will need the python interpreter
find_package(PythonInterp REQUIRED)

# we require python development headers
find_package(PythonLibs ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR} EXACT REQUIRED)

target_link_libraries(srwlpy srw)

# We need the undefined dynamic_lookup to overcome issues with a statically linked
# python. More details here: https://gitlab.kitware.com/vtk/vtk/-/issues/17214
set_target_properties(srwlpy PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")

target_include_directories(srwlpy PRIVATE
    ${PYTHON_INCLUDE_DIRS}
)

target_include_directories(srwlpy PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
)

# prevent cmake from creating a "lib" prefix
set_target_properties(srwlpy
  PROPERTIES
    PREFIX ""
)

if(WIN32)
  # python will not import dll but expects pyd
  set_target_properties(srwlpy
    PROPERTIES
      SUFFIX ".pyd"
    )
endif()
if(UNIX)
  # python will not import dll but expects pyd
  set_target_properties(srwlpy
    PROPERTIES
      SUFFIX ".so"
    )
endif()

install(TARGETS srwlpy LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})