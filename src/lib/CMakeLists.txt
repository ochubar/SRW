set(srw_lib_source_files
    srwlib.cpp
    srwlib.h
    srerror.cpp
    srerror.h
)

add_library(srw STATIC
    ${srw_lib_source_files}
    $<TARGET_OBJECTS:core>
    $<TARGET_OBJECTS:genesis>
    $<TARGET_OBJECTS:genmath>
    $<TARGET_OBJECTS:auxparse>
)

find_package(FFTW REQUIRED COMPONENTS FLOAT_LIB DOUBLE_LIB)
message(STATUS "Found FFTW3: ${FFTW_DOUBLE_LIB}")
message(STATUS "Found FFTW3F: ${FFTW_FLOAT_LIB}")


target_link_libraries(srw PUBLIC core m FFTW::Float FFTW::Double)

set(LINK_FLAGS ${LINK_FLAGS} "-Wl,-whole-archive")


target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../core>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../ext/auxparse>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../ext/genesis>
)
target_include_directories(srw PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../ext/genmath>
)

install(TARGETS srw LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})