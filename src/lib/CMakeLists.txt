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


# Add definitions needed to build SRW
target_compile_definitions(srw PUBLIC ${SRW_DEFINITIONS})

# Libraries in which SRW lib depends
target_link_libraries(srw core ${FFTW_DOUBLE_LIB} ${FFTW_FLOAT_LIB})

if (UNIX)
    # Math library is only needed by unix. When using Windows, VS does it for you.
    target_link_libraries(srw m)
endif()

if (USE_OPENMP)
    target_link_libraries(srw OpenMP::OpenMP_CXX)
    target_include_directories(srw PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
endif()

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