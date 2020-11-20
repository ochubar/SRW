set(srw_lib_source_files
    srwlib.cpp
    srwlib.h
)

add_library(srwlib ${srw_lib_source_files})

if (USE_FFTW3)
    target_link_libraries(srwlib auxparse genmath genesis core fftw fftwf)
else()
    find_package(FFTW CONFIG REQUIRED)
    get_property(_loc TARGET FFTW::fftw PROPERTY LOCATION)
    message(STATUS "Found FFTW: ${_loc} (found version ${FFTW_VERSION})")
    target_link_libraries(srwlib auxparse genmath genesis core FFTW::fftw)
endif()

target_include_directories(srwlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
)
target_include_directories(srwlib PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../core>
)