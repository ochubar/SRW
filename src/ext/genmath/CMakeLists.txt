set(genmath_source_files
    gmfft.cpp
    gmfft.h
    gmfit.cpp
    gmfit.h
    gminterp.cpp
    gminterp.h
    gmmeth.cpp
    gmmeth.h
    gmtrans.cpp
    gmtrans.h
)

add_library(genmath ${genmath_source_files})

if (USE_FFTW3)
    target_link_libraries(genmath auxparse fftw fftwf)
else()
    target_link_libraries(genmath auxparse fftw)
endif()


target_include_directories(genmath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})