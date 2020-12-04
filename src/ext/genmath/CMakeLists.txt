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

add_library(genmath OBJECT ${genmath_source_files})
target_link_libraries(genmath PUBLIC auxparse)
target_include_directories(genmath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(USE_OPENMP)
    target_include_directories(genmath PUBLIC ${OpenMP_CXX_INCLUDE_DIRS})
endif()