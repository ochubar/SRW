set(auxparse_source_files
    auxparse.cpp
    auxparse.h
)

add_library(auxparse OBJECT ${auxparse_source_files})
target_include_directories(auxparse PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
