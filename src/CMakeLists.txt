add_subdirectory(ext)
add_subdirectory(core)
add_subdirectory(lib)

if(BUILD_CLIENTS)
    add_subdirectory(clients)
endif()