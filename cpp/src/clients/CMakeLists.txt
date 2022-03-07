if(BUILD_CLIENT_C)
    add_subdirectory(c)
endif()

if(BUILD_CLIENT_PYTHON)
    add_subdirectory(python)
endif()

if(BUILD_CLIENT_IGOR)
    add_subdirectory(igor)
endif()