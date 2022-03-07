set(srw_clients_c_source_files
    srwlclient.cpp
)

add_executable(srwlclient ${srw_clients_c_source_files})

target_link_libraries(srwlclient srw)

target_include_directories(srwlclient PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
)

install(TARGETS srwlclient RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})