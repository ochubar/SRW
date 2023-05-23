find_package(FFTW QUIET COMPONENTS FLOAT_LIB DOUBLE_LIB)

if(FFTW_FOUND)
  message(STATUS "Found FFTW3: ${FFTW_DOUBLE_LIB}")
  message(STATUS "Found FFTW3F: ${FFTW_FLOAT_LIB}")

  set(
    FFTW_DOUBLE_LIB ${FFTW_DOUBLE_LIB}
    CACHE PATH "Path to FFTW"
    FORCE
  )

  set(
    FFTW_FLOAT_LIB ${FFTW_FLOAT_LIB}
    CACHE PATH "Path to FFTWF"
    FORCE
  )

else()
  message(STATUS "Suitable FFTW3 could not be located. Downloading and building!")
  include(ExternalProject)

  if(UNIX)
    ExternalProject_Add(fftw3_external
    URL
        http://www.fftw.org/fftw-3.3.8.tar.gz
    URL_HASH
        MD5=8aac833c943d8e90d51b697b27d4384d
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    DOWNLOAD_DIR $(CMAKE_CURRENT_LIST_DIR)
    CONFIGURE_COMMAND ""
    BUILD_COMMAND
        ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/configure --enable-float --with-pic --prefix=${STAGED_INSTALL_PREFIX} &&
        make -j8 &&
        make install &&
        make clean &&
        ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/configure --with-pic --prefix=${STAGED_INSTALL_PREFIX} &&
        make -j8 &&
        make install
    INSTALL_COMMAND ""
    PREFIX=${CMAKE_CURRENT_BINARY_DIR}
    )

    include(GNUInstallDirs)

    set(
    FFTW_ROOT ${STAGED_INSTALL_PREFIX}
    CACHE PATH "Path to internally built FFTW3Config.cmake"
    FORCE
    )

    set(
    FFTW_DOUBLE_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3.a
    CACHE PATH "Path to FFTW"
    FORCE
    )

    set(
    FFTW_FLOAT_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw3f.a
    CACHE PATH "Path to FFTWF"
    FORCE
    )

    # Libraries
    add_library(fftw3 STATIC IMPORTED)
    set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/libfftw3.a)
    add_library(fftw3f STATIC IMPORTED)
    set_property(TARGET fftw3f PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/libfftw3f.a)
  endif()
  if(WIN32)
    #ExternalProject_Add(fftw3_external
    #  URL
    #    http://www.fftw.org/fftw-3.3.8.tar.gz
    #  URL_HASH
    #    MD5=8aac833c943d8e90d51b697b27d4384d
    #  PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    #  DOWNLOAD_DIR $(CMAKE_CURRENT_LIST_DIR)
    #  CONFIGURE_COMMAND ""
    #  BUILD_COMMAND ""
    #  COMMAND cmake ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/ -DBUILD_SHARED_LIBS=OFF -DENABLE_FLOAT=ON -DBUILD_TESTS=OFF -DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
    #  COMMAND cmake --build . --config Release --target install
    #  COMMAND cmake ${CMAKE_CURRENT_BINARY_DIR}/src/fftw3_external/ -DBUILD_SHARED_LIBS=OFF -DENABLE_FLOAT=OFF -DBUILD_TESTS=OFF-DCMAKE_INSTALL_PREFIX=${STAGED_INSTALL_PREFIX}
    #  COMMAND cmake --build . --config Release --target install
    #  INSTALL_COMMAND ""
    #  )
    #
    #include(GNUInstallDirs)
    #
    #
    #set(
    #  FFTW_ROOT ${STAGED_INSTALL_PREFIX}
    #  CACHE PATH "Path to internally built FFTW3Config.cmake"
    #  FORCE
    #)
    #
    #set(
    #  FFTW_DOUBLE_LIB ${STAGED_INSTALL_PREFIX}/lib/fftw3.lib
    #  CACHE PATH "Path to FFTW"
    #  FORCE
    #)
    #
    #set(
    #  FFTW_FLOAT_LIB ${STAGED_INSTALL_PREFIX}/lib/fftw3f.lib
    #  CACHE PATH "Path to FFTWF"
    #  FORCE
    #)
    #
    ## Libraries
    #add_library(fftw3 STATIC IMPORTED)
    #set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/fftw3.lib)
    #add_library(fftw3f STATIC IMPORTED)
    #set_property(TARGET fftw3f PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/fftw3f.lib)

    set(
    FFTW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/..
    FORCE
    )

    set(
    FFTW_DOUBLE_LIB ${CMAKE_CURRENT_SOURCE_DIR}/../fftw3_64.lib
    CACHE PATH "Path to FFTW"
    FORCE
    )

    set(
    FFTW_FLOAT_LIB ${CMAKE_CURRENT_SOURCE_DIR}/../fftw3f_64.lib
    CACHE PATH "Path to FFTWF"
    FORCE
    )

    # Libraries
    add_library(fftw3 STATIC IMPORTED)
    set_property(TARGET fftw3 PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../fftw3_64.lib)
    add_library(fftw3f STATIC IMPORTED)
    set_property(TARGET fftw3f PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../fftw3f_64.lib)
  endif()

  add_dependencies(fftw3 fftw3_external)
  add_dependencies(fftw3f fftw3_external)

endif()
