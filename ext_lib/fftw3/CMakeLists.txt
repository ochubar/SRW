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
  ExternalProject_Add(fftw3_external
    URL
      http://www.fftw.org/fftw-3.3.8.tar.gz
    URL_HASH
      MD5=8aac833c943d8e90d51b697b27d4384d
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
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

  add_dependencies(fftw3 fftw3_external)
  add_dependencies(fftw3f fftw3_external)

endif()
