find_package(FFTW2 QUIET COMPONENTS FLOAT_LIB)

if(FFTW2_FOUND)
  message(STATUS "Found FFTW2F: ${FFTW_FLOAT_LIB}")

  set(
    FFTW_DOUBLE_LIB ${FFTW_DOUBLE_LIB}
    CACHE PATH "Path to FFTW"
    FORCE
  )

else()
  message(STATUS "Suitable FFTW could not be located. Downloading and building!")

  include(ExternalProject)
  ExternalProject_Add(fftw_external
    URL
      http://www.fftw.org/fftw-2.1.5.tar.gz
    URL_HASH
      MD5=8d16a84f3ca02a785ef9eb36249ba433
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    # This cp line is important. FFTW 2 must be build inplace otherwise
    # the make part will generate invalid library...
    COMMAND cp -r ${CMAKE_CURRENT_BINARY_DIR}/src/fftw_external/. .
    COMMAND ./configure --enable-float --with-pic
    COMMAND sed -i.bak -e "s/^CFLAGS = /CFLAGS = -fPIC /" Makefile
    COMMAND make -j8
    COMMAND mkdir -p ${STAGED_INSTALL_PREFIX}/lib
    COMMAND cp fftw/.libs/libfftw.a ${STAGED_INSTALL_PREFIX}/lib
    INSTALL_COMMAND ""
    )

  include(GNUInstallDirs)

  set(
    FFTW_ROOT ${STAGED_INSTALL_PREFIX}
    CACHE PATH "Path to internally built FFTWConfig.cmake"
    FORCE
  )

  set(
    FFTW_DOUBLE_LIB ${STAGED_INSTALL_PREFIX}/lib/libfftw.a
    CACHE PATH "Path to FFTW"
    FORCE
  )

  set(
    FFTW_FLOAT_LIB
    CACHE PATH "Path to FFTWF"
    FORCE
  )


  # Libraries
  add_library(fftw STATIC IMPORTED)
  set_property(TARGET fftw PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/libfftw.a)

  add_dependencies(fftw fftw_external)

endif()
