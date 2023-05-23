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

    if(UNIX)
        # https://cmake.org/cmake/help/latest/module/ExternalProject.html
        ExternalProject_Add(fftw_external
            URL
              http://www.fftw.org/fftw-2.1.5.tar.gz
            URL_HASH
              MD5=8d16a84f3ca02a785ef9eb36249ba433
            CMAKE_ARGS "${CMAKE_ARGS};-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}"
            PREFIX ${CMAKE_CURRENT_BINARY_DIR}
            DOWNLOAD_DIR $(CMAKE_CURRENT_LIST_DIR)
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            # This cp line is important. FFTW 2 must be build inplace otherwise
            # the make part will generate invalid library...
            COMMAND cp -r ${CMAKE_CURRENT_BINARY_DIR}/src/fftw_external/. .
            # Update config.sub and config.guess so fftw can be compiled with new target osx-arm64
            # Note (2021-10-22): the git.savannah.gnu.org website is not
            # accessible intermittently, so vendoring those two files:
            # COMMAND wget -O config.sub https://git.savannah.gnu.org/gitweb/?p=config.git&a=blob_plain&f=config.sub&hb=HEAD
            # COMMAND wget -O config.guess https://git.savannah.gnu.org/gitweb/?p=config.git&a=blob_plain&f=config.guess&hb=HEAD
            COMMAND cp -v ${CMAKE_CURRENT_SOURCE_DIR}/config.sub .
            COMMAND cp -v ${CMAKE_CURRENT_SOURCE_DIR}/config.guess .
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
    endif()
    if(WIN32)
        #ExternalProject_Add(fftw_external
        #    URL
        #      https://github.com/ochubar/SRW/raw/master/ext_lib/fftw64_f.lib
        #    URL_HASH
        #      MD5=3e34fe1af702ba2af176e83a39021250
        #    DOWNLOAD_NAME fftwf.lib
        #    DOWNLOAD_DIR ${STAGED_INSTALL_PREFIX}/lib
        #    DOWNLOAD_NO_EXTRACT true
        #    CONFIGURE_COMMAND ""
        #    BUILD_COMMAND ""
        #    INSTALL_COMMAND ""
        #)
        #
        #include(GNUInstallDirs)
        #
        #set(
        #    FFTW_ROOT ${STAGED_INSTALL_PREFIX}
        #    CACHE PATH "Path to internally built FFTWConfig.cmake"
        #    FORCE
        #)
        #
        #set(
        #    FFTW_DOUBLE_LIB ${STAGED_INSTALL_PREFIX}/lib/fftwf.lib
        #    CACHE PATH "Path to FFTW"
        #    FORCE
        #)
        #
        #set(
        #    FFTW_FLOAT_LIB
        #    CACHE PATH "Path to FFTWF"
        #    FORCE
        #)
        #
        # Libraries
        #add_library(fftw STATIC IMPORTED)
        #set_property(TARGET fftw PROPERTY IMPORTED_LOCATION ${STAGED_INSTALL_PREFIX}/lib/fftwf.lib)

        set(
            FFTW_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/ext_lib
            FORCE
        )

        if(CMAKE_SIZEOF_VOID_P EQUAL 8)
            set(
                FFTW_DOUBLE_LIB
                CACHE PATH "Path to FFTW"
                FORCE
            )

            set(
                FFTW_FLOAT_LIB ${CMAKE_CURRENT_SOURCE_DIR}/../fftw64_f.lib
                CACHE PATH "Path to FFTWF"
                FORCE
            )
        elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
            set(
                FFTW_DOUBLE_LIB
                CACHE PATH "Path to FFTW"
                FORCE
            )

            set(
                FFTW_FLOAT_LIB ${CMAKE_CURRENT_SOURCE_DIR}/../fftw_f.lib
                CACHE PATH "Path to FFTWF"
                FORCE
            )
        endif()


        # Libraries
        add_library(fftw STATIC IMPORTED)
        set_property(TARGET fftw PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../fftw_f.lib)
    endif()
    add_dependencies(fftw fftw_external)

endif()
