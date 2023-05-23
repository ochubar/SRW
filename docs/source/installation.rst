============
Installation
============

SRW uses leverages the OS portability features of CMake for its build system.
This package consists of a C++ library and client libraries for:

- C
- Python
- Igor PRO

------------
Dependencies
------------

SRW depends on the following external dependencies:

- CMake >= 3.12
- C/C++ compiler
- FFTW library (3.x for regular build or 2.x if selecting OpenMP)

Note: The proper FFTW library will be built as part of the C++ library if not found on the system.

Optional Dependencies
=====================

- Python 2.7 or 3.6+ (required for Python client)
- Igor Pro (required for Igor Pro client)
- Igor XOP (required for Igor Pro client)

---------------
Build Procedure
---------------

As mentioned above, SRW relies on CMake for generation of the files necessary for
building the main SRW and client libraries.

For users looking for the Python client, the procedure is very straightforward
and does not require source build. Please find the instructions below.

The Source Build is only recommended for users that need the main SRW library to
link with other code or the C and Igor client libraries.

Source Build
============

Before diving on the details for the source build of SRW library it is important
to first look into the options available via CMake to configure the build and install
of the package.

CMake Options
^^^^^^^^^^^^^

Here are the cmake options along with their usage and description for SRW library.

======================  ====== =============================================== ==========================================================
Option                  Type   Usage                                           Description
======================  ====== =============================================== ==========================================================
-DUSE_OPENMP            (bool) -DUSE_OPENMP=<ON|OFF>                           Whether or not to build OpenMP capable library
-DBUILD_CLIENT_C 		(bool) -DBUILD_CLIENT_C=<ON|OFF>                       Whether or not to build the C client library
-DBUILD_CLIENT_PYTHON   (bool) -DBUILD_CLIENT_PYTHON=<ON|OFF>                  Whether or not to build the Python client library
-DCMAKE_INSTALL_PREFIX  (str)  -DCMAKE_INSTALL_PREFIX=/home/user/software/SRW  The path in which to install the build artifacts
======================  ====== =============================================== ==========================================================


Linux/macOS Instructions
^^^^^^^^^^^^^^^^^^^^^^^^

When building cmake-based projects, it is common practice to do so in a folder other than the project root directory.

As a way to demonstrate how to build and install the SRW library we will simulate a case in which an user wants to
build and install the SRW library (with OpenMP) along with the C client only and have it installed at a specific folder under its home directory.

Here are the steps::

  $ cd SRW  # This is the folder in which you have the SRW code
  $ mkdir build
  $ cd build
  $ cmake .. -DUSE_OPENMP=ON -DBUILD_CLIENT_C=ON -DCMAKE_INSTALL_PREFIX=/home/user/my_install_folder/
  $ make
  $ make install

Windows
^^^^^^^

When building cmake-based projects, it is common practice to do so in a folder other than the project root directory.

As a way to demonstrate how to build and install the SRW library we will simulate a case in which an user wants to
build and install the SRW library (with OpenMP) along with the C client only and have it installed at a specific folder under its C: drive.

Here are the steps::

  $ cd SRW  # This is the folder in which you have the SRW code
  $ mkdir build
  $ cd build
  $ cmake .. -GNMake Makefiles -DUSE_OPENMP=ON -DBUILD_CLIENT_C=ON -DCMAKE_INSTALL_PREFIX=C:\my_install_folder
  $ cmake --build . --config Release --target install


Note: For Windows it is necessary to specify the *-G* option to select a build system that does not support multi-configuration builds otherwise it will install
the artifacts in a folder called Release or Debug depending on the *--config* selected.

Python Client
=============

When installing the Python client, it will automatically take care of building
the SRW library as part of the process.

Linux/macOS Instructions
^^^^^^^^^^^^^^^^^^^^^^^^

In order to install the regular SRW (non-OpenMP) do::

  $ cd python
  $ python setup.py install

In order to install the OpenMP capable SRW do::

  $ cd python
  $ MODE=omp python setup.py install

Windows
^^^^^^^

In order to install the regular SRW (non-OpenMP) do::

  $ cd python
  $ python setup.py install

In order to install the OpenMP capable SRW do::

  $ cd python
  $ set MODE=omp
  $ python setup.py install


Igor Pro Client
===============

TODO