Comments to "SRWLIB" for Python (November 2017)
====================================================

Basic content and installation of the SRWLIB package:
----------------------------------------------------

- srwlpy.pyd, srwlpy.so are SRW Python bindings (shared libraries) compiled for Windows and Linux respectively. Note that there is no guarantee that these files, located in ../srw_python/, are compatible with the versions of operating system and Python that you are using.
The compiled shared library files for Python 3.x and 2.7 for 64- and 32-bit Windows and Linux are available in ../srw_python/lib/. The file names of the compiled shared libraries are self-explanatory: e.g. srwlpy3_6_x64.pyd is a version for Python 3.6 for Windows 64 bit; srwlpy2_7_x86_64.so is for Python 2.7 for Linux 64 bit. Copy the compiled shared library file corresponding to your operating system and Python versions from ../srw_python/lib/ to ../srw_python/ and change its name to srwlpy.pyd (on Windows) or srwlpy.so (on Linux) e.g.:
on Windows:
copy ..\srw_python\lib\srwlpy3_6_x64.pyd ..\srw_python\srwlpy.pyd
on Linux:
cp ../srw_python/lib/srwlpy2_7_x86_64.so ../srw_python/srwlpy.so
where ".." is the path to the srw_python folder.
 
- srwlib.py is SRWLIB Python module definition file; it is compatible both with Python 3.2 (or higher) and 2.7 (or higher).

- SRWLIB_Example*.py files are SRWLIB Python example scripts, which are also compatible with Python 3.2 and 2.7.

- ../srw_python/lib/srw*.lib, ../srw_python/lib/libsrw*.a are static SRW libraries (C API) compiled for 32- (srw_win32.lib) and 64-bit (srw_x64.lib) Windows and 64-bit Linux (libsrw_x86_64.a). These files may not be necessary for the SRWLIB for Python (because *.pyd and *.so files should normally contain these static libraries linked).

- ../srw_python/lib/srwlib.h is the header file of SRW C API (it defines C structures and declares functions of SRWLIB).


To run SRWLIB examples:
----------------------------------------------------

cd ../srw_python
python SRWLIB_ExampleXX.py

The examples may read some input (ASCII) files and save some output (also mostly ASCII) files to the corresponding sub-folders located in ../srw_python/, e.g. ../srw_python/data_example_01/ for the Example # 1; see the example script files for details.


Optional configuring of Python and SRWLIB on Linux:
----------------------------------------------------

Make sure that path to Python 3.6 (or 2.7) is added to the PATH variable and "srw_python" directory to PYTHONPATH variable:
export PATH="$PATH:<absolute path to Python 3.6>" #this is not necessary if you install python using the distro's package manager
export PYTHONPATH="$PYTHONPATH:SRW_Dev/work/srw_python/" #temporarely solution
	or:
echo "export PYTHONPATH=$PYTHONPATH:SRW_Dev/work/srw_python/" >> ~/.bashrc # permanent solution for a single user
Setting up PYTHONPATH allows to import srwlpy module from any directory.

