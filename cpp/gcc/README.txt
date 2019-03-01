I. Compiling and testing SRW Library and its Python binding on Linux (without using Python "distutils")
------------------------------------------------------------------
I.1. Download and compile fftw-2.1.5 library as required for SRW.
	Download fftw-2.1.5.tar.gz from FFTW site (probably http://www.fftw.org/download.html) and place it to SRW_Dev/ext_lib:
	cd SRW_Dev/ext_lib
	tar -zxvf fftw-2.1.5.tar.gz
	cd fftw-2.1.5
	./configure --enable-float --with-pic
	Manually (using editor) add -fPIC option to CFLAGS in Makefile
	make
	make install
	cp fftw/.libs/libfftw.a SRW_Dev/ext_lib/

I.2. Compile the SRW library and Python binding.
	cd SRW_Dev/cpp/gcc
	Make sure Python 3.2 or higher (or Python 2.7) is installed. 
	In the SRW_Dev/cpp/gcc/Makefile, modify/correct PYPATH and PYFLAGS variables, i.e. specify path to Python header and library files. Depending on Linux environment, it may also be necessary to modify the name of compiler to be used, e.g.:
	CC  = gcc
	CXX = g++
	#CC  = cc
	#CXX = c++
	After this, execute the following:
	rm libsrw.a
	make all
	To compile SRW library supporting OpenMP based parallel calculations (e.g. for XFEL applications) use add "MODE=omp" after "make all":
	make all MODE=omp
	Then copy srwlpy.so to SRW_Dev/env/work/srw_python/:
	cp srwlpy.so ../../env/work/srw_python/

I.3. Checking the examples.
	Make sure the path to Python 3.x (or 2.7) is added to the PATH variable and "srw_python" to PYTHONPATH variable:
	export PATH="$PATH:<absolute path to Python 3.x>" # this is not necessary if you install python using the distro's package manager
	export PYTHONPAH="$PYTHONPATH:SRW_Dev/env/work/srw_python/" #temporary solution
	or
	echo "export PYTHONPATH=$PYTHONPATH:SRW_Dev/env/work/srw_python/" >> ~/.bashrc #permanent solution for a single user
	Setting up PYTHONPATH allows to import srwlpy module from any directory. Testing of the examples would preferably done in the "srw_python" directory:
	cd SRW_Dev/env/work/srw_python
	python SRWLIB_ExampleXX.py


II. Compiling and testing SRW Library and its Python binding on Mac OSX
------------------------------------------------------------------
	Try to follow the steps described in section III, describing options for compiling and testing SRW on Linux. 
	We were informed that the actions described in III.1.2.2 lead to successful compilation with gcc/g++ provided by Xcode 10.1, after the following modifications in SRW_Dev/cpp/gcc/Makefile:
	CC  = gcc
	CXX = g++
	#CC  = cc
	#CXX = c++
	...
	PYPATH=/Library/Frameworks/Python.framework/Versions/3.6
	PYFLAGS=-I$(PYPATH)/include/python3.6m -I$(PYPATH)/include/python3.6m -L$(PYPATH)/lib/python3.6/config-3.6m-darwin -lpython3.6m -ldl
	
	The correct path and flags can be obtained e.g. by executing from command line:
	python3-config --includes --ldflags
	and removing the option -framework

	With earlier versions of Xcode, the following manipulations, consisting in installation of "macports" and obtaining the whole gcc toolchain, were reported to be successful:
	sudo port install gcc47
	Modify the SRW_Dev/cpp/gcc/Makefile so that CC=<path to macports>/gcc and CXX=<path to macports>/g++, and proceed to the compilation as described in I.2.	