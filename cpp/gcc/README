Compiling and testing SRW Library and its Python binding on Linux:

1. Download and compile fftw-2.1.5 library as required for SRW (assuming that absolute path to SRW directory is "SRW_Dev"):
	Download fftw-2.1.5.tar.gz from FFTW site (probably http://www.fftw.org/download.html) and place it to SRW_Dev/ext_lib
	cd SRW_Dev/ext_lib
	tar -zxvf fftw-2.1.5.tar.gz
	cd fftw-2.1.5
	configure --enable-float --with-pic
	Manually (using editor) add -fPIC option to CFLAGS in Makefile
	make
	make install
	cp fftw/.libs/libfftw.a SRW_Dev/ext_lib/

2. Compile SRWLib library and Python binding:
	cd SRW_Dev/cpp/gcc
	Make sure that Python 3.2 or higher (or Python 2.7) is installed 
	In Makefile, specify correct PYFLAGS path to Python include file and library
	rm libsrw.a
	make all
	cp srwlpy.so ../../env/work/srw_python/

3. Check examples:
	Make sure that path to Python 3.2 (or 2.7) is added to the PATH variable and "srw_python" to PYTHONPATH variable:
	export PATH="$PATH:<absolute path to Python 3.2>" # this is not necessary if you install python using the distro's package manager
	export PYTHONPAH="$PYTHONPATH:SRW_Dev/env/work/srw_python/" # temporarely solution
	or
	echo "export PYTHONPATH=$PYTHONPATH:SRW_Dev/env/work/srw_python/" >> ~/.bashrc # permanent solution for a single user
	Setting up PYTHONPATH allows to import srwlpy module from any directory.
	cd SRW_Dev/env/work/srw_python
	python SRWLIB_ExampleXX.py