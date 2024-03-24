# This Makefile is to compile SRW with optional compilation of FFTW libraries
#
# The following options are available:
# - `make all` - will compile FFTW, C++ core and Python lib;
# - `make fftw` - will compile FFTW only;
# - `make` - will compile C++ core and Python lib;
# - `make clean` - will clean temporary files.
#
# Updated by Maksim Rakitin (NSLS-II, BNL) on May 2, 2016.

root_dir = $(realpath .)
ext_dir = $(root_dir)/ext_lib
gcc_dir = $(root_dir)/cpp/gcc
py_dir = $(root_dir)/cpp/py
fftw2_version = fftw-2.1.5
fftw2_dir = $(fftw2_version)
fftw2_file = fftw/$(fftw2_version).tar.gz
fftw3_version = fftw-3.3.8
fftw3_dir = $(fftw3_version)
fftw3_file = fftw3/$(fftw3_version).tar.gz
log_fftw = /dev/null
examples_dir = env/python/srwpy/examples
#example10_data_dir = $(examples_dir)/data_example_10
export MODE ?= 0

#HG20042024
ifeq ($(MODE), cuda)
export CUDA_PATH ?= /usr/local/cuda
export CUDA_MATHLIBS_PATH ?= /usr/local/cuda
endif

nofftw: core pylib

all: clean fftw core pylib

fftw:
	if [ ! -d "$(ext_dir)" ]; then \
	    mkdir $(ext_dir); \
	fi; \
	cd $(ext_dir); \
#	if [ ! -f "$(fftw2_file)" ]; then \
#	    wget https://raw.githubusercontent.com/ochubar/SRW/master/ext_lib/$(fftw2_file); \
#	fi; \
	if [ -d "$(fftw2_dir)" ]; then \
	    rm -rf $(fftw2_dir); \
	fi; \
	tar -zxf $(fftw2_file); \
	cd $(fftw2_dir); \
	./configure --enable-float --with-pic; \
	sed 's/^CFLAGS = /CFLAGS = -fPIC /' -i Makefile; \
	make -j8 && cp fftw/.libs/libfftw.a $(ext_dir); \
	cd $(ext_dir); \
	rm -rf $(fftw2_dir); \
	if [ -d "$(fftw3_dir)" ]; then \
	    rm -rf $(fftw3_dir); \
	fi; \
	tar -zxf $(fftw3_file); \
	cd $(fftw3_dir); \
	./configure --enable-float --with-pic; \
	sed 's/^CFLAGS = /CFLAGS = -fPIC /' -i Makefile; \
	make -j8 && cp .libs/libfftw3f.a ../; \
	cd $(ext_dir); \
	rm -rf $(fftw3_dir); \
	tar -zxf $(fftw3_file); \
	cd $(fftw3_dir); \
	./configure --with-pic; \
	sed 's/^CFLAGS = /CFLAGS = -fPIC /' -i Makefile; \
	make -j8 && cp .libs/libfftw3.a ../; \
	cd $(root_dir); \
	rm -rf $(ext_dir)/$(fftw3_dir);

core: 
	cd $(gcc_dir); make -j8 clean lib

pylib:
	cd $(py_dir); make python

clean:
	rm -f $(ext_dir)/libfftw3f.a $(ext_dir)/libfftw3.a $(gcc_dir)/libsrw.a $(gcc_dir)/srwlpy*.so; \
	rm -rf $(ext_dir)/$(fftw2_dir)/ $(ext_dir)/$(fftw3_dir)/ py/build/;
	if [ -d $(root_dir)/.git ]; then rm -f $(examples_dir)/srwlpy*.so && (git checkout $(examples_dir)/srwlpy*.so 2>/dev/null || :); fi;

.PHONY: all clean core fftw nofftw pylib

