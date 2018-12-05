# This Makefile is to compile SRW with optional compilation of FFTW-2.1.5.
#
# The following options are available:
# - `make all` - will compile FFTW, C++ core and Python lib;
# - `make fftw` - will compile FFTW only;
# - `make` - will compile C++ core and Python lib;
# - `make clean` - will clean temporary files.
#
# Updated by Maksim Rakitin (NSLS-II, BNL) on May 2, 2016.

root_dir = $(realpath .)
env_dir = $(root_dir)/env
ext_dir = $(root_dir)/ext_lib
gcc_dir = $(root_dir)/cpp/gcc
py_dir = $(root_dir)/cpp/py
fftw_version = fftw-2.1.5
fftw_dir = $(fftw_version)
fftw_file = $(fftw_version).tar.gz
log_fftw = /dev/null
examples_dir = $(env_dir)/work/srw_python
#example10_data_dir = $(examples_dir)/data_example_10
export MODE ?= 0

nofftw: core pylib

all: clean fftw core pylib

fftw:
	if [ ! -d "$(ext_dir)" ]; then \
	    mkdir $(ext_dir); \
	fi; \
	cd $(ext_dir); \
	if [ ! -f "$(fftw_file)" ]; then \
	    wget https://raw.githubusercontent.com/ochubar/SRW/master/ext_lib/$(fftw_file); \
	fi; \
	if [ -d "$(fftw_dir)" ]; then \
	    rm -rf $(fftw_dir); \
	fi; \
	tar -zxf $(fftw_file); \
	cd $(fftw_dir); \
	./configure --enable-float --with-pic; \
	sed 's/^CFLAGS = /CFLAGS = -fPIC /' -i Makefile; \
	make -j8 && cp fftw/.libs/libfftw.a $(ext_dir); \
	cd $(root_dir); \
	rm -rf $(ext_dir)/$(fftw_dir);

core: 
	cd $(gcc_dir); make -j8 clean lib

pylib:
	cd $(py_dir); make python

clean:
	rm -f $(ext_dir)/libfftw.a $(gcc_dir)/libsrw.a $(gcc_dir)/srwlpy*.so; \
	rm -rf $(ext_dir)/$(fftw_dir)/ py/build/;
	if [ -d $(root_dir)/.git ]; then rm -f $(examples_dir)/srwlpy*.so && (git checkout $(examples_dir)/srwlpy*.so 2>/dev/null || :); fi;

.PHONY: all clean core fftw nofftw pylib
