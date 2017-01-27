# This Makefile is to compile SRW with optional compilation of FFTW-2.1.5.
#
# The following options are available:
# - `make all` - will compile FFTW, C++ core and Python lib;
# - `make fftw` - will compile FFTW only;
# - `make` - will compile C++ core and Python lib;
# - `make test` - will execute `python SRWLIB_Example10.py` during 20 seconds;
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
example10_data_dir = $(examples_dir)/data_example_10

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

test:
	remove_tmp_dir=0; \
	if [ ! -d "$(example10_data_dir)" ]; then \
	    mkdir $(example10_data_dir); \
	    remove_tmp_dir=1; \
	fi; \
	cd $(examples_dir); \
	timeout 20 python SRWLIB_Example10.py; \
	code=$$?; \
	RED=1; \
	GREEN=2; \
	if [ $$code -eq 0 ]; then \
	    status='PASSED'; \
	    color=$${GREEN}; \
	    message=''; \
	elif [ $$code -eq 124 ]; then \
	    status='PASSED'; \
	    color=$${GREEN}; \
	    message=' (timeouted, expected)'; \
	else \
	    status='FAILED'; \
	    color=$${RED}; \
	    message=''; \
	fi; \
	echo -e -n "\n\tTest "; \
	tput setaf $${color}; \
	tput bold; \
	echo -e -n "$${status}"; \
	tput sgr0; \
	echo -e ". Code=$${code}$${message}\n"; \
	rm -f $(example10_data_dir)/{ex10_res_int_se.dat,ex10_res_int_prop_se.dat,ex10_res_int_prop_me.dat}; \
	if [ $$remove_tmp_dir -eq 1 ]; then \
	    cd $(root_dir); \
	    rm -rf $(example10_data_dir); \
	fi;

clean:
	rm -f $(ext_dir)/libfftw.a $(gcc_dir)/libsrw.a $(gcc_dir)/srwlpy*.so; \
	rm -rf $(ext_dir)/$(fftw_dir)/ py/build/;
	if [ -d $(root_dir)/.git ]; then rm -f $(examples_dir)/srwlpy*.so && (git checkout $(examples_dir)/srwlpy*.so 2>/dev/null || :); fi;

.PHONY: all clean core fftw nofftw pylib test
