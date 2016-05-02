root_dir = $(realpath .)
ext_dir = $(root_dir)/ext_lib
gcc_dir = $(root_dir)/cpp/gcc
py_dir = $(root_dir)/cpp/py
fftw_version = fftw-2.1.5
fftw_dir = $(fftw_version)
fftw_file = $(fftw_version).tar.gz
log_fftw = /dev/null

nofftw: core pylib

all: clean fftw core pylib

fftw:
	cd $(ext_dir); \
	if [ ! -f "$(fftw_file)" ]; then \
	    wget https://raw.githubusercontent.com/ochubar/SRW/master/ext_lib/$(fftw_file) > $(log_fftw) 2>&1; \
	fi; \
	if [ -d "$(fftw_dir)" ]; then \
	    rm -rf $(fftw_dir); \
	fi; \
	tar -zxf $(fftw_file); \
	cd $(fftw_dir); \
	./configure --enable-float --with-pic --prefix=$(ext_dir)/$(fftw_dir)/tmp >> $(log_fftw) 2>&1; \
	sed 's/^CFLAGS = /CFLAGS = -fPIC /' -i Makefile; \
	make && cp fftw/.libs/libfftw.a $(ext_dir);

core: 
	cd $(gcc_dir); make -j8 clean lib

pylib:
	cd $(py_dir); make python

clean:
	rm -f $(ext_dir)/libfftw.a $(gcc_dir)/libsrw.a $(gcc_dir)/srwlpy.so \
	rm -rf $(ext_dir)/$(fftw_dir)/ py/build/

.PHONY: all core clean fftw nofftw pylib
